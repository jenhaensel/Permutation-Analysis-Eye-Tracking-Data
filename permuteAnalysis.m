% Part of analysis script used for Haensel et al. (2020). Culture modulates face scanning 
% during dyadic interactions. See publication for data.

% To run the script, the following toolboxes are additionally required:
% - cosmo_mvpa (http://cosmomvpa.org/download.html)
% - fieldtrip (http://www.fieldtriptoolbox.org/download)

clear; clc; close all

%% Define participants and constants.
% Define Japanese participants
subjJP = [2;4;5;6;8;9;10;11;13;15;16;17;18;19;21;22;23;24;25;26;27;29;31;32;33;34;35];
subjNo_JP = sprintfc('A%02d_JP', subjJP);

% Define British/Irish participants
subjBR = [1;2;3;6;7;9;10;12;13;14;15;16;17;18;19;20;21;22;24;25;26;28;29;30;31;33;34;35;36];
subjNo_BR = sprintfc('A%02d_UK', subjBR);

% Set uncorrected p-value
uncorrected_pval = 0.01;


%% Load and summarise data
no_subjJP = length(subjNo_JP);
no_subjBR = length(subjNo_BR);

% Summarise data of Japanese participants
for N = 1:length(subjNo_JP)
    % Loads data from one participant in listening condition
    load(['mapListen_' num2str(subjNo_JP{N,1}) '.mat']);
    smoothGazeListen = smoothGazeFIX;
    clear smoothGazeFIX
    
    if N==1
        imgJP_ListenALL = zeros(size(smoothGazeListen,1), size(smoothGazeListen,2), no_subjJP);
    end
    % The following matrices contain all images for each participant
    imgJP_ListenALL(:,:,N) = smoothGazeListen;    
    
    subjID_JP(N,:) = N;
    
    clearvars -except subjNo_JP subjNo_BR no_subjJP no_subjBR N imgJP_ListenALL subjID_JP uncorrected_pval subjID_BR...
        imgBR_ListenALL
end


% Summarise data of British/Irish participants
for N = 1:length(subjNo_BR)
    % Loads data from one participant in listening condition
    load(['mapListen_' num2str(subjNo_BR{N,1}) '.mat']); 
    smoothGazeListen = smoothGazeFIX;
    clear smoothGazeFIX
    
    if N==1
        imgBR_ListenALL = zeros(size(smoothGazeListen,1), size(smoothGazeListen,2), no_subjBR);
    end
    % The following matrices contain all images for each participant
    imgBR_ListenALL(:,:,N) = smoothGazeListen; 
    
    subjID_BR(N,:) = N + subjID_JP(end);
    
    clearvars -except subjNo_JP subjNo_BR no_subjJP no_subjBR N imgJP_ListenALL subjID_JP ...
        imgBR_ListenALL subjID_BR uncorrected_pval 
end


%% Permutation analysis (CoSMoMVPA)
ft_defaults;

% Rearrange the dimensions for Cosmo
imgJP_Listen_perm = permute(imgJP_ListenALL, [3 1 2]);
imgBR_Listen_perm = permute(imgBR_ListenALL, [3 1 2]);

% Listening condition: Create data for Cosmo
% Create group 1 (JP):
ds1_listen = cosmo_flatten(imgJP_Listen_perm, {'freq', 'time'}, {1:size(imgJP_Listen_perm,2), 1:size(imgJP_Listen_perm,2)});
ds1_listen.sa.chunks = subjID_JP;
ds1_listen.sa.targets = ones(no_subjJP,1);

% Create group 2 (BR):
ds2_listen = cosmo_flatten(imgBR_Listen_perm, {'freq', 'time'}, {1:size(imgBR_Listen_perm,2), 1:size(imgBR_Listen_perm,2)});
ds2_listen.sa.chunks = subjID_BR;
ds2_listen.sa.targets = ones(no_subjBR,1)*2;

% merge the two dataset 
ds_listen = cosmo_stack({ds1_listen, ds2_listen}, 1);
cosmo_check_dataset(ds_listen);

% find the neighbours of each pixel
nbr=cosmo_cluster_neighborhood(ds_listen,'progress',true);

% Set options for monte carlo based clustering statistic
opt=struct();

% choose the cluster stat (options: maxsize, maxsum, tfce)
opt.cluster_stat = 'maxsize';

% if it is not TFCE, you have to choose an initial p value. 0.005 is good
% for fmri, it might be too much or too less here
if ~strcmp(opt.cluster_stat, 'tfce')
    opt.p_uncorrected = uncorrected_pval;
end

% choose the number of permutation
opt.niter=10000;

% do the t test. The function recognizes the type of test to do from
% ds.sa.targets and ds.sa.chunks
ds_z_listen = cosmo_montecarlo_cluster_stat(ds_listen, nbr, opt);

% the function returns z values, so any pixel lower than 1.96 (two tailed
% analysis) is not significant)
ds_z_listen.samples(abs(ds_z_listen.samples)<1.96)=0;

% make the data in a matrix form again
img_Z_listen = squeeze(cosmo_unflatten(ds_z_listen));

% visulise the results
figure('Name', 'permutation analysis'); imagesc(img_Z_listen);colormap('jet');axis image;colorbar;caxis([-2 2])
