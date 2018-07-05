%%%%%%%%% Load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
clear all;
close all;
clc;
% Set the filenames of the TIFF files of Green and Red Channels
fname1 = '/Volumes/Schultz_group_data/Crazy Eights/Ann/GCaMP6 imaging/2018.03.07/2P/20180307_17_09_31_2P/20180307_17_09_31_2P_XYT_ch_3.tif'; %green channel
fname2 = '/Volumes/Schultz_group_data/Crazy Eights/Ann/GCaMP6 imaging/2018.03.07/2P/20180307_17_09_31_2P/20180307_17_09_31_2P_XYT_ch_4.tif'; %red channel
% Get the information of the TIFF files
info1 = imfinfo(fname1);
info2 = imfinfo(fname2);
% Set the number of frames to take from the TIFF files
num_images1 = 930; 
num_images2 = 930;
% Initialise the arrays for green and red channels of type 'single'
green_channel=zeros(512,512,num_images1,'single');
red_channel=zeros(512,512,num_images2,'single');
% Load the images and save the values in the array
for k = 1:num_images1
    [green_channel(:,:,k), map_green] = imread(fname1, k, 'Info', info1);
    [red_channel(:,:,k), map_red] = imread(fname2, k, 'Info', info2);
end
%Clear the other variables to save memory space
clearvars -except green_channel red_channel;
% Get the dimensions of each channel
dim = size(red_channel);
N = dim(3); % Number of frames
% Add the folder where the other scripts are located
addpath([pwd, '/dependencies/'])
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%% Compute the summary images %%%%%%%%%%%%%%%%%%%%%%
%
meanIm_red = mean(red_channel,3);
meanIm_green = mean(green_channel,3);
corrIm_red = crossCorr(red_channel);
corrIm_green = crossCorr(green_channel);
% The summary image can be from the mean or correlation above or 
% a composite image of red and green channels
summary_image = imread('Composite3_new.jpg'); 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%% Initialisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Set the initial image 
initial_image = medfilt2(corrIm_green, [5 5]);

% Set alpha parameter
radius = 10; % expected cell radius in pixels
alpha = 0.1; % recommended value: 0.5
init_opt.blur_radius = 4; % default is 1; radius of the blurring applied
                          % to the input summary image
% Secondary Metric
init_opt.secondary_metric = meanIm_red./meanIm_green;
retuned_opt.secondary_metric = init_opt.secondary_metric;
init_opt.second_alpha = 0.1;
% Initialise by getting the candidate ROI seeds
phi_0 = initialise(initial_image, radius, alpha, init_opt);

exp_ROIs = 150; %expected number of ROIs in the frame

retuned_alpha = alpha;
initial_masks_num = size(phi_0,3);
retuned_masks_num = initial_masks_num;
while (retuned_masks_num > exp_ROIs)
    retuned_alpha = retuned_alpha + 0.05;
    retuned_opt.second_alpha = retuned_alpha;
    disp(['Current Alpha is ', num2str(retuned_alpha)]);
    if retuned_alpha > 0.95;
        break
    end
    phi_0 = initialise(initial_image, radius, retuned_alpha, retuned_opt);
    retuned_masks_num = size(phi_0,3);
end

%}

% Set lambda parameter
retune_lambda  = 1; %Set to 0 if you want a predefine lambda, 1 if you want
                   %to retune it

% Algorithm parameters
seg_opt.lambda              = 30; % this depends on the data
seg_opt.mergeCorr           = 0.95; % correlation coefficient threshold
                                    % above which two neigbouring ROIs will
                                    % be merged
seg_opt.mergeDuring         = 1;
seg_opt.plot_progress       = 1;

if retune_lambda  
	lambda         = tune_lambda(phi_0, green_channel, red_channel, radius,...
                                 seg_opt, corrIm_green, meanIm_red);
    seg_opt.lambda = lambda;
end

tic;

%%%%%%%%%% Segmentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
[masks, cell_tsG, nhbd_tsG, cell_tsR, nhbd_tsR] = segment_ratio2(phi_0, ...
    green_channel, red_channel, radius, seg_opt);
runtime = toc
mask_num = size(masks,3); % Number of detected ROIs

% Display initial results
clear opts;
% Plot masks on summary image
opts.plot_ids = 1; %Set to 1 if you want to view the ID number of the ROIs
plotContoursOnSummaryImage(summary_image, masks, opts);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% Time Series Extraction %%%%%%%%%%%%%%%%%%%%%%%
%{
CaSeriesGreen = zeros(size(masks,3),N);
CaSeriesRed = zeros(size(masks,3),N);
sm_span = 7; %span of smoothing; if not declared, the default is 5
for i=1:mask_num
    CaSeriesGreen(i,:) = smooth(cell_tsG(i,:),sm_span);
    CaSeriesRed(i,:) = smooth(cell_tsR(i,:),sm_span);
end
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% Plot Ratiometric Signal %%%%%%%%%%%%%%%%%%%%%%
%
figure;
plotCa(cell_tsG, cell_tsR,7);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}