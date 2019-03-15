%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NeuroSEE: An automated Neuronal Source Extraction
%             and Exploration toolbox
%   
%   Author: Seigfred Prado   
%   Supervisor: Simon Schultz
%   Acknowledgment: Stephanie Reynolds, Pier Luigi Dragotti
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Hybrid Level Set Segmentation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   INPUTS:
%       greench: green channel matrix
%       redch: red channel matrix
%       cellRadius: expected radius of the cell in pixels
%       maxCells: estimated maximum number of cells in a frame
%       session: session number or date/time
%   OUTPUTS:
%       cell_tsG: green channel time series
%       cell_tsR: red channel time series
%       masks: masks of the segmented ROIs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [cell_tsG, cell_tsR, masks] = neuroSEE_segment(greench,redch,cellRadius,maxCells,session)

%% Initialise the summary images

% Compute the summary images, i.e. mean and correlation images
meanIm_red = mean(redch,3);
meanIm_green = mean(greench,3);
corrIm_red = crossCorr(redch);
corrIm_green = crossCorr(greench);

summary_image = meanIm_green./meanIm_red; 
% Set the initial images 
green_initial_image = medfilt2(corrIm_green, [5 5]);
red_initial_image = summary_image;
red_initial_image = red_initial_image - min(red_initial_image(:));
red_initial_image = red_initial_image/max(red_initial_image(:));

% Set initial parameters
if nargin<3
    cellRadius = 10;
end
if nargin<4
    maxCells = 100;
end
if nargin<5
    session = 'test';
end

radius = cellRadius; 
alpha = 0.05; % starts with a low value then increase that during tuning
init_opt.blur_radius = 4; % default is 1; radius of the blurring applied
                          % to the input summary image

% Secondary Metric
% Divide the red over the green average image so that the only remaining
% signal is from neuronal nuclei, thus eliminating field illumination
% inhomogeneity and non-cellular structures like blood vessels
init_opt.secondary_metric = red_initial_image;
%init_opt.secondary_metric = meanIm_red./meanIm_green;
retuned_opt.secondary_metric = init_opt.secondary_metric;
init_opt.second_alpha = 0.5; % same as alpha

%% Initialise by getting the candidate ROI seeds using the green channel
% This is activity-based initialisation
phi_0 = initialise(green_initial_image, radius, alpha, init_opt);
exp_ROIs = maxCells; 

%% Retune the alpha parameter
retuned_alpha = alpha;
initial_masks_num = size(phi_0,3);
retuned_masks_num = initial_masks_num;
while (retuned_masks_num > exp_ROIs)
    retuned_alpha = retuned_alpha + 0.05;
    retuned_opt.second_alpha = retuned_alpha;
    disp(['Alpha Tuning: Current Alpha is ', num2str(retuned_alpha)]);
    phi = initialise(green_initial_image, radius, retuned_alpha, retuned_opt);
    retuned_masks_num = size(phi,3);
end
disp(['Final alpha value is ', num2str(retuned_alpha)]);
phi_0_g=phi;

%% Initialise by getting the candidate ROI seeds using the red channel
% This is morphology-based initialisation
I_norm_ori = red_initial_image;
min_sigma = 4;
max_sigma = 8;
Ntheta = 8;
C = 0;
tm = 0;
p2 = 2*max_sigma*3+1;
bw = adaptivethreshold(I_norm_ori,p2,C,tm);

[new_coordinate_x,new_coordinate_y,absigma_leave]=SeedDetection(I_norm_ori,min_sigma,max_sigma,bw,Ntheta);

MinArea = 6;


ParametersOutput=fastSegmentation(I_norm_ori,new_coordinate_x,new_coordinate_y,absigma_leave,MinArea,Ntheta,min_sigma);

phi_0_r = ones(512,512,length(ParametersOutput.xypos));
for k=1:length(ParametersOutput.Pixels)
    for j=1:length(ParametersOutput.Pixels{k})
        phi_0_r(ParametersOutput.Pixels{k}(2,j),ParametersOutput.Pixels{k}(1,j),k) = -1;
    end 
end

phi_0 = cat(3,phi_0_g,phi_0_r); 


%% Retune the lambda parameter
% Initialise the level set algorithm parameters
seg_opt.lambda = 60; % this depends on the data
seg_opt.mergeCorr = 0.95; % correlation coefficient threshold
                                    % above which two neigbouring ROIs will
                                    % be merged
seg_opt.mergeDuring = 1;
retune_lambda = 0; %make this 1 if you want to retune lambda

if retune_lambda
    % Tune lambda based on the data
    [retuned_lambda, phi] = lambda_tune(phi_0, greench, redch, radius,...
                                     seg_opt, corrIm_green, red_initial_image);
    seg_opt.lambda = retuned_lambda;
else
    phi = phi_0;
end


tic;

%% Segmentation proper
seg_opt.corrIm = corrIm_green;
seg_opt.meanIm = red_initial_image;
%seg_opt = rmfield(seg_opt,'plot_progress');

[masks, cell_tsG, nhbd_tsG, cell_tsR, nhbd_tsR] = segment_cells(phi, ...
    greench, redch, radius, seg_opt);
runtime = toc;
mask_num = size(masks,3); % Number of detected ROIs


%% Display initial results
clear opts;
% Plot masks on summary image
opts.plot_ids = 1; %Set to 1 if you want to view the ID number of the ROIs
plotContoursOnSummaryImage(summary_image, masks, opts);

%% Saving the .mat file of the segmentation workspace variables
save(session, '-v7.3');
end 