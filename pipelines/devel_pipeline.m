clear;

%% Basic setup
test = 1;       % set to 1 if testing, this will use one of smaller files in ../test
display = 1;    % set to 1 to display results (all results are saved in 
                %    individual file directory with a pdf summary regardless)
force = 0;      % set to 1 to overwrite saved processed files. This will 
                %    force pipeline to redo all steps incl. raw to tif
                %    conversion. If force = 0, processing step is skipped
                %    if output of said step already exists in individual
                %    file directory.

% gcp;           % start parallel pool
addpath(genpath('../ABLE'));
addpath(genpath('../extract_tracking'));
addpath(genpath('../motion_correction'));
addpath(genpath('../utilities'));

% Data location
data_locn = '/Volumes/thefarm2/live/CrazyEights/AD_2PCa/GCaMP6s/';
if ~exist(data_locn,'dir')
    data_locn = '/Volumes/RDS/project/thefarm2/live/CrazyEights/AD_2PCa/GCaMP6s/';
end
    
if test
    currdir = pwd;
    if currdir(length(currdir)-8:end)=='pipelines'
        data_locn = [currdir(1:length(currdir)-9) 'test'];
    elseif currdir(length(currdir)-5:end)=='latest'
        data_locn = [currdir '/test'];
    end
end

%% Read files
% This should be later changed to accept multiple files/directories

file = '20181016_10_11_35';

%% Load raw and save tif files 
% or load tif files if they exist

filedir_2P = [data_locn '/' file(1:8) '/2P/' file '_2P/'];
[imG,imR] = load_file(filedir_2P,file,force);

%% Dezipper and do motion correction
% Motion corrected tif stacks are saved along with summary of motion
% corrected images (fig and pdf)

scale = 1; % image downsampling factor
Nimg_ave = 14; % number of images to be averaged for calculating pixel shift (zippering)
[imG, imR] = neuroSEE_motionCorrect(imG, imR, filedir_2P, file, Nimg_ave, scale, display, force);

%% Use ABLE to extract ROIs

if test
    cellrad = 10;   % expected radius of a cell (pixels)
    maxcells = 70;  % estimated number of cells in FOV
end
[cell_tsG, cell_tsR, masks] = neuroSEE_segment( imG, mean(imR,3), filedir_2P, file, cellrad, maxcells, display, force );

%% Run FISSA to extract neuropil-corrected time-series


%% Extract tracking data

% Find tracking file that corresponds to 2P recording
filedir_track = [data_locn '/' file(1:8) '/Neurotar/'];
fname_track = findMatchingTrackingFile(filedir,file);

% Import tracking file 



csv_import(fname_track);

%% Generate place field maps