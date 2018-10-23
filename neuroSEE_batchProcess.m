%% Initialisation

addpath([pwd,'/dependencies/'])
addpath([pwd,'/dependencies/motion_correct/'])
addpath([pwd,'/dependencies/gLog_filtering/'])

% Write the directory of where the TIFF videos will be accessed from
dir = '/Volumes/Schultz_group_data/Crazy Eights/Ann/GCaMP6 imaging/2018.10.13/2P/';

% Enumerate the filenames of the TIFF videos to be processed
video_names = {'20181013_11_24_00_2P/20181013_11_24_00_2P_XYT_1.tif',...
    '20181013_11_33_44_2P/20181013_11_28_44_2P_XYT_1.tif',...
    '20181013_11_41_34_2P/20181013_11_38_34_2P_XYT_1.tif',...
    '20181013_11_41_34_2P/20181013_11_38_34_2P_XYT_2.tif',...
    '20181013_11_45_40_2P/20181013_11_45_40_2P_XYT_1.tif',...
    '20181013_11_45_40_2P/20181013_11_45_40_2P_XYT_2.tif',...
    '20181013_11_50_23_2P/20181013_11_50_23_2P_XYT_1.tif',...
    '20181013_11_50_23_2P/20181013_11_50_23_2P_XYT_2.tif'};

% Enumerate the tracking data filenames
track_dir = '/Volumes/Schultz_group_data/Crazy Eights/Ann/GCaMP6 imaging/2018.10.16/Neurotar/';
track_names = {'Track_2018-10-16-10-11-40/SavedTrack-2018-10-16-10-11-40.mat'
    };

% Write the directory of where you want to save the processed files
path_name = [pwd, '/m62/'];

% Write the details of the session and the mouse
sessiondate = '161018'; %DDMMYY format
mouse = 'm62';

% Initialise parallel pool
parpool;
poolobj = gcp;
pctRunOnAll addpath([pwd,'/dependencies/'])
pctRunOnAll addpath([pwd,'/dependencies/motion_correct/'])
pctRunOnAll addpath([pwd,'/dependencies/gLog_filtering/'])

%% Motion Correction
neuroSEE_motioncorrect(dir, video_names, path_name, sessiondate, mouse);

%% ROI Segmentation

% Get a stack of the first 500 frames of each of the videos 
green_stack = [];
red_stack = [];
for i = 1:length(video_names)
    fname = [path_name,sessiondate,'_',mouse,'_X',num2str(i),'_mcorrected.mat'];
    load(fname,'green_corrected','red_corrected');
    green_stack = cat(3,green_stack,green_corrected(:,:,1:500));
    red_stack = cat(3,red_stack,red_corrected(:,:,1:500));
end

clearvars green_corrected red_corrected;

% Perform segmentation on these stacks to get the masks
[~,~,masks] = neuroSEE_segment(green_stack,red_stack,cellRadius,maxCells,[]);

% Initialise data structures, i.e. names, fields, etc
m62.fam1fam2.neuron = cell(size(masks,3),1); %change the structure name based on the mouse name, environment, etc.

% Upon getting the masks, extract the timeseries from each of the videos
tsG = [];
tsR = [];
for i = 1:length(video_names)
    fname = [path_name,sessiondate,'_',mouse,'_X',num2str(i),'_mcorrected.mat'];
    load(fname,'green_corrected','red_corrected');
    timeSeriesG = [];
    timeSeriesR = [];
    mods = (0:1:(size(green_corrected,3)-1))*size(green_corrected,1)*size(green_corrected,2);
    video_reshapedG = reshape(green_corrected,[size(green_corrected,1)*size(green_corrected,2)*size(green_corrected,3),1,1]);
    video_reshapedR = reshape(red_corrected,[size(red_corrected,1)*size(red_corrected,2)*size(red_corrected,3),1,1]);
    for ii = 1:size(masks,3)
        currentMask = masks(:,:,ii);
        [timeSeriesG(ii,:), timeSeriesR(ii,:)] = extractTimeSeries(currentMask, mods, video_reshapedG, video_reshapedR, size(green_corrected,3));
    end  
    % concatenate the timeseries:
    tsG = cat(2,tsG,timeSeriesG);
    tsR = cat(2,tsR,timeSeriesR);
end

% Edit the structure names below !!!!
for ii = 1:size(masks,3)
    m62.fam1fam2.neuron{ii}.CaTimeSeriesGreen = tsG(ii,:); %timeseries from green channel
    m62.fam1fam2.neuron{ii}.CaTimeSeriesRed = tsR(ii,:); %timeseries from red channel 
    m62.fam1fam2.neuron{ii}.deltaRoverR = ratiometric_Ca(tsG(ii,:),tsR(ii,:)); %deltaR/R
    m62.fam1fam2.neuron{ii}.spikes = neuroSEE_nnd(m62.fam1fam2.neuron{ii}.deltaRoverR,0.05,2); % non-negative deconv
end

%save('m62.mat', 'm62', '-v7.3');    %// UNCOMMENT THIS IF YOU WANT TO
                                     %// SAVE PRIOR TO PLACE FIELD ANALYSIS


%% Process tracking data
neuroSEE_trackProcess(track_dir,track_names,path_name,sessiondate,mouse);

track_xy = [];
track_phi = [];
for i = 1:length(video_names)
    fname = [path_name,sessiondate,'_',mouse,'_X',num2str(i),'_track.mat'];
    load(fname,'xy','phi');
    track_xy = cat(1,track_xy,xy);
    track_phi = cat(1,track_phi,phi);
end

% Edit the structure names below !!!!
m62.fam1fam1inv.spatialx = track_xy(:,1);
m62.fam1fam1inv.spatialy = track_xy(:,2);
m62.fam1fam1inv.spatialphi = track_phi;

%% Place Field Analysis

% Initialise parameters
n = 36; %spatial bins
CaFs = 30.92; %calcium imaging sampling rate
trackFs = 100; %tracking data sampling rate
mode = 3; % 1 if 1d place maps, 2 if 2d place maps, 3 if 1d but reduced radial bins
Rmap = [];
occMap = [];
placeMap = [];
xy = [m62.fam1fam1inv.spatialx m62.fam1fam1inv.spatialy];
phi = m62.fam1fam1inv.spatialphi;

for ii = 1:length(m62.fam1fam1inv.neuron)
    R = m62.fam1fam1inv.neuron{ii}.deltaRoverR;
    spikes = m62.fam1fam1inv.neuron{ii}.spikes;
    if mode == 1 || mode == 3
        [Rmap(ii,:), occMap, placeMap(ii,:), ~] = neuroSEE_PFmap(R,spikes,xy,phi,n,CaFs,trackFs,mode);
    else
        [Rmap(:,:,ii), occMap_temp(:,:,ii), placeMap(:,:,ii), ~] = neuroSEE_PFmap(R,spikes,xy,phi,n,CaFs,trackFs,mode);
    end
end

figure; imagesc(placeMap); colormap(jet)