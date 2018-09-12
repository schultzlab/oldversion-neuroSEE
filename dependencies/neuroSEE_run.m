% Modify the directory, filename and session name here:
dir = '/Volumes/Schultz_group_data/Crazy Eights/Ann/GCaMP6 imaging/2018.07.18/2P/';
filename='20180718_10_15_00_2P/20180718_10_15_00_2P_XYT_1.tif';
session='01';

disp('Loading imaging file...');

FileName = strcat(dir,filename);
addpath([pwd, '/dependencies/'])
addpath([pwd, '/dependencies/motion_correct'])
addpath([pwd, '/dependencies/gLog_filtering'])
addpath([pwd, '/dependencies/moments'])

tic; 
video2p = read_file(FileName); 
toc;

green_channel = video2p(:,:,1:2:end-1);
red_channel = video2p(:,:,2:2:end);

green_channel = single(green_channel); % convert to single precision 
red_channel = single(red_channel); 
T = size(green_channel,ndims(green_channel));

green_channel = green_channel - min(green_channel(:));
red_channel = red_channel - min(red_channel(:));

%% Motion Artefact Correction
disp('Initialising parameters for motion artefact correction...');
options_rigid = NoRMCorreSetParms('d1',size(red_channel,1),'d2',size(red_channel,2),'bin_width',200,'max_shift',15,'us_fac',50,'init_batch',200);

tic;
disp('Performing motion artefact correction...');
[red_corrected, green_corrected,shifts_rigid,template_rigid,options_rigid] = normcorre_seig(red_channel,green_channel,options_rigid); 
disp('Finished correction.');
toc;
%}
%% Segmentation
cellRadius = 10;
maxCells = 100;

disp('Performing segmentation...');
[green_timeseries,red_timeseries,masks] = neuroSEE_segment(green_corrected,red_corrected,cellRadius,maxCells,session);
disp('Finished segmenting cells.');
%% Ratiometric Signals
R = ratiometric_Ca(green_timeseries, red_timeseries);
figure; 
plotCa2(R);
%% Place Field Analysis
spikes = neuroSEE_nnd(R);
n = 50;
CaFs = 30.92;
trackFs = 100;
mode = 1;

[Rmap, occMap, placeMap, nxy] = neuroSEE_PFmap(R,spikes,xy,phi,n,CaFs,trackFs,mode);

figure;
subplot(1,2,1)
imagesc(occMap)
c = colorbar();
c.Label.String = 'cumulative samples';
title('Occupancy Map')


subplot(1,2,2)
plot(nxy(:,1),nxy(:,2),'Color',[0,0.7,0.9])
title('Behavioural Track')
set(gca,'box','off');
set(gca,'visible','off');
camroll(-90)

nC = size(R,1);
figure('units','normalized','outerposition',[0 0 0.8 0.8])
title('Place Maps');
for i=1:nC
    %subplot(2+rem(nC,2),floor(nC/2),ip)
    subplot(ceil(nC/2),2,i)
    imagesc(Rmap(:,:,i)./occMap, [0 0.9])
    c = colorbar();
    c.Label.String = 'dF/F per sample';
    title(['cell ',num2str(i)])
end