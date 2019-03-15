%% Initialise names
%{
dir = '/Volumes/Schultz_group_data/Crazy Eights/Ann/GCaMP6 imaging/2018.07.18/2P/';

video_names = {'20180718_14_34_06_2P/20180718_14_34_06_2P_XYT_1.tif',...
    '20180718_14_34_06_2P/20180718_14_34_06_2P_XYT_2.tif',...
    '20180718_15_18_50_2P/20180718_15_18_50_2P_XYT_1.tif',...
    '20180718_15_18_50_2P/20180718_15_18_50_2P_XYT_2.tif',...
    '20180718_15_36_01_2P/20180718_15_36_01_2P_XYT_1.tif',...
    '20180718_15_36_01_2P/20180718_15_36_01_2P_XYT_2.tif',...
    '20180718_15_41_53_2P/20180718_15_41_53_2P_XYT_1.tif',...
    '20180718_15_41_53_2P/20180718_15_41_53_2P_XYT_2.tif',...
    '20180718_15_50_54_2P/20180718_15_50_54_2P_XYT_1.tif',...
    '20180718_15_50_54_2P/20180718_15_50_54_2P_XYT_2.tif'};

track_names={'/Volumes/Schultz_group_data/Crazy Eights/Ann/GCaMP6 imaging/2018.07.18/Neurotar tracker/Track_2018-07-18-14-34-41/SavedTrack-2018-07-18-14-34-41.mat',...
            '/Volumes/Schultz_group_data/Crazy Eights/Ann/GCaMP6 imaging/2018.07.18/Neurotar tracker/Track_2018-07-18-15-18-58/SavedTrack-2018-07-18-15-18-58.mat',...
            '/Volumes/Schultz_group_data/Crazy Eights/Ann/GCaMP6 imaging/2018.07.18/Neurotar tracker/Track_2018-07-18-15-36-20/SavedTrack-2018-07-18-15-36-20.mat',...
            '/Volumes/Schultz_group_data/Crazy Eights/Ann/GCaMP6 imaging/2018.07.18/Neurotar tracker/Track_2018-07-18-15-42-09/SavedTrack-2018-07-18-15-42-09.mat',...
            '/Volumes/Schultz_group_data/Crazy Eights/Ann/GCaMP6 imaging/2018.07.18/Neurotar tracker/Track_2018-07-18-15-51-10/SavedTrack-2018-07-18-15-51-10.mat'};

path_name = '/Volumes/Schultz_group_data/Seigfred/Processed Files/m50/';
%}

dir = '/Volumes/Schultz_group_data/Crazy Eights/Ann/GCaMP6 imaging/2018.08.31/2P/';

video_names = {'20180831_15_10_49_2P/20180831_15_10_49_2P_XYT_1.tif',...
    '20180831_15_10_49_2P/20180831_15_10_49_2P_XYT_2.tif',...
    '20180831_15_15_24_2P/20180831_15_15_24_2P_XYT_1.tif',...
    '20180831_15_15_24_2P/20180831_15_15_24_2P_XYT_2.tif',...
    '20180831_15_21_03_2P/20180831_15_21_03_2P_XYT_1.tif',...
    '20180831_15_21_03_2P/20180831_15_21_03_2P_XYT_2.tif',...
    '20180831_15_26_05_2P/20180831_15_26_05_2P_XYT_1.tif',...
    '20180831_15_26_05_2P/20180831_15_26_05_2P_XYT_2.tif',...
    '20180831_15_30_39_2P/20180831_15_30_39_2P_XYT_1.tif',...
    '20180831_15_30_39_2P/20180831_15_30_39_2P_XYT_2.tif',...
    '20180831_15_35_14_2P/20180831_15_35_14_2P_XYT_1.tif',...
    '20180831_15_35_14_2P/20180831_15_35_14_2P_XYT_2.tif',...
    '20180831_15_40_03_2P/20180831_15_40_03_2P_XYT_1.tif',...
    '20180831_15_40_03_2P/20180831_15_40_03_2P_XYT_2.tif',...
    '20180831_15_45_18_2P/20180831_15_45_18_2P_XYT_1.tif',...
    '20180831_15_45_18_2P/20180831_15_45_18_2P_XYT_2.tif',...
    '20180831_16_10_58_2P/20180831_16_10_58_2P_XYT_1.tif',...
    '20180831_16_10_58_2P/20180831_16_10_58_2P_XYT_2.tif',...
    '20180831_16_23_08_2P/20180831_16_23_08_2P_XYT_1.tif',...
    '20180831_16_23_08_2P/20180831_16_23_08_2P_XYT_2.tif',...
    '20180831_16_28_21_2P/20180831_16_28_21_2P_XYT_1.tif',...
    '20180831_16_28_21_2P/20180831_16_28_21_2P_XYT_2.tif'};
%{
track_names={'/Volumes/Schultz_group_data/Crazy Eights/Ann/GCaMP6 imaging/2018.07.18/Neurotar tracker/Track_2018-07-18-14-34-41/SavedTrack-2018-07-18-14-34-41.mat',...
            '/Volumes/Schultz_group_data/Crazy Eights/Ann/GCaMP6 imaging/2018.07.18/Neurotar tracker/Track_2018-07-18-15-18-58/SavedTrack-2018-07-18-15-18-58.mat',...
            '/Volumes/Schultz_group_data/Crazy Eights/Ann/GCaMP6 imaging/2018.07.18/Neurotar tracker/Track_2018-07-18-15-36-20/SavedTrack-2018-07-18-15-36-20.mat',...
            '/Volumes/Schultz_group_data/Crazy Eights/Ann/GCaMP6 imaging/2018.07.18/Neurotar tracker/Track_2018-07-18-15-42-09/SavedTrack-2018-07-18-15-42-09.mat',...
            '/Volumes/Schultz_group_data/Crazy Eights/Ann/GCaMP6 imaging/2018.07.18/Neurotar tracker/Track_2018-07-18-15-51-10/SavedTrack-2018-07-18-15-51-10.mat'};
%}
path_name = '/Users/svp15/neurosee/Pipeline/m57/';
%{
%% Process tracking data
j=0;
for i=1:5
    j = j+1;
    try
        load(track_names{i});
        if mod(size(X,1),2)==0
            x1 = X(2:size(X,1)/2+1,1);
            x2 = X(size(X,1)/2+1:end,1);
            y1 = Y(2:size(Y,1)/2+1,1);
            y2 = Y(size(Y,1)/2+1:end,1);
            xy = [x1 y1];
            phi_1 = phi(2:size(X,1)/2+1,1);
            phi_2 = phi(size(X,1)/2+1:end,1);
            phi = phi_1;
            save([path_name,'180718_m50_X', num2str(j), '_track'], '-v7.3');
            j = j+1;
            xy = [x2 y2];
            phi = phi_2;
            save([path_name,'180718_m50_X', num2str(j), '_track'], '-v7.3');
        else
            x1 = X(2:size(X,1)/2+1,1);
            x2 = X(size(X,1)/2+2:end,1);
            y1 = Y(2:size(Y,1)/2+1,1);
            y2 = Y(size(Y,1)/2+2:end,1);
            xy = [x1 y1];
            phi_1 = phi(2:size(X,1)/2+1,1);
            phi_2 = phi(size(X,1)/2+2:end,1);
            phi = phi_1;
            save([path_name,'180718_m50_X', num2str(j), '_track'], '-v7.3');
            j = j+1;
            xy = [x2 y2];
            phi = phi_2;
            save([path_name,'180718_m50_X', num2str(j), '_track'], '-v7.3');
        end
    catch
        disp(['ERROR: Failed processing data for <', track_names{i}(end-22:end), '>.']);
    end
end
%
%% Motion-correct the imaging data
for i=1:22
%for i=4:length(video_names)
    session = ['310818_m57_X', num2str(i)];
    try
        tic;
        disp('Loading imaging file...');
        FileName = strcat(dir,video_names{i});
        video2p = read_file(FileName); 
        
        green_channel = video2p(:,:,1:2:end-1);
        red_channel = video2p(:,:,2:2:end);

        green_channel = single(green_channel); % convert to single precision 
        red_channel = single(red_channel); 
        T = size(green_channel,ndims(green_channel));

        green_channel = green_channel - min(green_channel(:));
        red_channel = red_channel - min(red_channel(:));

        % Motion Artefact Correction
        disp('Initialising parameters for motion artefact correction...');
        options_rigid = NoRMCorreSetParms('d1',size(red_channel,1),'d2',size(red_channel,2),'bin_width',200,'max_shift',15,'us_fac',50,'init_batch',200);

        disp('Performing motion artefact correction...');
        [red_corrected, green_corrected,shifts_rigid,template_rigid,options_rigid] = normcorre_seig(red_channel,green_channel,options_rigid); 
        disp('Finished correction.');
        
        nnY = quantile(red_channel(:),0.005);
        mmY = quantile(red_channel(:),0.995);

        [cY,mY,vY] = motion_metrics(red_channel,10);
        [cM1,mM1,vM1] = motion_metrics(red_corrected,10);
        T = length(cY);

        toc;
        savedir = [path_name, session];
        save([savedir,'_mcorrected'],'-regexp', '^(?!(video2p|red_channel|green_channel)$).','-v7.3'); 
        disp(['Successfully saved motion-corrected file for <', session, '>.']);
    catch
        disp(['ERROR: Failed processing imaging file for <', session, '>.']);
    end
    %clearvars -except dir video_names path_name session i;
end
%}
%% ROI Segmentation
%{
for j = [1 3:6 8:18 21:22]
%for j=1:length(video_names)
    path_name = '/Users/svp15/neurosee/Pipeline/m57/';
    %path_name = '/Volumes/Schultz_group_data/Seigfred/Processed Files/m50/';
    session_name = ['310818_m57_X', num2str(j), '_timeseries'];
    try
        session = ['310818_m57_X', num2str(j)];
        tic;
        disp('Loading imaging file...');
        FileName = [path_name,'310818_m57_X', num2str(j),'_mcorrected.mat'];
        load(FileName, 'green_corrected', 'red_corrected');
        %clearvars -except green_corrected red_corrected session session_name j path_name;
        
        % Segmentation
        cellRadius = 10;
        maxCells = 100; 

        disp('Performing segmentation...');
        [green_timeseries,red_timeseries,masks] = neuroSEE_segment(green_corrected,red_corrected,cellRadius,maxCells,session);
        disp('Finished segmenting cells.');

        savedir = [path_name, session_name];
        save(savedir, 'green_timeseries','red_timeseries','masks', '-v7.3');
        disp(['Successfully saved processed file for <', session_name, '>.']);
        %clear all;
    catch
        disp(['ERROR: Failed processing imaging file for <', session_name, '>.']);
    end
end
%}


%{
%% Cell natching across different sessions
% Set the directory name where the files are accessed
path_name = '/Volumes/Schultz_group_data/Seigfred/Processed Files/m50/';
% Set the file which will serve as the basis of cell matching
filename1 = [path_name, '180718_m50_X3_timeseries.mat'];
load(filename1,'masks');
masks_X1 = masks;

features_X1 = cell(1,size(masks_X1,3));
for i=1:size(masks_X1,3)
    temp = regionprops(masks_X1(:,:,i),'Area','Centroid');
    features_X1{i} = [temp(1).Area temp(1).Centroid];   
end
% For the succeeding imaging files:
for t=4:10
    try
    filename2 = [path_name, '180718_m50_X', num2str(t), '_timeseries.mat'];

    load(filename2, 'masks');
    masks_X2 = masks;
    %clearvars -except features_X1 masks_X1 masks_X2 filename1 filename2;
    
    features_X2 = cell(1,size(masks_X2,3));
    for i=1:size(masks_X2,3)
        temp = regionprops(masks_X2(:,:,i),'Area','Centroid');
        features_X2{i} = [temp(1).Area temp(1).Centroid];   
    end
    
    if size(masks_X2,3)<size(masks_X1,3)
        masks_X2_matched = zeros(size(masks_X2));
        ind = zeros(1,size(masks_X2,3));
    else
        masks_X2_matched = zeros(size(masks_X1));
        ind = zeros(1,size(masks_X1,3));
    end
    
    for i=1:length(ind)
        tempdiff=zeros(1,size(masks_X2,3));
        for j=1:size(masks_X2,3)
            tempdiff(j) =  pdist2(features_X1{i}(1,2:3),features_X2{j}(1,2:3));
        end
        ind(i) = find(tempdiff==min(tempdiff));
        masks_X2_matched(:,:,i) = masks_X2(:,:,ind(i));    
    end
    
    filename3 = [filename2(1:end-4), '_matched'];
    save(filename3, '-v7.3');
    %clearvars -except features_X1 masks_X1 filename1;
    catch
        disp(['ERROR: Failed matching for <', filename2, '>.']);
    end
end
%}

%% Place field analysis
%{
base_name = '/Volumes/Schultz_group_data/Seigfred/Processed Files/m50/180718_m50_X';
load([base_name, '3_timeseries.mat'],'green_timeseries','red_timeseries');
m=20:size(green_timeseries,1)-20;
for k=3:10
    session_name = [base_name, num2str(k), '_PF'];
    try
        if k == 3
            load([base_name, num2str(k),'_track.mat']);
            R = ratiometric_Ca(green_timeseries(m,:), red_timeseries(m,:));
        else
            load([base_name, num2str(k),'_timeseries.mat'],'green_timeseries','red_timeseries');
            load([base_name, num2str(k),'_track.mat']);
            load([base_name, num2str(k),'_timeseries_matched.mat'],'ind');
            R = ratiometric_Ca(green_timeseries(ind(m),:), red_timeseries(ind(m),:));
        end
        
        spikes = neuroSEE_nnd(R, 0.022);
        n = 50;
        CaFs = 30.92;
        trackFs = 100;
        mode = 1;
        
        [Rmap, occMap, placeMap, nxy] = neuroSEE_PFmap(R,spikes,xy,phi,n,CaFs,trackFs,mode);
        
        save(session_name,'R','spikes','xy','Rmap', 'occMap', 'placeMap', 'nxy', 'phi', '-v7.3');
        %clearvars -except m k;
    catch
        disp(['ERROR: Failed processing data for <', session_name, '>.']);
    end
end
%}