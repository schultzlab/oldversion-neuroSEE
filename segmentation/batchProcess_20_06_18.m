% Day 2

dir = '/Volumes/Schultz_group_data/Crazy Eights/Ann/GCaMP6 imaging/2018.03.06/2P/';


try
    filename = '20180306_16_57_28_2P/20180306_16_57_28_2P_XYT.tif';
    sumImage = 'Composite2nd.jpg';
    tic;
    neuroSEE_batch(dir,filename,sumImage);
    toc;
catch
    disp(['ERROR: Failed processing imaging file <', filename, '>.']);
end

try
    filename = '20180306_17_04_06_2P/20180306_17_04_06_2P_XYT.tif';
    sumImage = 'Composite2nd.jpg';
    tic;
    neuroSEE_batch(dir,filename,sumImage);
    toc;
catch
    disp(['ERROR: Failed processing imaging file <', filename, '>.']);
end


try
    filename = '20180306_17_09_09_2P/20180306_17_09_09_2P_XYT.tif';
    sumImage = 'Composite2nd.jpg';
    tic;
    neuroSEE_batch(dir,filename,sumImage);
    toc;
catch
    disp(['ERROR: Failed processing imaging file <', filename, '>.']);
end

