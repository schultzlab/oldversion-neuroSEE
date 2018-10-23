function neuroSEE_motioncorrect(dir, filenames, path_name, sessiondate, mouse)

parfor i=1:length(filenames)
    session = [sessiondate,'_', mouse, '_X', num2str(i+20)];
    try
        tic;
        disp(['Loading imaging file for <', session, '>.']);
        FileName = strcat(dir,filenames{i});
        video2p = read_file(FileName); 
        
        green_channel = video2p(:,:,1:2:end-1);
        red_channel = video2p(:,:,2:2:end);

        green_channel = single(green_channel); % convert to single precision 
        red_channel = single(red_channel); 
        T = size(green_channel,ndims(green_channel));

        green_channel = green_channel - min(green_channel(:));
        red_channel = red_channel - min(red_channel(:));

        % Motion Artefact Correction
        disp(['Initialising parameters for motion artefact correction for <', session, '>.']);
        options_rigid = NoRMCorreSetParms('d1',size(red_channel,1),'d2',size(red_channel,2),'bin_width',200,'max_shift',15,'us_fac',50,'init_batch',200);

        disp(['Performing motion artefact correction for <', session, '>.']);
        [red_corrected, green_corrected,~,~,~] = normcorre_seig(red_channel,green_channel,options_rigid); 
        disp(['Finished correction for <', session, '>.']);
        toc;

        parsave(green_corrected,red_corrected,[path_name, session, '_mcorrected']);
        disp(['Successfully saved motion-corrected file for <', session, '>.']);
        
    catch
        disp(['ERROR: Failed processing imaging file for <', session, '>.']);
    end
end