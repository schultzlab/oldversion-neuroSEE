% Written by Ann Go
%
% INPUTS
%   filedir : directory where file is
%   file :  part of filename of 2P image in the format
%             yyyymmdd_hh_mm_ss
%   force : set to 1 to overwrite existing tif files (not motion corrected)
%   suffix : =[]         to load original raw or tif files
%            ='mcorr'    to load motion corrected tif files   
% OUTPUTS
%   imG : matrix of green channel image stack
%   imR : matrix of red channel image stack

function [imG,imR] = load_file(filedir,file,force,suffix)
    if nargin<4, suffix  = [];    end
    if nargin<3, force  = 0;    end
    

    str = sprintf('Loading %s images...', suffix);
    cprintf(str)
    fname_tif_gr = [filedir file '_2P_XYT_green' suffix '.tif'];
    fname_tif_red = [filedir file '_2P_XYT_red' suffix '.tif'];
    fname_raw = [filedir file '_2P_XYT.raw'];


    % Check if tif stacks for green and red channel already exist. 
    % If yes: load tif stack/s
    %    no: load raw and create tif stack/s
    % But create tif stacks regardless if force = 1

    if force % load raw and create tif stacks
        options.skipN = 2; %skip every other frame
        imG = read_file( fname_raw, 1, Inf, options );
            writeTifStack( imG, fname_tif_gr );
        imR = read_file( fname_raw, 2, Inf, options );
            writeTifStack( imR, fname_tif_red );
        % This is somehow faster than reading entire raw file, saving half of
        % it as imG, the other half as imR and then creating tif stacks.
    else
        % Find out if either green or red tif stack exists
        yn_tif_gr = exist(fname_tif_gr,'file');
        yn_tif_red = exist(fname_tif_red,'file');

        % If both exist, load tif stacks
        if yn_tif_gr && yn_tif_red
            imG = read_file( fname_tif_gr );
            imR = read_file( fname_tif_red );
            % If the sizes of green and red stacks don't match, load raw
            if ~all(size(imG)==size(imR))
                options.skipN = 2; %skip every other frame
                imG = read_file( fname_raw, 1, Inf, options );
                    writeTifStack( imG, fname_tif_gr );
                imR = read_file( fname_raw, 2, Inf, options );
                    writeTifStack( imR, fname_tif_red );
            end
            
        % If only green tif exists, load red from raw and create red tif 
        elseif yn_tif_gr && ~yn_tif_red
            imG = read_file( fname_tif_gr );
            options.skipN = 2; %skip every other frame
            imR = read_file( fname_raw, 2, Inf, options );
                writeTifStack( imR, fname_tif_red );
            % If the sizes of green and red stacks don't match, load raw
            if ~all(size(imG)==size(imR))
                options.skipN = 2; %skip every other frame
                imG = read_file( fname_raw, 1, Inf, options );
                    writeTifStack( imG, fname_tif_gr );
                imR = read_file( fname_raw, 2, Inf, options );
                    writeTifStack( imR, fname_tif_red );
            end

        % If only red tif exists, load green from raw and create green tif   
        elseif ~yn_tif_gr && yn_tif_red
            imR = read_file( fname_tif_red );
            options.skipN = 2; %skip every other frame
            imG = read_file( fname_raw, 1, Inf, options );
                writeTifStack( imG, fname_tif_gr );
            % If the sizes of green and red stacks don't match, load raw
            if ~all(size(imG)==size(imR))
                options.skipN = 2; %skip every other frame
                imG = read_file( fname_raw, 1, Inf, options );
                    writeTifStack( imG, fname_tif_gr );
                imR = read_file( fname_raw, 2, Inf, options );
                    writeTifStack( imR, fname_tif_red );
            end

        % If no tif stacks exists, load all from raw and create tif stacks
        else
            options.skipN = 2; %skip every other frame
            imG = read_file( fname_raw, 1, Inf, options );
                writeTifStack( imG, fname_tif_gr );
            imR = read_file( fname_raw, 2, Inf, options );
                writeTifStack( imR, fname_tif_red );
        end
    end

    newstr = sprintf('%s Images loaded\n', suffix);
    refreshdisp(  newstr,str );
