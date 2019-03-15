% Written by Ann Go
%
% This function implements Katie's motionCorrectToNearestPixel fxn if 
% force = 1 or if either green or red motion corrected tif files don't yet
% exist in filedir.
%
% INPUTS
%   imG : matrix of green image stack
%   imR : matrix of red image stack
%   filedir : directory of image stacks
%   file : part of file name of image stacks in the format
%           yyyymmdd_hh_mm_ss
%   Nimg_ave : number of images to be averaged when correcting for pixel
%                shift (zippering)
%   scale : downsampling factor for motion correction
%   display : if =1, mean motion corrected images for green and red channels
%               will be displayed. Image stack are saved in filedir
%               regardless.
%   force : if =1, motion correction will be done even though motion
%             corrected images already exist
% OUTPUTS
%   imG : matrix of motion corrected green image stack
%   imR : matrix of motion corrected red image stack
%
% Other outputs of motionCorrectToNearestPixel
% [stack_g, stack_r, output, fh] = motionCorrectToNearestPixel(...)
%   output : cell array containing
%               output.green & output.red
%            each contains
%               output.<colour>.template : template used for registration
%               output.<colour>.meanframe : mean frame for original stack 
%               output.<colour>.meanregframe : mean frame for registered stack
%               output.<colour>.shift : matrix of x&y shift for each frame
%   fh : figure summarising comparison between mean frames for original and
%           registered stacks for green and red channels

function [imG, imR] = neuroSEE_motionCorrect(imG, imR, filedir, file, Nimg_ave, scale, display, force)
    if nargin<8, force = 0;      end
    if nargin<7, display = 1;    end
    if nargin<6, scale = 1;      end
    if nargin<5, Nimg_ave = 10;  end

    % If asked to force overwrite, run motion correction right away
    if force
        [imG, imR] = motionCorrectToNearestPixel(double(imG), double(imR), scale, display, Nimg_ave, filedir, file);
    else
        fname_tif_gr_mcorr = [filedir file '_2P_XYT_green_mcorr.tif'];
        fname_tif_red_mcorr = [filedir file '_2P_XYT_red_mcorr.tif'];
        fname_mat_mcorr = [filedir file '_2P_mcorr_output.mat'];
        fname_fig = [filedir file '_2P_mcorr_summary.fig'];
        yn_gr_mcorr = exist(fname_tif_gr_mcorr,'file');
        yn_red_mcorr = exist(fname_tif_red_mcorr,'file');
        yn_mat_mcorr = exist(fname_mat_mcorr,'file');
        yn_fig_mcorr = exist(fname_fig,'file');

        % If any of motion corrected tif stacks or motion correction output
        % mat doesn't exist, run motion correction
        if any([~yn_gr_mcorr,~yn_red_mcorr,~yn_mat_mcorr])
            [imG, imR] = motionCorrectToNearestPixel(double(imG), double(imR), scale, display, Nimg_ave, filedir, file);
        else
            % If they do exist, load motion corrected tif stacks
            [imG, imR] = load_file(filedir,file,force,'_mcorr');
            % If summary fig exists, open it if asked to display
            if yn_fig_mcorr
                if display, openfig(fname_fig,'new'); end
                 
            else
                % If summary fig doesn't exist, create it   
                output = load(fname_mat_mcorr);
                out_g = output.green;
                out_r = output.red;
                fh = figure; 
                subplot(221), 
                    imagesc( out_g.meanframe ); 
                    axis image; colorbar; axis off;
                    title( 'Mean frame for raw green' );
                subplot(222), 
                    imagesc( out_g.meanregframe ); 
                    axis image; colorbar; axis off; 
                    title( 'Mean frame for corrected green' );
                subplot(223), 
                    imagesc( out_r.meanframe ); 
                    axis image; colorbar; axis off; 
                    title( 'Mean frame for raw red' );
                subplot(224), 
                    imagesc( out_r.meanregframe ); 
                    axis image; colorbar; axis off;
                    title( 'Mean frame for corrected red' );
                fname_fig = [filedir file '_2P_mcorr_summary'];
                    savefig( fh, fname_fig, 'compact' );
                    saveas( fh, fname_fig, 'pdf' );
                if ~display, close( fh ); end
            end
        end
    end
