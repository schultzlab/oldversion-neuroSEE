% Written by Ann Go
%
% This function implement the ABLE fxn (renamed version of neuroSEE_segment
% which Katie adapted from Seig's code) when force = 1 or if figure with 
% ROIs don't yet exist in filedir.
%
% INPUTS
%   imG : matrix of green image stack
%   imR : matrix of red image stack
%   filedir : directory of image stacks
%   file : part of file name of image stacks in the format
%           yyyymmdd_hh_mm_ss
%   cellrad : expected cell radius in pixels
%   maxcells : expected maximum number of cells in FOV
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

function [cell_tsG, cell_tsR, masks] = neuroSEE_segment(stack_g, mean_r, filedir, file, cellrad, maxcells, display, force)
    if nargin<8, force = 0;      end
    if nargin<7, display = 1;    end
    if nargin<6, maxcells = 200; end
    if nargin<5, cellrad = 10;   end

    % If asked to force overwrite, run ABLE right away
    if force
        [cell_tsG, cell_tsR, masks] = ABLE( stack_g, mean_r, filedir, file, cellrad, maxcells, display );
    else
        fname_ts = [filedir file '_2P_segment_output.mat'];
        fname_fig = [filedir file '_2P_ROIs'];
        yn_ts = exist(fname_ts,'file');
        yn_fig = exist(fname_fig,'file');
        
        % If timeseries mat file doesn't exist, run ABLE
        if ~yn_ts
            [cell_tsG, cell_tsR, masks] = ABLE( stack_g, mean_r, filedir, file, cellrad, maxcells, display );
        else
            % If it exists, load it 
            segmentOutput = load(fname_ts);
            cell_tsG = segmentOutput.cell_tsG;
            cell_tsR = segmentOutput.cell_tsR;
            masks = segmentOutput.masks;
            % If ROI figure exists, open it if asked to display
            if yn_fig
                if display, openfig(fname_fig,'new'); end
            % If it doesn't exist, create & save figure
            else
               plotopts.plot_ids = 1; % set to 1 to view the ID number of the ROIs on the plot
               fig = plotContoursOnSummaryImage(mean_R, masks, plotopts);
               fname_fig = [filedir file '_2P_ROIs'];
               savefig(fig,fname_fig,'compact');
               saveas(fig,fname_fig,'pdf');
               if ~display, close(fig); end
            end
        end
    end
