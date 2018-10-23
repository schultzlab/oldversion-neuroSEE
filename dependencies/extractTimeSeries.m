function [tsG, tsR] = extractTimeSeries(mask, mods, video_reshapedG, video_reshapedR, t_len)
 
    mask([1,end], :) = 0;
    mask(:, [1,end]) = 0;
    loc              = find(mask);
    vid_loc          = bsxfun(@plus,loc, mods);
    raw_tsG           = video_reshapedG(vid_loc(:));
    raw_tsG           = reshape(raw_tsG, length(loc),t_len);
    tsG               = mean(raw_tsG,1, 'double'); %% much faster to calculate as double
    %%%added codes by svp
    raw_tsR = video_reshapedR(vid_loc(:));
    raw_tsR = reshape(raw_tsR, length(loc),t_len);
    tsR = mean(raw_tsR,1, 'double');
    %%%
end