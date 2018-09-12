function R = ratiometric_Ca(green_timeseries, red_timeseries)

C = size(green_timeseries,1); %number of cells 
Tl = size(green_timeseries,2); %number of sample points in the time series
t = [1:Tl]*1/30;
T = max(t);
sm=7;
for c=1:C
    green_timeseries(c,:) = smooth(green_timeseries(c,:),sm);
    red_timeseries(c,:) = smooth(red_timeseries(c,:),sm);
end
R_ = green_timeseries./red_timeseries;

R0 = mode(R_');
R0_ = repmat(R0',[1 Tl]);
R = (R_-R0_)./R0_;

end