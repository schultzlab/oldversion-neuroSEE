%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NeuroSEE: An automated Neuronal Source Extraction
%             and Exploration toolbox
%   
%   Author: Seigfred Prado   
%   Supervisor: Simon Schultz
%   Acknowledgment: Stephanie Reynolds, Pier Luigi Dragotti
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Place Field Mapping
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   INPUTS:
%       R: ratiometric time series whose size is [cells x samples]
%       spikes: non-negative derivative spike estimates
%       xy: x-y coordinates as a function of time
%       n: number of spatial bins
%       CaFs: sampling frequency of two-photon calcium imaging
%       trackFs: sampling frequency of the tracker
%       mode: 1 if 1d, 2 if 2d
%   OUTPUTS:
%       Rmap: rate map
%       occMap: occupancy map
%       placeMap: place field map
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Rmap, occMap, placeMap, nxy] = neuroSEE_PFmap(R,spikes,xy,phi,n,...
    CaFs,trackFs,mode)

thr = 30;
thrV = thr/trackFs;
nR = bsxfun(@rdivide, bsxfun(@minus, R, min(R,[],2)),(max(R,[],2)-min(R,[],2)));
nC = size(nR,1);

if CaFs>trackFs
    %Downsample Ca trace to tracking
    downSbins = round(linspace(1,length(R),length(xy)+1));
    nRd = zeros(size(nR,1),length(xy));
    % nRd: downsampled, normalised ratiometric time series
    for a=1:length(xy)
        nRd(:,a) = mean(spikes(:,downSbins(a):downSbins(:,a+1)),2);
    end
else
    %Downsample tracking to Ca trace
    downSbins = round(linspace(1,length(xy),length(R)+1));
    xyd = zeros(length(R),2);
    % xyd: downsampled tracking series
    for a=1:length(R)
        xyd(a,:) = mean(xy(downSbins(a):downSbins(a+1),:),1);
    end
    xy = xyd;
    nRd = nR;
end

%% Consider only samples when the mouse is active
active = getActiveSamples(xy,thrV);
xyA = xy(active,:);
nRa = nRd(:,active);
nRa = double(nRa);

% normalise tracking to unit square
x = xyA(:,1);
y = xyA(:,2);
x = (x-min(x))/(max(x)-min(x));
y = (y-min(y))/(max(y)-min(y));
nxy = [x,y];

if mode==1 
    theta = phi;
    %theta = atan(xy(:,2)./xy(:,1));
    % make sure theta has positive angles
    %{
    for k = 1:length(theta)
        if theta(k)<0
            theta(k) = theta(k) + 2*pi;
        end
    end
    %}
    thetaA = theta(active,:);
    % Normalise theta values to unit 1
    theta_n = (thetaA-min(thetaA))/(max(thetaA)-min(thetaA));
    % Obtain occupational map - create computational grid
    theta_grids = linspace(min(theta_n),max(theta_n)*1.0001,n+1); % 2d bins
    l = floor((theta_n(:,1)-theta_grids(1))/(theta_grids(2)-theta_grids(1)))+1; % discretising location
    
    % enters a 1 in a sparse matrix at each xy position
    occMap = full(sparse(1,l,1,1,n));
  
    % Obtain ratemaps
    Rmap = zeros(1,n,nC);
    placeMap = zeros(1,n,nC);
    for ri=1:nC
        % enters dR/R in a sparse matrix at each xy position
        Rmap(:,:,ri) = full(sparse(1,l,nRa(ri,:),1,n));
        pm = Rmap(:,:,ri)./occMap;
        pm(isnan(pm)) = 0;
        placeMap(:,:,ri) = pm;
    end
    
elseif mode==2

% Obtain occupational map - create computational grid
xh = linspace(min(x),max(x)*1.0001,n+1); % 2d bins
yh = linspace(min(y),max(y)*1.0001,n+1); % 2d bins
i = floor((nxy(:,1)-xh(1))/(xh(2)-xh(1)))+1; % discretising location
j = floor((nxy(:,2)-yh(1))/(yh(2)-yh(1)))+1;
occMap = full(sparse(i,j,1,n,n));

%% Obtain rate maps and place maps
Rmap = zeros(n,n,nC);
placeMap = zeros(n,n,nC);
    for ri=1:nC
        % enters dR/R in a sparse matrix at each xy position
        Rmap(:,:,ri) = full(sparse(i,j,nRa(ri,:),n,n));
        pm = Rmap(:,:,ri)./occMap;
        pm(isnan(pm)) = 0;
        placeMap(:,:,ri) = pm;
    end
end
end