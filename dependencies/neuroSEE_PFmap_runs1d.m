function [Rmap, occMap, placeMap, nxy] = neuroSEE_PFmap_runs1d(R,spikes,xy, theta,n,...
    CaFs,trackFs,mode)

thr = 0.02;
thrV = thr/trackFs;
nspikes = bsxfun(@rdivide, bsxfun(@minus, spikes, min(spikes,[],2)),(max(spikes,[],2)-min(spikes,[],2))); %normalise
nC = size(R,1); %number of cells
nRd = nspikes;

%% Consider only samples when the mouse is active
active = getActiveSamples(xy,thrV);
xyA = xy(active,:);
thetaA = theta(active,:);
nRa = nRd(:,active);
nRa = double(nRa);

% normalise tracking to unit square
x = xyA(:,1);
y = xyA(:,2);
x = (x-min(x))/(max(x)-min(x));
y = (y-min(y))/(max(y)-min(y));
nxy = [x,y];

if mode==1 
    n = 360;
    
    % enters a 1 in a sparse matrix at each theta position
    occMap = full(sparse(1,thetaA,1,1,n));
  
    % Obtain ratemaps
    Rmap = zeros(nC,n);
    placeMap = zeros(nC,n);
    for ri=1:nC
        % enters dR/R in a sparse matrix at each xy position
        Rmap(ri,:) = full(sparse(1,thetaA,nRa(ri,:),1,n));
        pm = Rmap(ri,:)./occMap;
        pm(isnan(pm)) = 0;
        placeMap(ri,:) = pm;
    end
    
elseif mode==2

% Obtain occupational map - create computational grid
xh = linspace(min(x),max(x)*1.0001,n+1); % 2d bins
yh = linspace(min(y),max(y)*1.0001,n+1); % 2d bins
i = floor((nxy(:,1)-xh(1))/(xh(2)-xh(1)))+1; % discretising location
j = floor((nxy(:,2)-yh(1))/(yh(2)-yh(1)))+1;
occMap = full(sparse(i,j,1,n,n));

Rmap = zeros(n,n,nC);
placeMap = zeros(n,n,nC);
    for ri=1:nC
        % enters dR/R in a sparse matrix at each xy position
        Rmap(:,:,ri) = full(sparse(i,j,nRa(ri,:),n,n));
        pm = Rmap(:,:,ri)./occMap;
        pm(isnan(pm)) = 0;
        placeMap(:,:,ri) = pm;
    end

elseif mode==3
    % Normalise theta values to unit 1
    theta_n = (thetaA-min(thetaA))/(max(thetaA)-min(thetaA));
    % Obtain occupational map - create computational grid
    theta_grids = linspace(min(theta_n),max(theta_n)*1.0001,n+1);
    l = floor((theta_n(:,1)-theta_grids(1))/(theta_grids(2)-theta_grids(1)))+1; % discretising location
    
    % enters a 1 in a sparse matrix at each theta position
    occMap = full(sparse(1,l,1,1,n));
  
    % Obtain ratemaps
    Rmap = zeros(nC,n);
    placeMap = zeros(nC,n);
    for ri=1:nC
        % enters dR/R in a sparse matrix at each xy position
        Rmap(ri,:) = full(sparse(1,l,nRa(ri,:),1,n));
        pm = Rmap(ri,:)./occMap;
        pm(isnan(pm)) = 0;
        placeMap(ri,:) = pm;
    end
end
end