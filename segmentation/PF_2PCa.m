function [placeMap, occMap, Rmap, nxy] = PF_2PCa(day,green,red,xy,DECONV,PLOT,thr)
    %% Set up
    % manage inputs
    if nargin<4; DECONV=0; end
    if nargin<5; PLOT=0; end
    if nargin<6; thr=30; end
    % set directories and parameters
    homeDir = cd;
    resDir = [homeDir,'/m40'];
    addpath(genpath([homeDir,'/CaGPPF'])); % add path for GP estimation
    n = 10; % no spatial bins
    ng = 50; % no spatial bins for GP estimation
    CaFs = 30.92; % ca trace freq
    trackFs = 15; % tracking freq
    thrV = thr/trackFs; % to be discusse as the tracking is not in cm/s yet

    %% Work
    C = size(green,1);
    Tl = size(green,2);
    t = [1:Tl]*1/30;
    T = max(t);
    sm=7;

    % first smooth
    for c=1:C
        green(c,:) = smooth(green(c,:),sm);
        red(c,:) = smooth(red(c,:),sm);
    end
    R_ = green./red;

    R0 = mode(R_');
    R0_ = repmat(R0',[1 Tl]);
    dRonR = (R_-R0_)./R0_;
    
    R = dRonR;
 
    nR = bsxfun(@rdivide, bsxfun(@minus, R, min(R,[],2)),(max(R,[],2)-min(R,[],2)));
    nC = size(nR,1);
    % non-negative derivative of Ca signal
    mAlength = 5;
    thr = 0.04;
    Fdot = zeros(size(nR));
    if ~DECONV
        for ri=1:nC
            nRs = conv(nR(ri,:),ones(mAlength,1)/mAlength,'same');
            temp = diff(nRs);
            temp = [temp,temp(end)];
            temp(temp<thr) = 0;
            Fdot(ri,:) = temp;
        end
    else
        % non negative deconvolution (OASIS-AR1, Paninski)
        deconvR = zeros(size(R));
        g1 = 0.95; lam1 = 2.4;
        for ri=1:nC
            [~,deconvR(ri,:)] = oasisAR1(nR(ri,:),g1,lam1);
        end
        Fdot = deconvR;
    end
    % downsample Ca trace to tracking
    downSbins = round(linspace(1,length(R),length(xy)+1));
    nRd = zeros(size(nR,1),length(xy));
    for a=1:length(xy)
        nRd(:,a) = mean(Fdot(:,downSbins(a):downSbins(:,a+1)),2);
    end
    % Consider only samples when the mouse is active
    active = getActiveSamples(xy,thrV);
    xyA = xy(active,:);
    nRa = nRd(:,active);
    % normalise tracking to unit square
    x = xyA(:,1);
    y = xyA(:,2);
    x = (x-min(x))/(max(x)-min(x));
    y = (y-min(y))/(max(y)-min(y));
    nxy = [x,y];
    % Obtain occupational map - create computational grid
    xh = linspace(min(x),max(x)*1.0001,n+1); % 2d bins
    yh = linspace(min(y),max(y)*1.0001,n+1); % 2d bins
    i = floor((nxy(:,1)-xh(1))/(xh(2)-xh(1)))+1; % discretising location
    j = floor((nxy(:,2)-yh(1))/(yh(2)-yh(1)))+1;
    % enters a 1 in a sparse matrix at each xy position
    occMap = full(sparse(i,j,1,n,n));
    mask = occMap>0;
    % Obtain ratemaps
    Rmap = zeros(n,n,nC);
    placeMap = zeros(n,n,nC);
    for ri=1:nC
        % enters dR/R in a sparse matrix at each xy position
        Rmap(:,:,ri) = full(sparse(i,j,nRa(ri,:),n,n));
        pm = Rmap(:,:,ri)./occMap;
        pm(isnan(pm)) = 0;
        placeMap(:,:,ri) = pm;
    end

    %% GP place map inference

    %     kernel = 0; % SE kernel
    %     GPest = zeros(ng,ng,nC);
    %     for cid=1:1
    %         [pf,hyp2,counts,ll] = GPPF_CaData(nRa(cid,:)',nxy,ng,ng,kernel,[]);
    %         GPest(:,:,cid) = pf.mtuning*trackFs;
    %     end

    %% Plotting
    if PLOT
        figure()
        plot(nxy(:,1),nxy(:,2))
        title('Tracking')

        figure % occupational map
        imagesc(occMap)
        title(['Day ',num2str(day),' - Occupational Map'])
        c = colorbar();
        c.Label.String = 'cumulative samples';
        % rate maps
        figure('units','normalized','outerposition',[0 0 0.8 0.8])
        title(['Day ',num2str(day),' - Rate Maps']);
        for ip=1:nC
            subplot(2+rem(nC,2),floor(nC/2),ip)
            imagesc(Rmap(:,:,ip))
            c = colorbar();
            c.Label.String = 'cumulative dF/F';
            title(['cell ',num2str(ip)])
        end
        subtitle(['Day ',num2str(day),' - Rate Maps']);
        
        % place maps
        figure('units','normalized','outerposition',[0 0 0.8 0.8])
        title(['Day ',num2str(day),' - Place Maps']);
        for ip=1:nC
            subplot(2+rem(nC,2),floor(nC/2),ip)
            imagesc(Rmap(:,:,ip)./occMap)
            c = colorbar();
            c.Label.String = 'dF/F per sample';
            title(['cell ',num2str(ip)])
        end
        subtitle(['Day ',num2str(day),' - Place Maps']);
       

        %     figure
        %     imagesc(GPest(:,:,cid))
    end

end

%% Extra functions

function active = getActiveSamples(tracking,thrV)
    vx = diff(tracking(:,1));
    vy = diff(tracking(:,2));
    active = sqrt(vx.^2+vy.^2)>thrV;
end