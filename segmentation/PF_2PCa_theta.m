function placeMap = PF_2PCa_theta(day,green_ch,red_ch,track_xy,track_theta,DECONV,PLOT,thr)
    %% Set up
    % manage inputs
    if nargin<5; DECONV=0; end
    if nargin<6; PLOT=1; end
    if nargin<7; thr=30; end
    
    n = 60; % no spatial bins

    CaFs = 30.92; % ca trace freq
    trackFs = 15; % tracking freq
    thrV = thr/trackFs; % to be discusse as the tracking is not in cm/s yet

    %% Computing dR/R from the green and red channel Ca time series
    C = size(green_ch,1); %number of cells
    Tf = size(green_ch,2); %number of frames
    t = [1:Tf]*1/30;
    T = max(t);
    sm=7;

    % first smooth
    for c=1:C
        green_ch(c,:) = smooth(green_ch(c,:),sm);
        red_ch(c,:) = smooth(red_ch(c,:),sm);
    end
    R_ = green_ch./red_ch;

    R0 = mode(R_');
    R0_ = repmat(R0',[1 Tf]);
    dRonR = (R_-R0_)./R0_;
    
    R = dRonR;
    
    % Normalise the dR/R Ca traces
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
    downSbins = round(linspace(1,length(R),length(track_theta)+1));
    nRd = zeros(nC,length(track_theta));
    
    for a=1:length(track_theta)
        nRd(:,a) = mean(Fdot(:,downSbins(a):downSbins(:,a+1)),2);
    end
    
    % Consider only samples when the mouse is active
    active = getActiveSamples(track_xy,thrV);
    thetaA = track_theta(active,:);
    nRa = nRd(:,active);
    
    xyA = track_xy(active,:);
    nRa = nRd(:,active);
    % normalise tracking to unit square
    x = xyA(:,1);
    y = xyA(:,2);
    x = (x-min(x))/(max(x)-min(x));
    y = (y-min(y))/(max(y)-min(y));
    nxy = [x,y];
    
    % Normalise theta values to unit 1
    theta = (thetaA-min(thetaA))/(max(thetaA)-min(thetaA));
   
    % Obtain occupational map - create computational grid
    theta_grids = linspace(min(theta),max(theta)*1.0001,n+1); % 2d bins
    i = floor((theta(:,1)-theta_grids(1))/(theta_grids(2)-theta_grids(1)))+1; % discretising location

    % enters a 1 in a sparse matrix at each xy position
    occMap = full(sparse(1,i,1,1,n));
  
    % Obtain ratemaps
    Rmap = zeros(1,n,nC);
    placeMap = zeros(1,n,nC);
    for ri=1:nC
        % enters dR/R in a sparse matrix at each xy position
        Rmap(:,:,ri) = full(sparse(1,i,nRa(ri,:),1,n));
        pm = Rmap(:,:,ri)./occMap;
        pm(isnan(pm)) = 0;
        placeMap(:,:,ri) = pm;
    end

    %% Plotting
    if PLOT
        figure()
        plot(nxy(:,1),nxy(:,2))
        title('Tracking')

        figure('Units','normalized','OuterPosition',[0 0 0.8 0.2]) % occupational map
        imagesc(occMap)
        title(['Day ',num2str(day),' - Occupational Map'])
        c = colorbar();
        c.Label.String = 'cumulative samples';
        
        % rate maps
        figure('Units','normalized','OuterPosition',[0 0 0.8 0.8])
        title(['Day ',num2str(day),' - Rate Maps']);
        for ip=1:nC
            subplot(nC,1,ip)
            imagesc(Rmap(:,:,ip))
            c = colorbar();
            c.Label.String = 'cumulative dR/R';
            title(['Cell ',num2str(ip)])
        end
        subtitle(['Day ',num2str(day),' - Rate Maps']);
        
        % place maps
        figure('Units','normalized','OuterPosition',[0 0 0.8 0.8])
        title(['Day ',num2str(day),' - Place Maps']);
        for ip=1:nC
            subplot(nC,1,ip)
            %subplot(floor(nC/2),2+rem(nC,2),ip)
            imagesc(placeMap(:,:,ip))
            c = colorbar();
            c.Label.String = 'dR/R per sample';
            title(['Cell ',num2str(ip)])
        end
        subtitle(['Day ',num2str(day),' - Place Maps']);
       
    end

end

%% Extra functions

function active = getActiveSamples(tracking,thrV)
    vx = diff(tracking(:,1));
    vy = diff(tracking(:,2));
    active = sqrt(vx.^2+vy.^2)>thrV;
end