%Place Field Mapping


green = [Green(:,1:9780) Green(:,13501:end)];
red = [Red(:,1:9780) Red(:,13501:end)];
xy = [xy_summary(1:5403,:); xy_summary(7205:end,:)];

DECONV=0;
thr=30;

n = 40; % no spatial bins
CaFs = 30.92; % ca trace freq
trackFs = 15; % tracking freq
thrV = thr/trackFs;

%% dR/R
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
R = (R_-R0_)./R0_;
dRonR = R;
dRonR_max = max(dRonR);
maxmax = max(dRonR_max);
dRonR_min = min(dRonR);


%% Normalise Ca traces
nR = bsxfun(@rdivide, bsxfun(@minus, R, min(R,[],2)),(max(R,[],2)-min(R,[],2)));
nC = size(nR,1);


%% non-negative derivative of Ca signal
mAlength = 5;
thr = 0.009;
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

%% downsample Ca trace to tracking
downSbins = round(linspace(1,length(R),length(xy)+1));
nRd = zeros(size(nR,1),length(xy));
for a=1:length(xy)
    nRd(:,a) = mean(Fdot(:,downSbins(a):downSbins(:,a+1)),2);
end

%% Consider only samples when the mouse is active
active = getActiveSamples(xy,thrV);
xyA = xy(active,:);
nRa = nRd(:,active);

%% normalise tracking to unit square
x = xyA(:,1);
y = xyA(:,2);
x = (x-min(x))/(max(x)-min(x));
y = (y-min(y))/(max(y)-min(y));
nxy = [x,y];

%% Obtain occupational map - create computational grid
xh = linspace(min(x),max(x)*1.0001,n+1); % 2d bins
yh = linspace(min(y),max(y)*1.0001,n+1); % 2d bins
i = floor((nxy(:,1)-xh(1))/(xh(2)-xh(1)))+1; % discretising location
j = floor((nxy(:,2)-yh(1))/(yh(2)-yh(1)))+1;
occMap = full(sparse(i,j,1,n,n));

%% Obtain ratemaps
Rmap = zeros(n,n,nC);
placeMap = zeros(n,n,nC);
for ri=1:nC
    % enters dR/R in a sparse matrix at each xy position
    Rmap(:,:,ri) = full(sparse(i,j,nRa(ri,:),n,n));
    pm = Rmap(:,:,ri)./occMap;
    pm(isnan(pm)) = 0;
    placeMap(:,:,ri) = pm;
end

%
%% Plot behavioural track
figure()
plot(nxy(:,1),nxy(:,2))
title('Tracking')

%% Plot occupancy map
figure 
imagesc(occMap)
title('Occupancy Map')
c = colorbar();
c.Label.String = 'cumulative samples';

%% Plot rate map
figure('units','normalized','outerposition',[0 0 0.8 0.8])
title('Rate Maps');
for ip=1:nC
    %subplot(2+rem(nC,2),floor(nC/2),ip)
    subplot(5,2,ip)
    imagesc(Rmap(:,:,ip))
    c = colorbar();
    c.Label.String = 'cumulative dF/F';
    title(['cell ',num2str(ip)])
end
subtitle('Rate Maps');

%% Plot place maps

figure('units','normalized','outerposition',[0 0 0.8 0.8])
title('Place Maps');
for ip=1:nC
    %subplot(2+rem(nC,2),floor(nC/2),ip)
    subplot(5,2,ip)
    imagesc(Rmap(:,:,ip)./occMap)
    c = colorbar();
    c.Label.String = 'dF/F per sample';
    title(['cell ',num2str(ip)])
end
subtitle('Place Maps');



%% Plot summary
figure('units','normalized','outerposition',[0 0 0.7 1])
subplot(4,3,1)
imagesc(occMap)
c = colorbar();
c.Label.String = 'cumulative samples';
title('Occupancy Map')

subplot(4,3,2)
imagesc(Rmap(:,:,4))
c = colorbar();
c.Label.String = 'cumulative dR/R';
title('Rate Map: Cell 4')

subplot(4,3,3)
imagesc(placeMap(:,:,4))
c = colorbar();
c.Label.String = 'dR/R per sample';
title('Place Map: Cell 4')

subplot(4,3,4)
plot(nxy(1:300,1),nxy(1:300,2),'Color',[0,0.7,0.9])
title('Behavioural Track')
hold on; plot(nxy(577:end,1),nxy(577:end,2),'Color',[0,0.7,0.9]);
hold on; plot(nxy(301:438,1),nxy(301:438,2),'Color',[0,0.7,0.9]);
hold on; plot(nxy(439:576,1),nxy(439:576,2),'Color',[0,0.7,0.9]);

subplot(4,3,5)
bar(1:100,histocc,'FaceColor',[0,0.7,0.9]); xlim([1 100]);
title('Histogram of Bin Occupancy')
xlabel('Spatial Bins')

subplot(4,3,6)
histogram(speeds_active,'FaceColor',[0,0.7,0.9]); xlim([0 10])
title('Histogram of Speeds')
xlabel('Speed (cm/s)')

subplot(4,3,7)
plot(t(1:2340),R(4,1:2340),'k');
xlabel('Time (s)');
ylabel('dR/R');
title('Ca Traces: Day 1, Cell 4');

subplot(4,3,8)
plot(t(6061:9780)-202,R(4,6061:9780),'k');
xlabel('Time (s)');
ylabel('dR/R');
title('Ca Traces: Day 2, Cell 4');

subplot(4,3,9)
plot(t(13501:17220)-450,R(4,13501:17220),'k');
xlabel('Time (s)');
ylabel('dR/R');
title('Ca Traces: Day 3, Cell 4');

subplot(4,3,10)
plot(t(1:2340),Fdot1(4,1:2340),'k');
xlabel('Time (s)');
%ylabel('dR/R');
title('Spikes: Day 1, Cell 4');

subplot(4,3,11)
plot(t(6061:9780)-202,Fdot2(4,6061:9780),'k');
xlabel('Time (s)');
%ylabel('dR/R');
title('Spike: Day 2, Cell 4');

subplot(4,3,12)
plot(t(13501:17220)-450,Fdot3(4,13501:17220),'k');
xlabel('Time (s)');
%ylabel('dR/R');
title('Spikes: Day 3, Cell 4');

%
%set(gcf,'papertype','A4');
%print(fig,'-dpdf');
%
%