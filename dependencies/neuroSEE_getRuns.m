function [runs] = neuroSEE_getRuns(R, spikes, xy, theta)

%downsample tracking data
downSbins_xy = round(linspace(1,length(xy),length(R)+1));
downSbins_theta = round(linspace(1,length(xy),length(R)));

for a=1:length(R)
    xyd(a,:) = mean(xy(downSbins_xy(a):downSbins_xy(a+1),:),1);
end

theta = theta(downSbins_theta,1); 
xy = xyd;

theta = floor(theta); %discretise the theta positions
theta(theta==0)=360; %replace all 0s with 360

% Translate the theta values from 1 to 360 starting from the initial theta
% value:
theta_test=zeros(length(theta),1);
ind = [theta(1):-1:1 360:-1:theta(1)+1]'; 
for i=1:length(theta)
theta_test(i,1) = find(ind==theta(i));
end
theta = theta_test;

diff_theta = abs(diff(theta)); %difference of the consecutive theta values
rep_ind = find(diff_theta>100); %find the indices of where the set of theta values repeats

%translate xy values into normalised ones according to the translated theta
%values
xy(:,1) = cosd(theta);
xy(:,2) = sind(theta);

runs = cell(4,length(rep_ind)+1);

for j = 1:length(rep_ind)+1
    if j == 1
        runs{1,j} = R(:,1:rep_ind(j));
        %runs{2,j} = spikes(:,1:rep_ind(j));
        runs{2,j} = neuroSEE_nnd(runs{1,j});
        runs{3,j} = xy(1:rep_ind(j),:);
        runs{4,j} = theta(1:rep_ind(j),:);
        r = rep_ind(j)+1;
    elseif j == length(rep_ind)+1
        runs{1,j} = R(:,r:size(R,2));
        %runs{2,j} = spikes(:,r:size(spikes,2));
        runs{2,j} = neuroSEE_nnd(runs{1,j});
        runs{3,j} = xy(r:size(xy,1),:);
        runs{4,j} = theta(r:size(theta,1),:);
    else
        runs{1,j} = R(:,r:rep_ind(j));
        %runs{2,j} = spikes(:,r:rep_ind(j));
        runs{2,j} = neuroSEE_nnd(runs{1,j});
        runs{3,j} = xy(r:rep_ind(j),:);
        runs{4,j} = theta(r:rep_ind(j),:);
        r = rep_ind(j)+1;
    end
end

indices = [];
for k = 1:size(runs,2)
    if size(runs{1,k},2) < 150
        runs{1,k} = [];
        runs{2,k} = [];
        runs{3,k} = [];
        runs{4,k} = [];
    else 
        indices =  [indices k];
    end
end
runs = runs(:,indices);
end