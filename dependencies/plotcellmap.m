%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NeuroSEE: An automated Neuronal Source Extraction
%             and Exploration toolbox
%   
%   Author: Seigfred Prado   
%   Supervisor: Simon Schultz
%   Acknowledgment: Stephanie Reynolds, Pier Luigi Dragotti
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotcellmap(summaryImage, cellIDs, cellMasks, N, deltaR, offsetframes)

% summaryImage: the image used for initialisation; could be correlation or
% mean image
% cellIDs: a vector that contains the IDs of the cells to be mapped
% N: number of frames

numCells = length(cellIDs);
figure;
opts.m = numCells;
opts.n = 4;
opts.p = 1;
opts.subplot = 1;
masks = cellMasks;

%T          = 1/8;
t          = 1:N/offsetframes;
line_width = 2.5;

rId = zeros(size(masks,3),2);
cId = zeros(size(masks,3),2);
for i = 1:size(masks,3)
    [row,col]=find(masks(:,:,i));
    rId(i,1) = min(row);
    rId(i,2) = max(row);
    cId(i,1) = min(col);
    cId(i,2) = max(col);
end

r_scale = zeros(size(masks,3),2);
c_scale = zeros(size(masks,3),2);
for i = 1:size(masks,3)
   x = round((50-(rId(i,2)-rId(i,1)))/2);
   y = round((50-(cId(i,2)-cId(i,1)))/2); 
   r_scale(i,1)=rId(i,1)-x;
   if r_scale(i,1)<0
       r_scale(i,1)=0;
       r_scale(i,2)=49;
   else
       r_scale(i,2)=rId(i,2)+(50-(rId(i,2)-rId(i,1)))-x-1;
   end
   c_scale(i,1)=cId(i,1)-y;
   if c_scale(i,1)<0
       c_scale(i,1)=0;
       c_scale(i,2)=49;
   else
       c_scale(i,2)=cId(i,2)+(50-(cId(i,2)-cId(i,1)))-y-1;
   end
   
   if r_scale(i,2)>size(masks,1)
       r_scale(i,2) = size(masks,1);
       r_scale(i,1) = size(masks,1)-49;
   end

   if c_scale(i,2)>size(masks,2)
       c_scale(i,2) = size(masks,2);
       c_scale(i,1) = size(masks,2)-49;
   end
end

for i=1:numCells
    if i==1
        opts.p=1;
        ii = 2:4;
    else
        opts.p = opts.p + 4;
        ii = opts.p+1:opts.p+3;
    end
    plotContoursOnSummaryImage(summaryImage, masks(:,:,cellIDs(i)),opts);
    axis([c_scale(cellIDs(i),1),c_scale(cellIDs(i),2),r_scale(cellIDs(i),1),r_scale(cellIDs(i),2)])
    subplot(numCells,4,ii); plot(t,deltaR(cellIDs(i),:),'k','LineWidth', line_width);
    set(gca,'Visible','off')
    set(gca,'color','w')
    xlim([0 3720])
    ylim([-1 3])
end
end
