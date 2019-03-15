function plotCa2(R,gain)

if nargin<2
    gain=1.0;
end

C = size(R,1);
Tl = size(R,2);
t = [1:Tl]*1/30;
T = max(t);

dRonR =R;
dRonR_max = max(dRonR,[],2);
maxmax = max(dRonR_max);
dRonR_min = min(dRonR,[],2);

dy = 1.0/(C+1); y0=[1:C]*dy; % set axes positions

%range_dFonF = max_dFonF-min_dFonF;
%gain = 1.0; % turn this up to reduce whitespace between traces

figure;
clf
for c=1:C
   %subplot(C,1,c);
   plot(t,R(c,:)*gain+c,'k');
   axis([0 max(t) dRonR_min(c) maxmax]);
   hold on
end
hold off
axis([0 T 0 C+maxmax*gain]);
set(gca,'box','off');
set(gca,'visible','off');

% now plot scale bar
bar_T = 10; % sec
bar_dFonF = 1.0; % dFonF
axes('position',[0.7 0.05 bar_T/T*0.775 bar_dFonF*gain/(C+maxmax*gain)*0.815]);
plot([0 1],[0 0],'k','linewidth',2);
hold on
plot([0 0],[0 01],'k','linewidth',2);
hold off
text(0.5,0.05,'1.0 \Delta R/R');
text(0.5,0,'10 sec');
set(gca,'box','off');
set(gca,'visible','off');

