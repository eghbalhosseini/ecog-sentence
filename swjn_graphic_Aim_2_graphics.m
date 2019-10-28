clear all
close all
data_path='/Users/eghbalhosseiniasl1/MyData/ecog-sentence';

%
save_path='/Users/eghbalhosseiniasl1/MyData/ECoG-sentence/crunched/';
analysis_path='/Users/eghbalhosseiniasl1/MyData/ECoG-sentence/analysis/swjn_Aim_2_graphics';
%
if 1
    fprintf('adding basic ecog tools to path \n');
    addpath('~/MyCodes/basic-ecog-tools/');
    addpath(genpath('~/MyCodes/ecog-sentence/'));
    addpath(genpath('~/MyCodes/basic-ecog-tools/activeBrain'));
    addpath(genpath('~/MyCodes/basic-ecog-tools/ecog-filters'));
    addpath(genpath('~/MyCodes/basic-ecog-tools/mex'));
end


%% mean change 
x=[1:16];
x_1=[1:16]-.1;
x_2=[1:16]+.1;
y_1=[.3+.3*rand(1,8),1.2+.3*rand(1,8)] ;
y_2=[1.2+.3*rand(1,8),.3+.3*rand(1,8)] ;
%y_2=mean(y_1)+[.2*rand(1,16)] ;
npts=4;
c_1=[1,1,1,1];
f=figure;
set(f,'position',[1000 1022 1303 316])
ax=axes('position',[.1,.1,.7,.4]);
a=stem(x_1,y_1,'DisplayName','pMI','Color',[0,0,0],'LineWidth',2);
ax.XTick=[1:16];
hold on; 
a=stem(x_2,y_2,'DisplayName','pMI','Color',[.5,.5,.5],'LineWidth',2);
ax.XAxis.TickLength=[0,0];
%thr=.25;
%plot(ax.XLim,ax.YLim*0+thr,'--','DisplayName','threshold','Color',[.5,.5,.5])

ax.XLabel.String='word position'
ax.YLabel.String='PMI'
ax.YLim=[0,1.6]
ax.YTick=[];
ax.Box='off';
ax.XLim=[0,17];
%y_nan=y;
%y_nan(y_nan<thr)=0;
%legend('show','Location','northwest')
%text(ax.XLim(1),ax.YLim(2)+.1,'A','HorizontalAlignment','right','FontSize',12,'FontWeight','bold')
ax=axes('position',[.1,.55,.7,.35]);

ax.YTick=[];
ax.XTick=[];
y_auto_1=conv(y_1,c_1);
y_auto_2=conv(y_2,c_1);
a=plot(x,y_auto_1(x),'o-','DisplayName','constant kernel','color',[1,0,.5],'linewidth',2);
hold on;
a=plot(x,y_auto_2(x),'o-','DisplayName','constant kernel','color',[1,.3,.5],'linewidth',2);
%text(ax.XLim(1),ax.YLim(2)+.3,'B','HorizontalAlignment','right','FontSize',12,'FontWeight','bold')
ax.XLim=[0,17];
ax.XAxis.Visible='off';
box off
ax.YLabel.String='\gamma';
ax.YTick=[];

a=plot(x,cumsum(y,2),'o--','DisplayName','no kernel','color',[.5,.5,.5],'linewidth',2);
ax.XTick=x;
text(ax.XLim(1),ax.YLim(2)+.4,'D','HorizontalAlignment','right','FontSize',12,'FontWeight','bold')
legend('show','Location','northwest')
arr.MaxHeadSize=0;
        if ~exist(strcat(analysis_path))
            mkdir(strcat(analysis_path))
        end 

% 

%print(f, '-djpeg', strcat(analysis_path,'/norm_sum.jpeg'));

%% 
