clear all
close all
data_path='/Users/eghbalhosseiniasl1/MyData/ecog-sentence';

analysis_path='/Users/eghbalhosseiniasl1/MyData/ECoG-sentence/analysis/swjn_graphic_incremental_mutual_info';
%
if 1
    fprintf('adding basic ecog tools to path \n');
    addpath('~/MyCodes/basic-ecog-tools/');
    addpath(genpath('~/MyCodes/ecog-sentence/'));
    addpath(genpath('~/MyCodes/basic-ecog-tools/activeBrain'));
    addpath(genpath('~/MyCodes/basic-ecog-tools/ecog-filters'));
    addpath(genpath('~/MyCodes/basic-ecog-tools/mex'));
end


%%
x=2:8;
y=rand(1,7) ;
close all 
f=figure;
aspect_ration=9.32./4.13;
y_dim=250;
set(f,'position',[591 455 aspect_ration*y_dim y_dim]);
ax=axes('position',[.05,.2,.3,.4]);
ax.FontSize=12;
ax.FontWeight='bold'
hold on
a=stem(x,y,'DisplayName','pMI','Color',[0,0,0],'LineWidth',2)
ax.XTick=x;
thr=.4;
plot(ax.XLim,ax.YLim*0+thr,'--','DisplayName','threshold','Color',[.5,.5,.5])
ax.XTick=x;
ax.XLabel.String='word position'
ax.YLabel.String='amplitude'
ax.YLim=[0,1.5];
ax.XLim=[1,9];
y_nan=y;
y_nan(y_nan<thr)=0;
legend('show','Location','northwest');
ax.YTick=[];
%text(ax.XLim(1)-.5,ax.YLim(2)+.15,'A','HorizontalAlignment','right','FontSize',14,'FontWeight','bold')
ax=axes('position',[.45,.2,.4,.4]);
hold on
a=plot(x,cumsum(y),'o-','DisplayName','without threshold','color',[0,0,0],'linewidth',2)
ax.XTick=x;
ax.YLabel.String='gamma power'
ax.XLabel.String='word position'
ax.YTick=[];
ax.XLim=[1,9];
hold on
a=plot(x,nancumsum(y_nan,2),'o-','DisplayName','with threshold','color',[.5,.5,.5],'linewidth',2)
ax.XTick=x;
%text(ax.XLim(1)-.5,ax.YLim(2)+.5,'B','HorizontalAlignment','right','FontSize',14,'FontWeight','bold')
legend('show','Location','northwest')
ax.FontSize=12;
ax.FontWeight='bold'
arr.MaxHeadSize=0;
        if ~exist(strcat(analysis_path))
            mkdir(strcat(analysis_path))
        end 

% 
print(f, '-djpeg', strcat(analysis_path,'/pmi_threshold_presentation.jpeg'));
%% graph for threshold 
thr=exp(-.5*(1:1:10));
f=figure;
aspect_ration=9.32./4.13;
y_dim=250;
set(f,'position',[591 455 aspect_ration*y_dim y_dim]);
ax=axes('position',[.05,.2,.3,.4]);
ax.FontSize=12;
ax.FontWeight='bold'
hold on
a=bar(thr);
a.FaceColor=[1,.2,.2]
a.LineWidth=2
ax.YTick=[];
ax.YLabel.String='R^2'
ax.XTick=[1,5]
ax.XTickLabel={'th=0','th>0'};
% 

thr=gampdf([1:1:10],2.5,2);
ax=axes('position',[.45,.2,.3,.4]);
ax.FontSize=12;
ax.FontWeight='bold'
hold on
a=bar(thr);
a.FaceColor=[1,.2,.2]
a.LineWidth=2
ax.YTick=[];
ax.YLabel.String='R^2'
ax.XTick=[1,5]
ax.XTickLabel={'th=0','th>0'};
print(f, '-djpeg', strcat(analysis_path,'/pmi_threshold_Rsqured.jpeg'));

%% 
[X,Y] = meshgrid(cosd(theta)',sind(theta)');
f=figure;
ax=axes('position',[.1,.1,.8,.8]);
hold on
arr=arrayfun(@(x,y) quiver3(0*x,0*x,0*x,x,y,x*0,0,'k'), X(1,1:end),transpose(Y(1:end,1)))
arrayfun(@(x) set(x,'linewidth',1),arr)
arr=quiver3(0,0,0,0,0,1,0,'k');
arrayfun(@(x) set(x,'linewidth',1),arr);
set(ax,'View',[151.8282 19.1377])
ax.YAxis.Visible='off';
ax.XAxis.Visible='off';
ax.ZAxis.Visible='off';
% 
U=rand(1,3);
U=U./norm(U);
V=rand(1,3);
V=V./norm(V);
U=1*[.3,.6,.8];
V=1*[.6,.3,.5];
V=V./norm(V);
U=U./norm(U);
W=U*dot(U,V)/(norm(U)*norm(V));
arr=quiver3(0,0,0,W(1),W(2),W(3),0,'Color',[.2,.2,.2]);
text(U(1)/2-.02,U(2)/2,U(3)/2.5,'|V_{w_j}|','HorizontalAlignment','left','fontsize',16)
arr.MaxHeadSize=0;
arrayfun(@(x) set(x,'linewidth',3),arr)
arr=quiver3(0,0,0,U(1),U(2),U(3),0,'r');
arr.MaxHeadSize=.5;
text(U(1)-.05,U(2),U(3),'V_{w_i}','HorizontalAlignment','left','fontsize',16)
arrayfun(@(x) set(x,'linewidth',2),arr)
arr=quiver3(0,0,0,V(1),V(2),V(3),0,'b');
text(V(1)/2-.05,V(2)/2,V(3)/2,'V_{w_j}','HorizontalAlignment','left','fontsize',16)
arr.MaxHeadSize=.8;
arrayfun(@(x) set(x,'linewidth',2),arr)

W=W-V;
arr=quiver3(V(1),V(2),V(3),W(1),W(2),W(3),0,'Color',[.2,.2,.2],'LineStyle','--');
arr.MaxHeadSize=0;
        if ~exist(strcat(analysis_path))
            mkdir(strcat(analysis_path))
        end 
print(f, '-djpeg', strcat(analysis_path,'/norm_sum.jpeg'));

%% 
