clear all
close all
data_path='/Users/eghbalhosseiniasl1/MyData/ecog-sentence';

%
save_path='/Users/eghbalhosseiniasl1/MyData/ECoG-sentence/crunched/';
analysis_path='/Users/eghbalhosseiniasl1/MyData/ECoG-sentence/analysis/swjn_graphic_auto_regressive_process';
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
theta=[0:15:90];
x=cosd(theta)
y=sind(theta)
z=theta*0+0;
 
[X,Y] = meshgrid(cosd(theta)',sind(theta)');
f=figure;
ax=axes('position',[.1,.1,.8,.8]);
hold on
arr=arrayfun(@(x,y) quiver3(0*x,0*x,0*x,x,y,x*0,0,'color',[.4,.4,.4]), X(1,1:end),transpose(Y(1:end,1)))
arrayfun(@(x) set(x,'linewidth',2),arr)
arr=quiver3(0,0,0,0,0,1,0,'color',[.4,.4,.4]);
arrayfun(@(x) set(x,'linewidth',2),arr);
set(ax,'View',[151.8282 19.1377])

% 

npts = 8;
t = linspace(0,2*pi,npts);
z = (sin(t/3)+.1*t)/1.5+.05;
x=0.1*t;
y=1-0.1613*t;
xyz = [x; y; z];
plot3(xyz(1,:),xyz(2,:),xyz(3,:),'o','LineWidth',2,'Color',[.2,.2,.2]);
text(xyz(1,:),xyz(2,:),xyz(3,:),[repmat('V_{w_',npts,1), num2str((1:npts)'),repmat('}',npts,1)],...
    'VerticalAlignment','bottom','fontsize',16,'HorizontalAlignment','center')
hold on
fnplt(cscvn(xyz),'k-',2);

ax.YAxis.Visible='off';
ax.XAxis.Visible='off';
ax.ZAxis.Visible='off';


arrayfun(@(x,y,z) quiver3(0,0,0,x,y,z,0,'Color',[1,0,0],'LineWidth',2,'MaxHeadSize',.5),x,y,z)
% 

if ~exist(strcat(analysis_path))
            mkdir(strcat(analysis_path))
end 
print(f, '-djpeg', strcat(analysis_path,'/vector_sum.jpeg'));
%% 
x=2:8;
y=rand(1,7) ;
npts=3;
c_1=[1,1,1];
c_2=exp(-.3*[1,2,3]);
f=figure;
set(f,'position',[595   595   778   210])
ax=axes('position',[.05,.3,.15,.5]);
hold on
a=stem(x,y,'DisplayName','pMI','Color',[0,0,0],'LineWidth',2)
ax.XTick=x;
%thr=.25;
%plot(ax.XLim,ax.YLim*0+thr,'--','DisplayName','threshold','Color',[.5,.5,.5])
ax.XLabel.String='word position'
ax.YLabel.String='PMI'
ax.YLim=[0,1.5]
ax.XLim=[1,9];
ax.YTick=[];
ax.FontSize=16

text(ax.XLim(1),ax.YLim(2)+.2,'A','HorizontalAlignment','right','FontSize',20,'FontWeight','bold')
ax=axes('position',[.25,.65,.05,.1]);
a=stem(1:length(c_1),c_1,'Color',[1,0,.5],'LineWidth',2)
box off
ax.YTick=[];
ax.XTick=[];
text(ax.XLim(1),ax.YLim(2)+.6,'B','HorizontalAlignment','right','FontSize',20,'FontWeight','bold')
ax=axes('position',[.25,.3,.05,.15]);
a=stem(1:length(c_2),c_2,'Color',[1,.5,0],'LineWidth',2)
box off
ax.XTickLabel={'c_1','c_2','c_3'};
ax.YTick=[];
text(ax.XLim(1),ax.YLim(2)+.5,'C','HorizontalAlignment','right','FontSize',20,'FontWeight','bold')
ax.FontSize=16

ax=axes('position',[.35,.3,.2,.5]);
hold on
y_auto_1=conv(y,c_1);
y_auto_2=conv(y,c_2);
a=plot(x,y_auto_1(x),'o-','DisplayName','constant','color',[1,0,.5],'linewidth',2)
a=plot(x,y_auto_2(x),'o-','DisplayName','exponential','color',[1,.5,0],'linewidth',2)
ax.XTick=x;
ax.YLabel.String='\gamma'
ax.YTick=[];
ax.XLabel.String='word position'
hold on
a=plot(x,cumsum(y,2),'o--','DisplayName','no kernel','color',[.5,.5,.5],'linewidth',2)
ax.XTick=x;
text(ax.XLim(1),ax.YLim(2)+.6,'D','HorizontalAlignment','right','FontSize',20,'FontWeight','bold')
%legend('show','Location','northwest')
ax.FontSize=16;
arr.MaxHeadSize=0;
ax.XLim=[1,9];
        if ~exist(strcat(analysis_path))
            mkdir(strcat(analysis_path))
        end 

% 

print(f, '-djpeg', strcat(analysis_path,'/norm_sum.jpeg'));

%% 
