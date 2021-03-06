clear all
close all
data_path='/Users/eghbalhosseiniasl1/MyData/ecog-sentence';

%
save_path='/Users/eghbalhosseiniasl1/MyData/ECoG-sentence/crunched/';
analysis_path='/Users/eghbalhosseiniasl1/MyData/ECoG-sentence/analysis/swjn_graphic_probe_analysis';
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
close all 
theta=[0:15:90];
x=cosd(theta)
y=sind(theta)
z=theta*0+0;
 
[X,Y] = meshgrid(cosd(theta)',sind(theta)');

f=figure;
aspect_ration=9.32./4.13;
y=600;
set(f,'position',[-798   598 aspect_ration*y y]);


ax=axes('position',[.5,.1,.5,.7]);
hold on
arr=arrayfun(@(x,y) quiver3(0*x,0*x,0*x,x,y,x*0,0,'color',[.4,.4,.4]), X(1,1:end),transpose(Y(1:end,1)))
arrayfun(@(x) set(x,'linewidth',1.5),arr)
arr=quiver3(0,0,0,0,0,1,0,'color',[.4,.4,.4]);
arrayfun(@(x) set(x,'linewidth',2),arr);
set(ax,'View',[110.2611 8.8477]);

% 
x=[.8,.85,.5,.2,.23,.62,.55,.52];
y=[.1,.21,.15,.08,.19,.5,.58,.66];
z=[.42,.41,.5,.56,.56,.44,.47,.47];
npts = 8;
sent_text=strsplit("The dog chased the cat all day long",' ')
%t = linspace(0,2*pi,npts);
%z = (sin(t/3)+.1*t)/1.5+.05;
%x=0.1*t;
%y=1-0.1613*t;
xyz = [x; y; z];
h0=plot3(xyz(1,:),xyz(2,:),xyz(3,:),'o','LineWidth',2,'Color',[.2,.2,.2]);
h0.MarkerFaceColor=[0,0,1];
h0.MarkerSize=7
text(xyz(1,:)+.02,xyz(2,:)+.02,xyz(3,:)+.02,sent_text,...
    'VerticalAlignment','bottom','fontsize',14,'HorizontalAlignment','left','FontWeight','bold')
hold on
fnplt(cscvn(xyz),'k-',2);
child=get(ax,'Children');
child(1).Color=[.4,.4,.4];
uistack(child(1),'bottom');
%
ax.YAxis.Visible='off';
ax.XAxis.Visible='off';
ax.ZAxis.Visible='off';


arrayfun(@(x,y,z) quiver3(0,0,0,x,y,z,0,'Color',[0,.5,1],'LineWidth',2,'MaxHeadSize',.3),x,y,z)


% 
norm_vect=[mean(x),mean(y),mean(z)];
n_x=norm_vect(1);
n_y=norm_vect(2);
n_z=norm_vect(3);
ax = atan2d(sqrt(n_y^2+n_z^2),n_x);
ay = atan2d(sqrt(n_z^2+n_x^2),n_y);
az = atan2d(sqrt(n_x^2+n_y^2),n_z);

[x_p y_p] = meshgrid(0:0.1:1); % Generate x and y data

z_p = x_p*0; % Solve for z data
h=surf(x_p,y_p,z_p+n_z) %Plot the surface

direction = [1 0 0];
rotate(h,direction,-3)
direction = [0 1 0];
rotate(h,direction,13)
% 
h.FaceColor=[.6,.2,1];
h.FaceAlpha=.2;
h.EdgeColor='none'
%direction = [0 0 1];
%rotate(h,direction,az)

if ~exist(strcat(analysis_path))
            mkdir(strcat(analysis_path))
end 
print(f, '-djpeg', strcat(analysis_path,'/probe_presentation.jpeg'));
%% 
x=2:8;
y=rand(1,7) ;
npts=3;
c_1=[1,1,1];
c_2=exp(-.3*[1,2,3]);
f=figure;
aspect_ration=9.32./4.13;
size_plot=500;
set(f,'position',[591 455 aspect_ration*size_plot size_plot]);
ax=axes('position',[.05,.1,.3,.6]);
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
ax.FontSize=12;
ax.FontWeight='bold'

%text(ax.XLim(1),ax.YLim(2)+.2,'A','HorizontalAlignment','right','FontSize',20,'FontWeight','bold')
ax=axes('position',[.45,.4,.1,.2]);
a=stem(1:length(c_1),c_1,'Color',[1,0,.5],'LineWidth',2)
box off
ax.YTick=[];
ax.XTick=[];
%text(ax.XLim(1),ax.YLim(2)+.6,'B','HorizontalAlignment','right','FontSize',20,'FontWeight','bold')
ax=axes('position',[.45,.1,.1,.2]);
a=stem(1:length(c_2),c_2,'Color',[1,.5,0],'LineWidth',2)
box off
ax.XTickLabel={'c_1','c_2','c_3'};
ax.YTick=[];
%text(ax.XLim(1),ax.YLim(2)+.5,'C','HorizontalAlignment','right','FontSize',20,'FontWeight','bold')
ax.FontSize=12
ax.FontWeight='bold'

ax=axes('position',[.65,.1,.3,.6]);
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
%text(ax.XLim(1),ax.YLim(2)+.6,'D','HorizontalAlignment','right','FontSize',20,'FontWeight','bold')
%legend('show','Location','northwest')
ax.FontSize=12
ax.FontWeight='bold'
arr.MaxHeadSize=0;
ax.XLim=[1,9];
        if ~exist(strcat(analysis_path))
            mkdir(strcat(analysis_path))
        end 

% 

print(f, '-djpeg', strcat(analysis_path,'/norm_sum_presentation.jpeg'));

%% 
