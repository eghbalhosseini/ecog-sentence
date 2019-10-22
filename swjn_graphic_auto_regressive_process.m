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
