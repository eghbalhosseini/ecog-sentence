clear all
close all
data_path='/Users/eghbalhosseiniasl1/MyData/ecog-sentence';
subject_id='AMC';
d= dir([data_path,strcat('/**/',subject_id,'*_crunched.mat')]);
fprintf(' %d .mat files were found \n', length(d));
gamma_band_index=4;

%
subject_ids=table2cell(unique(cell2table(cellfun(@(x) x(regexpi(x,'AMC')+[0:5]), {d.name},'UniformOutput',false)')));
%

save_path='/Users/eghbalhosseiniasl1/MyData/ECoG-sentence/crunched/';
analysis_path='/Users/eghbalhosseiniasl1/MyData/ECoG-sentence/analysis/swjn_graphic_distributional_compositions';
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
set(f,'position',[-1743 834 970 970]);
ax=axes('position',[.1,.1,.35,.35]);
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
W=V-U;
arr=quiver3(U(1),U(2),U(3),W(1),W(2),W(3),0,'Color',[.2,.2,.2]);
text(U(1)+W(1)/2,U(2)+W(2)/2,U(3)+W(3)/2+.12,'V_{w_i}-V_{w_j}','HorizontalAlignment','center','fontsize',16)
arr.MaxHeadSize=.5;
arrayfun(@(x) set(x,'linewidth',3),arr)
arr=quiver3(0,0,0,U(1),U(2),U(3),0,'r');
arr.MaxHeadSize=.5;
text(U(1)/2-.05,U(2)/2,U(3)/2,'V_{w_i}','HorizontalAlignment','left','fontsize',16)
arrayfun(@(x) set(x,'linewidth',2),arr)
arr=quiver3(0,0,0,V(1),V(2),V(3),0,'b');
text(V(1)/2-.05,V(2)/2,V(3)/2,'V_{w_j}','HorizontalAlignment','left','fontsize',16)
arr.MaxHeadSize=.8;
arrayfun(@(x) set(x,'linewidth',2),arr)
a=annotation(f,'textbox',[.1 .41 .7 .05],'String',...
    'A',...
    'FontSize',15,'FontWeight','bold','LineStyle','none');
 %%%%%%%%%%%%%%%%%%%%%%%
[X,Y] = meshgrid(cosd(theta)',sind(theta)');
ax=axes('position',[.45,.1,.35,.35]);
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

a=annotation(f,'textbox',[.45 .41 .7 .05],'String',...
    'B',...
    'FontSize',15,'FontWeight','bold','LineStyle','none');

W=W-V;
arr=quiver3(V(1),V(2),V(3),W(1),W(2),W(3),0,'Color',[.2,.2,.2],'LineStyle','--');
arr.MaxHeadSize=0;
        if ~exist(strcat(analysis_path))
            mkdir(strcat(analysis_path))
        end 
print(f, '-djpeg', strcat(analysis_path,'/distributional_dist.jpeg'));

%% 
