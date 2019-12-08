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

data_path='/Users/eghbalhosseiniasl1/MyData/ecog-sentence';

d_pmi= dir([data_path,'/**/*pmi_sentences.csv']);
fprintf(' %d pmi files were found \n', length(d_pmi));
i=1;
[pmi_cell]=generate_pmi_table_eh(strcat(d_pmi(i).folder,'/',d_pmi(i).name));
%% 
a=triu(ones(8),1)
a(a==0)=nan;
b=cellfun(@(x) x.*a, pmi_cell(:,2),'UniformOutput',false);
c=[b{:}];
d=reshape(c,[],1);
d(isnan(d))=[];
%% 
f=figure;
h=histfit(d',50,'kernel')
ylabel('frequency')
xlabel('PMI')
box off
x_fit=h(2).XData;
y_fit=h(2).YData;
print(f, '-djpeg', strcat(analysis_path,'/PMI_dist.jpeg'));
%% 
close all 
x=x_fit;
y=1:1000;
[X,Y]=meshgrid(x,y);

Z=normpdf(X,0,20)+.001*rand(size(X,1),size(X,2));
Z=Z./repmat(max(Z,[],2),1,size(Z,2));
f=figure;
set(f,'position',[827 604 1679 741]);
ax=axes('position',[0.05,0.05,.8,.8]);
a=surf(X,Y,Z);
a.EdgeColor='none'
daspect([10,50,.5])
set(gca,'view',[46.7185   26.1282]);

% make a trajectory

N = 12 ; % number of steps
dims=1;
% positions, starting at (0,0,...,0)
x_ind = transpose(cumsum(full(sparse(1:N, randi(dims,1,N), [0 2*randi([0 1],1,N-1)-1], N, dims)))) ; 
y_ind=ceil(linspace(10,12*N,N));

z_traj=arrayfun(@(x,y) Z(y,x),x_ind+50, y_ind)+0.05;
hold on 
traj=plot3(x_ind,y_ind,z_traj,'color','k','LineWidth',1,'Marker','.'); 
traj.MarkerSize=10;
grid off
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
ax.XAxisLocation = 'origin';
ax.YAxis.Visible='off';
ax.XAxis.Visible='off';
ax.ZAxis.Visible='off';

arr=quiver3(0,0,0,0,80,0,0,'k');
arr.MaxHeadSize=.5;
arr=quiver3(0,0,0,0,0,1,0,'k');

arr=quiver3(0,0,0,30,0,0,0,'k');
arr=quiver3(0,0,0,-30,0,0,0,'k');


if ~exist(strcat(analysis_path))
    mkdir(strcat(analysis_path))
end
print(f, '-djpeg', strcat(analysis_path,'/PMI.jpeg'));

%% 


for p=1:100
    
    
    if p==1
        ax=axes('position',[0.05,0.15,.15,.35]);
        x = 1:2;
        y = 1:50;
        [X,Y] = meshgrid(x,y);
        C=ones([size(X,1),size(X,2),3]);
        a=surf(X,Y,5+0*X,C);
        a.FaceAlpha=1;
        hold on 
       
        
        grid off
        ax.YAxis.Visible='off';
        ax.XAxis.Visible='off';
        text(ax.XLim(1),ax.YLim(2)+.15,'A','HorizontalAlignment','right','FontSize',12,'FontWeight','bold')
        ar = annotation('arrow');
        ar.Position=[.13,.16,.15,0.045];
        ar = annotation('arrow');
        ar.Position=[.14,.48,.15,-0.024];
        
        % 
        ax=axes('position',[0.15,0.2,.3,.3]);
        similarity_1=squareform(pdist(vectors_1,'correlation'));
        similarity_2=squareform(pdist(vectors_2,'correlation'));
        x = 1:length(vectors_1);
        [X,Y] = meshgrid(x,x);
        plot3([0,1000,1000,0,0],[0,0,1000,1000,0],[2,2,2,2,2],'Color','k');
        hold on ;
        plot3([0,1000,1000,0,0]+50,[0,0,1000,1000,0]-50,-1+[0,0,0,0,0],'Color','k');
        hold on ;
        a=surf(X,Y,similarity_1);
        colors = cbrewer('div', 'RdBu', 128);
        colors = flipud(colors); % puts red on top, blue at the bottom
        colormap(colors);
        
        a.EdgeColor='none';
        set(gca, 'ydir', 'reverse','box','off');
        set(ax,'view',[0,90]);
        hold on ;
        a=surf(X+50,Y-50,similarity_1);
        a.EdgeColor='none';
        daspect([1,1,1]);
        ax.YAxis.Visible='off';
        ax.XAxis.Visible='off';
        text(ax.XLim(1),ax.YLim(1)+.15,'B','HorizontalAlignment','right','FontSize',12,'FontWeight','bold')
        ar = annotation('arrow');
        ar.Position=[.40,.26,.046,-0.030];
        ar = annotation('arrow');
        ar.Position=[.413,.4,.033,+0.024];
        
        % 
        vectors=vectors_1;
        [idx_clust,~,~] = spectralcluster(vectors,k);
        [~, SCORE] = pca(vectors);
        half_clust=(unique(idx_clust));
        half_clust=half_clust(randperm(k,floor(k/2)));
        b=arrayfun(@(x) find(idx_clust==x), half_clust, 'UniformOutput', false);
        idx_sel=cell2mat(cellfun(@(x) x(randi(length(x))),b,'UniformOutput',false));
        similarity=squareform(pdist(vectors(idx_sel,:),'correlation'));
        % 
        ax=axes('position',[0.45,0.15,.15,.15]);
        h1=gscatterEH([SCORE(:,1),SCORE(:,2)],idx_clust,ColorTag);
        arrayfun(@(x) set(x,'MarkerFaceAlpha',.3), h1.Children)
        arrayfun(@(x) set(x,'MarkerEdgeColor','none'), h1.Children)
        hold on 
        h2=gscatterEH([SCORE(idx_sel,1),SCORE(idx_sel,2)],idx_clust(idx_sel),ColorTag);
        arrayfun(@(x) set(x,'linewidth',1), h2.Children);
        ax.XTick=[];
        ax.YTick=[];
        text(ax.XLim(1),ax.YLim(2)+.15,'D','HorizontalAlignment','right','FontSize',12,'FontWeight','bold')
        ar = annotation('arrow');
        ar.Position=[.61,.225,.03,0];
        

        
        % 
        ax=axes('position',[0.65,0.15,.15,.15]);
        x = 1:length(similarity);
        [X,Y] = meshgrid(x,x);
        plot3([1,max(x),max(x),1,1],[1,1,max(x),max(x),1],[2,2,2,2,2],'Color','k');
        hold on 
        a=surf(X,Y,similarity);
        colors = cbrewer('div', 'RdBu', 128);
        colors = flipud(colors); % puts red on top, blue at the bottom
        colormap(colors);
        a.EdgeColor='none';
        set(gca, 'ydir', 'reverse','box','off');
        set(ax,'view',[0,90]);
        axis tight
        ax.YAxis.Visible='off';
        ax.XAxis.Visible='off';
        daspect([1,1,1])
        text(ax.XLim(1)-2,ax.YLim(1)-2,'F','HorizontalAlignment','right','FontSize',12,'FontWeight','bold')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        vectors=vectors_2;
        [idx_clust,~,~] = spectralcluster(vectors,k);
        [~, SCORE] = pca(vectors);
        half_clust=(unique(idx_clust));
        half_clust=half_clust(randperm(k,floor(k/2)));
        b=arrayfun(@(x) find(idx_clust==x), half_clust, 'UniformOutput', false);
        idx_sel=cell2mat(cellfun(@(x) x(randi(length(x))),b,'UniformOutput',false));
        similarity=squareform(pdist(vectors(idx_sel,:),'correlation'));
        % 
        ax=axes('position',[0.45,0.35,.15,.15]);
        h1=gscatterEH([SCORE(:,1),SCORE(:,2)],idx_clust,ColorTag);
        arrayfun(@(x) set(x,'MarkerFaceAlpha',.3), h1.Children)
        arrayfun(@(x) set(x,'MarkerEdgeColor','none'), h1.Children)
        hold on 
        h2=gscatterEH([SCORE(idx_sel,1),SCORE(idx_sel,2)],idx_clust(idx_sel),ColorTag);
        arrayfun(@(x) set(x,'linewidth',1), h2.Children);
        ax.XTick=[];
        ax.YTick=[];
        text(ax.XLim(1),ax.YLim(2)+.15,'C','HorizontalAlignment','right','FontSize',12,'FontWeight','bold')
        ar = annotation('arrow');
        ar.Position=[.61,.425,.03,0];
        
        % 
        ax=axes('position',[0.65,0.35,.15,.15]);
        x = 1:length(similarity);
        [X,Y] = meshgrid(x,x);
        plot3([1,max(x),max(x),1,1],[1,1,max(x),max(x),1],[2,2,2,2,2],'Color','k');
        hold on 
        a=surf(X,Y,similarity);
        colors = cbrewer('div', 'RdBu', 128);
        colors = flipud(colors); % puts red on top, blue at the bottom
        colormap(colors);
        a.EdgeColor='none';
        set(gca, 'ydir', 'reverse','box','off');
        set(ax,'view',[0,90]);
        axis tight
        ax.YAxis.Visible='off';
        ax.XAxis.Visible='off';
        text(ax.XLim(1)-2,ax.YLim(1)-2,'E','HorizontalAlignment','right','FontSize',12,'FontWeight','bold')
        
    end 
    
end 

corr_=sort([corr_result{:,3}]);
ax=axes('position',[0.85,0.15,.1,.35]);


a=scatter(corr_,1:length(corr_));
xLims=ax.XLim;
ax.XLim=[-1,1]*1.05*max(xLims);
ax.YLim=[-5,105];
a.SizeData=15;
a.MarkerEdgeColor=[.3,.3,.3];
set(gca, 'ydir', 'reverse','box','off');

hold on 
a=plot(ax.XLim*0+0.02,ax.YLim,'LineStyle','--','color','k');
a=plot(ax.XLim*0-0.02,ax.YLim,'LineStyle','--','color','k');
passed=(corr_<0.01 & corr_>-0.01);
a=scatter(corr_(passed),find(passed),'filled');
a.SizeData=25;
a.MarkerEdgeColor=[0,0,0];
a.MarkerFaceColor=[1,0,0];

ax.YAxis.Visible='off';
ax.XTick=[ax.XLim(1),0,ax.XLim(2)];
ax.XTickLabel={'-','0','+'};
ax.FontSize=14;
ax.XLabel.String='Correlation';
text(ax.XLim(1),ax.YLim(1)-2,'G','HorizontalAlignment','right','FontSize',12,'FontWeight','bold')

        if ~exist(strcat(analysis_path))
            mkdir(strcat(analysis_path))
        end 

print(f, '-djpeg', strcat(analysis_path,'/sentence_selection.jpeg'));
%% 
ar = annotation('arrow');
ar.Position=[.61,.225,.03,0];
c = ar.Color;


%% 
figure;
k = 100;
[idx1_clust2,V,D] = spectralcluster(vectors_2,k);
[COEFF, SCORE] = pca(vectors_2);

half_clust_2=(unique(idx1_clust2));
half_clust_2=half_clust_2(randperm(k,floor(k/2)));
b=arrayfun(@(x) find(idx1_clust2==x), half_clust_2, 'UniformOutput', false);
idx_sel_2=cell2mat(cellfun(@(x) x(randi(length(x))),b,'UniformOutput',false));
ColorTag = cbrewer('qual', 'Paired', k);

h1=gscatterEH([SCORE(:,1),SCORE(:,2)],idx1_clust2,ColorTag);
arrayfun(@(x) set(x,'MarkerFaceAlpha',.3), h1.Children)
arrayfun(@(x) set(x,'MarkerEdgeColor','none'), h1.Children)
hold on 

h2=gscatterEH([SCORE(idx_sel_2,1),SCORE(idx_sel_2,2)],idx1_clust1(idx_sel_2),ColorTag);
arrayfun(@(x) set(x,'linewidth',2), h2.Children)

%% 

figure
colors = cbrewer('div', 'RdBu', 128);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);
imagesc(similarity_1,[0,2]);
daspect([1,1,1])
set(gca, 'ydir', 'reverse','box','off');



