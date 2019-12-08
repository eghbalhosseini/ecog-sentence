clear all
close all
data_path='/Users/eghbalhosseiniasl1/MyData/ecog-sentence';

%
save_path='/Users/eghbalhosseiniasl1/MyData/ECoG-sentence/crunched/';
analysis_path='/Users/eghbalhosseiniasl1/MyData/ECoG-sentence/analysis/swjn_Aim_3_graphics';
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
emb = fastTextWordEmbedding;
%%
random_set=randperm(5000);
all_vectors_w2v=cell2mat(arrayfun(@(x) word2vec(emb,emb.Vocabulary(x)),...
    random_set,'UniformOutput',false)');

%% 
vectors_1=all_vectors_w2v(randperm(5000,50),:);
vectors_2=all_vectors_w2v(randperm(5000,50),:);

% 
corr_result={};
k = 20;
ColorTag_1 = cbrewer('div', 'BrBG', k);
ColorTag_2 = cbrewer('div', 'RdYlGn', k);
close all 
f=figure;
aspect_ration=9.32./4.13;
y=600;
set(f,'position',[-1762 668 aspect_ration*y y]);
    
vectors=vectors_1;
[idx_clust_1,~,~] = spectralcluster(vectors,k);

vectors=vectors_2;
[idx_clust_2,~,~] = spectralcluster(vectors,k);
    
if ~exist(strcat(analysis_path))
            mkdir(strcat(analysis_path))
        end 

for p=1:100
    p
    idx_clust=idx_clust_1;
    half_clust=(unique(idx_clust));
    half_clust=half_clust(randperm(k,floor(k/2)));
    b=arrayfun(@(x) find(idx_clust==x), half_clust, 'UniformOutput', false);
    idx_sel=cell2mat(cellfun(@(x) x(randi(length(x))),b,'UniformOutput',false));
    idx_sel_1=idx_sel;
    % 
    idx_clust=idx_clust_2;
    half_clust=(unique(idx_clust));
    half_clust=half_clust(randperm(k,floor(k/2)));
    b=arrayfun(@(x) find(idx_clust==x), half_clust, 'UniformOutput', false);
    idx_sel=cell2mat(cellfun(@(x) x(randi(length(x))),b,'UniformOutput',false));
    idx_sel_2=idx_sel;
%   
    vector_1_spec=[vectors_1(idx_sel_1,:);vectors_1(idx_sel_2,:)];
    vector_2_spec=[vectors_2(idx_sel_1,:);vectors_2(idx_sel_2,:)];
    similarity_spec_1=squareform(pdist(vector_1_spec,'correlation'));
    similarity_spec_2=squareform(pdist(vector_2_spec,'correlation'));
%
    Y_1=flatten_RDM(similarity_spec_1,1);
    Y_2=flatten_RDM(similarity_spec_2,1);
    % 
    A=corrcoef(Y_1,Y_2);
    corr_result=[corr_result;{[vectors_1(idx_sel_1,:);vectors_1(idx_sel_2,:)]},...
        {[vectors_2(idx_sel_1,:);vectors_2(idx_sel_2,:)]},{A(1,2)}];

    
    if p==1
        ax=axes('position',[0.01,0.3,.15,.4]);
        x = 1:2;
        y = 1:50;
        [X,Y] = meshgrid(x,y);
        C=reshape(repmat([1,.5,.5],50,1,2),[],2,3);
        C=ones([size(X,1),size(X,2),3]);
        a=surf(X,Y,5+0*X,C);
        a.FaceAlpha=1;
        hold on 
        a.FaceAlpha=1;
        
        set(ax,'view',[0,90])
        set(ax,'ylim',[1,50])
        daspect([1,1.5,1.5]);
        grid off
        ax.YAxis.Visible='off';
        ax.XAxis.Visible='off';
        % 
        
        ax=axes('position',[0.05,0.05,.15,.4]);
        x = 1:2;
        y = 1:50;
        [X,Y] = meshgrid(x,y);
        C=reshape(repmat([1,.5,.5],50,1,2),[],2,3);
        %C=ones([size(X,1),size(X,2),3]);
        a=surf(X,Y,5+0*X,C);
        a.FaceAlpha=1;
        hold on 
        a.FaceAlpha=1;
        
        set(ax,'view',[0,90])
        set(ax,'ylim',[1,50])
        daspect([1,1.5,1.5]);
        grid off
        ax.YAxis.Visible='off';
        ax.XAxis.Visible='off';
        % 
        ax=axes('position',[0.05,0.5,.15,.4]);
        C=reshape(repmat([204,229,255]/255,50,1,2),[],2,3);
        a=surf(X,Y,5+0*X,C)
        a.FaceAlpha=1;
        
        set(ax,'view',[0,90])
        set(ax,'ylim',[0,50])
        daspect([1,1.5,1.5]);
        grid off
        ax.YAxis.Visible='off';
        ax.XAxis.Visible='off';
        % 
        similarity_1=squareform(pdist(vectors_1,'correlation'));
        similarity_2=squareform(pdist(vectors_2,'correlation'));
        x = 1:length(vectors_1);
        
        % 
        ax=axes('position',[0.15,0.05,.175,.175*aspect_ration]);
        a=imagesc(x,x,similarity_1)
        colors = cbrewer('div', 'RdBu', 128);
        colors = flipud(colors); % puts red on top, blue at the bottom
        colormap(colors);
        
        daspect([1,1,1]);
        set(gca, 'ydir', 'reverse','box','on');
        ax.XTick=[];
        ax.YTick=[];
        
        % 
        ax=axes('position',[0.15,0.5,.175,.175*aspect_ration]);
        a=imagesc(x,x,similarity_2)
        colors = cbrewer('div', 'RdBu', 128);
        colors = flipud(colors); % puts red on top, blue at the bottom
        colormap(colors);
        
        daspect([1,1,1]);
        set(gca, 'ydir', 'reverse','box','on');
        ax.XTick=[];
        ax.YTick=[];
        % 
        vectors=vectors_1;
        [idx_clust_1,~,~] = spectralcluster(vectors,k);
        [~, SCORE_1] = pca(vectors);
        half_clust=(unique(idx_clust_1));
        half_clust=half_clust(randperm(k,floor(k/2)));
        b=arrayfun(@(x) find(idx_clust_1==x), half_clust, 'UniformOutput', false);
        idx_sel=cell2mat(cellfun(@(x) x(randi(length(x))),b,'UniformOutput',false));
        idx_sel1=idx_sel;
        %
        vectors=vectors_2;
        [idx_clust_2,~,~] = spectralcluster(vectors,k);
        [~, SCORE_2] = pca(vectors);
        half_clust=(unique(idx_clust_2));
        half_clust=half_clust(randperm(k,floor(k/2)));
        b=arrayfun(@(x) find(idx_clust_2==x), half_clust, 'UniformOutput', false);
        idx_sel=cell2mat(cellfun(@(x) x(randi(length(x))),b,'UniformOutput',false));
        idx_sel2=idx_sel;
        
        %%%%%%%%%%%%%%spectral clustering  
        ax=axes('position',[0.37,0.05,.175,.175*aspect_ration]);
        h1=gscatterEH([SCORE_1(:,1),SCORE_1(:,2)],idx_clust_1,ColorTag_1);
        arrayfun(@(x) set(x,'MarkerFaceAlpha',.4), h1.Children)
        arrayfun(@(x) set(x,'SizeData',100), h1.Children)
        arrayfun(@(x) set(x,'MarkerEdgeColor','none'), h1.Children)
        hold on 
        h2=gscatterEH([SCORE_1(idx_sel1,1),SCORE_1(idx_sel1,2)],idx_clust_1(idx_sel1),ColorTag_1);
        arrayfun(@(x) set(x,'linewidth',1), h2.Children);
        arrayfun(@(x) set(x,'SizeData',120), h2.Children)
        ax.XTick=[];
        ax.YTick=[];
       
        % 
        ax=axes('position',[0.37,0.5,.175,.175*aspect_ration]);
        h1=gscatterEH([SCORE_2(:,1),SCORE_2(:,2)],idx_clust_2,ColorTag_2);
        arrayfun(@(x) set(x,'MarkerFaceAlpha',.3), h1.Children)
        arrayfun(@(x) set(x,'SizeData',100), h1.Children)
        arrayfun(@(x) set(x,'MarkerEdgeColor','none'), h1.Children)
        
        hold on 
        h2=gscatterEH([SCORE_2(idx_sel2,1),SCORE_2(idx_sel2,2)],idx_clust_2(idx_sel2),ColorTag_2);
        arrayfun(@(x) set(x,'linewidth',1), h2.Children);
        arrayfun(@(x) set(x,'SizeData',120), h2.Children)
        ax.XTick=[];
        ax.YTick=[];
       % 
       % 
               ax=axes('position',[0.51,0.3,.15,.4]);
        x = 1:2;
        y = 1:20;
        [X,Y] = meshgrid(x,y);
        C=reshape(repmat([1,.5,.5],20,1,2),[],2,3);
        C=ones([size(X,1),size(X,2),3]);
        a=surf(X,Y,5+0*X,C);
        a.FaceAlpha=1;
        hold on 
        a.FaceAlpha=1;
        
        set(ax,'view',[0,90])
        set(ax,'ylim',[1,20])
        daspect([1,1.5,1.5]);
        grid off
        ax.YAxis.Visible='off';
        ax.XAxis.Visible='off';
        % 
       ax=axes('position',[0.55,0.05,.15,.4]);
        x = 1:2;
        y = 1:20;
        [X,Y] = meshgrid(x,y);
        C=reshape(repmat([1,.5,.5],20,1,2),[],2,3);
        %C=ones([size(X,1),size(X,2),3]);
        a=surf(X,Y,5+0*X,C);
        a.FaceAlpha=1;
        hold on 
        a.FaceAlpha=1;
        
        set(ax,'view',[0,90])
        set(ax,'ylim',[1,20])
        daspect([1,1.5,1.5]);
        grid off
        ax.YAxis.Visible='off';
        ax.XAxis.Visible='off';
        % 
        ax=axes('position',[0.55,0.48,.15,.41]);
        C=reshape(repmat([204,229,255]/255,20,1,2),[],2,3);
        a=surf(X,Y,5+0*X,C)
        a.FaceAlpha=1;
        
        set(ax,'view',[0,90])
        set(ax,'ylim',[0,20])
        daspect([1,1.5,1.5]);
        grid off
        ax.YAxis.Visible='off';
        ax.XAxis.Visible='off';
        % 
        
        
        vector_1_spec=[vectors_1(idx_sel1,:);vectors_1(idx_sel2,:)];
        vector_2_spec=[vectors_2(idx_sel1,:);vectors_2(idx_sel2,:)];
        similarity_spec_1=squareform(pdist(vector_1_spec,'correlation'));
        similarity_spec_2=squareform(pdist(vector_2_spec,'correlation'));
        x = 1:size(vector_1_spec,1);
        
        ax=axes('position',[0.66,0.05,.175,.175*aspect_ration]);
        a=imagesc(x,x,similarity_spec_2)
        colors = cbrewer('div', 'RdBu', 128);
        colors = flipud(colors); % puts red on top, blue at the bottom
        colormap(colors);
        
        daspect([1,1,1]);
        set(gca, 'ydir', 'reverse','box','on');
        ax.XTick=[];
        ax.YTick=[];
        
        
        ax=axes('position',[0.66,0.5,.175,.175*aspect_ration]);
        a=imagesc(x,x,similarity_spec_1)
        colors = cbrewer('div', 'RdBu', 128);
        colors = flipud(colors); % puts red on top, blue at the bottom
        colormap(colors);
        set(gca, 'ydir', 'reverse','box','on');
       ax.XTick=[];
       ax.YTick=[];
        
        daspect([1,1,1]);
        
    end 
    
end 
%
corr_=sort([corr_result{:,3}]);
ax=axes('position',[0.9,0.1,.07,.8]);
a=scatter(corr_,1:length(corr_));
xLims=ax.XLim;
ax.XLim=[-1,1]*1.05*max(xLims);
ax.YLim=[-5,105];
a.SizeData=25;
a.MarkerEdgeColor=[0,0,0];
a.LineWidth=1.5
a.MarkerFaceColor=[.9,.9,.9];
set(gca, 'ydir', 'normal','box','off');
hold on 
ax.YAxis.Visible='off';
ax.XTick=[ax.XLim(1),0,ax.XLim(2)];
ax.XTickLabel={'-','0','+'};
ax.FontSize=12;
ax.Title.String='Corr';
print(f, '-djpeg', strcat(analysis_path,'/sentence_selection_presentation_2.jpeg'));

threshold=0.03;
a=plot(ax.XLim*0+threshold,ax.YLim,'LineStyle','--','color','k');
a=plot(ax.XLim*0-threshold,ax.YLim,'LineStyle','--','color','k');
passed=(corr_<threshold & corr_>-threshold);
print(f, '-djpeg', strcat(analysis_path,'/sentence_selection_presentation_3.jpeg'));
a=scatter(corr_(passed),find(passed),'filled');
a.SizeData=25;
a.MarkerEdgeColor=[0,0,0];
a.MarkerFaceColor=[1,0,0];
a.LineWidth=1.5

ax.YAxis.Visible='off';
ax.XTick=[ax.XLim(1),0,ax.XLim(2)];
ax.XTickLabel={'-','0','+'};
ax.FontSize=12;
ax.Title.String='Corr';
%text(ax.XLim(1),ax.YLim(2)+10,'F','HorizontalAlignment','right','FontSize',20,'FontWeight','bold')

        if ~exist(strcat(analysis_path))
            mkdir(strcat(analysis_path))
        end 

print(f, '-djpeg', strcat(analysis_path,'/sentence_selection_presentation_4.jpeg'));

%%
%% schematic of 
vectors_1=all_vectors_w2v(randperm(100,20),:);
vectors_2=all_vectors_w2v(randperm(100,20),:);
vectors_3=all_vectors_w2v(randperm(100,20),:);
similarity_1=squareform(pdist(vectors_1,'correlation'));
similarity_2=squareform(pdist(vectors_2,'correlation'));
similarity_3=squareform(pdist(vectors_3,'correlation'));
%
close all 
f=figure;
aspect_ration=9.32./4.13;
y=500;
set(f,'position',[591 455 aspect_ration*y y]);
%
ax=axes('position',[0.01,0.55,.15,.4]);
x = 1:2;
y = 1:50;
[X,Y] = meshgrid(x,y);
C=ones([size(X,1),size(X,2),3]);
a=surf(X,Y,5+0*X,C);
a.FaceAlpha=1;
hold on
set(ax,'view',[0,90])
set(ax,'ylim',[0,55])
daspect([1,1.5,1]);
grid off
ax.YAxis.Visible='off';
ax.XAxis.Visible='off';

ax=axes('position',[0.26,0.55+.01*aspect_ration,.15,.15*aspect_ration]);
a=imagesc(x,x,similarity_1)
colors = cbrewer('div', 'RdBu', 128);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);
set(gca, 'ydir', 'reverse','box','off');
daspect([1,1,1]);
set(gca, 'ydir', 'reverse','box','on');
ax.XTick=[];
ax.YTick=[];


ax=axes('position',[0.25,0.55,.15,.15*aspect_ration]);
a=imagesc(x,x,similarity_2)
colors = cbrewer('div', 'RdBu', 128);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);
set(gca, 'ydir', 'reverse','box','on');
ax.XTick=[];
ax.YTick=[];
%
ax=axes('position',[0.03,0.05+.01*aspect_ration,.15,.15*aspect_ration]);

AMC026_elec_fsaverage=load(['/Users/eghbalhosseiniasl1/MyData/ecog-sentence/subjects_raw/AMC026/IMAGING/','AMC026_elec_fsavg.mat']);
AMC026_elec_fsaverage=AMC026_elec_fsaverage.elec_fsavg_frs;
fspial_lh = ft_read_headshape('/Applications/freesurfer/subjects/fsaverage/surf/lh.pial');
fspial_lh.coordsys = 'fsaverage';
h=ft_plot_mesh(fspial_lh);
view([-100 10]);
material dull;
lighting gouraud;
camlight;
%
ax=axes('position',[0.25,0.1,.15,.15*aspect_ration]);
a=imagesc(x,x,similarity_3)
colors = cbrewer('div', 'RdBu', 128);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);
set(gca, 'ydir', 'reverse','box','on');
ax.XTick=[];
ax.YTick=[];
% 

Y_1=flatten_RDM(similarity_1,1);
Y_2=flatten_RDM(similarity_2,1);
Y_3=flatten_RDM(similarity_3,1);

A=corrcoef(Y_1,Y_3);
B=corrcoef(Y_2,Y_3);

ax=axes('position',[0.50,0.35,.20,.15*aspect_ration]);
bl=bar(1,[A(1,2)],'DisplayName','F_1');
bl.FaceColor=[1,.3,.2]
hold on 
bl.LineWidth=2
bl1=bar(2,[B(1,2)],'DisplayName','F_2');
bl1.LineWidth=2
bl1.FaceColor=[.4,1,.6]
ax.XTick=[];
ax.YTick=[];
set(ax,'box','off')
ax.YLabel.String='Corr'
ax.XAxis.LineWidth=2;
ax.YAxis.LineWidth=2;
ax.YLim(2)=1.1*ax.YLim(2);
ax.FontSize=12
ax.Title.String={'correlation between',' brain responses ', 'and features'}

print(f, '-djpeg', strcat(analysis_path,'/rsa_schematic.jpeg'));
%% 



    
