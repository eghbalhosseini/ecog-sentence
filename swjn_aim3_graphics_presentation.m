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
random_set=randperm(10000);
all_vectors_w2v=cell2mat(arrayfun(@(x) word2vec(emb,emb.Vocabulary(x)),...
    random_set,'UniformOutput',false)');

%% 
vectors_1=all_vectors_w2v(randperm(10000,1000),:);
vectors_2=all_vectors_w2v(randperm(10000,1000),:);

%% 
corr_result={};
k = 100;
ColorTag = cbrewer('qual', 'Paired', k);
close all 
f=figure;
aspect_ration=9.32./4.13;
y=600;
set(f,'position',[591 455 aspect_ration*y y]);
    
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
        ax=axes('position',[0.01,0.1,.15,.8]);
        x = 1:2;
        y = 1:50;
        [X,Y] = meshgrid(x,y);
        C=ones([size(X,1),size(X,2),3]);
        a=surf(X,Y,5+0*X,C);
        a.FaceAlpha=1;
        hold on 
        a=surf(X+.75,Y+2.5,4+0*X,C)
        a.FaceAlpha=1;
        
        set(ax,'view',[0,90])
        set(ax,'ylim',[0,55])
        daspect([1,1.5,1]);
        grid off
        ax.YAxis.Visible='off';
        ax.XAxis.Visible='off';
        %text(ax.XLim(1)-.5,ax.YLim(2)-4,'A','HorizontalAlignment','right','FontSize',20,'FontWeight','bold')
        %ar = annotation('arrow');
        %ar.Position=[.12,.25,.1,0.045];
        %ar = annotation('arrow');
        %ar.Position=[.125,.48,.1,-0.024];
        
        similarity_1=squareform(pdist(vectors_1,'correlation'));
        similarity_2=squareform(pdist(vectors_2,'correlation'));
        x = 1:length(vectors_1);
        [X,Y] = meshgrid(x,x);
        % 
        ax=axes('position',[0.12+.01,0.15+.01*aspect_ration,.3,.3*aspect_ration]);
        a=imagesc(x,x,similarity_2)
        colors = cbrewer('div', 'RdBu', 128);
        colors = flipud(colors); % puts red on top, blue at the bottom
        colormap(colors);
        set(gca, 'ydir', 'reverse','box','off');
        daspect([1,1,1]);
        set(gca, 'ydir', 'reverse','box','on');
        ax.XTick=[];
        ax.YTick=[];
        
        
        
        % 
        ax=axes('position',[0.12,0.15,.3,.3*aspect_ration]);
        a=imagesc(x,x,similarity_1)
        colors = cbrewer('div', 'RdBu', 128);
        colors = flipud(colors); % puts red on top, blue at the bottom
        colormap(colors);
        set(gca, 'ydir', 'reverse','box','off');
        daspect([1,1,1]);
        set(gca, 'ydir', 'reverse','box','on');
        ax.XTick=[];
        ax.YTick=[];
        print(f, '-djpeg', strcat(analysis_path,'/sentence_selection_presentation_0.jpeg'));
        
        
        ar = annotation('arrow');
        ar.Position=[.42,.3,.045,0];
        ar.LineWidth=2
        ar = annotation('arrow');
        ar.Position=[.43,.7,.036,0];
        ar.LineWidth=2
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
        ax=axes('position',[0.47,0.13,.15,.15*aspect_ration]);
        h1=gscatterEH([SCORE_1(:,1),SCORE_1(:,2)],idx_clust_1,ColorTag);
        arrayfun(@(x) set(x,'MarkerFaceAlpha',.3), h1.Children)
        arrayfun(@(x) set(x,'MarkerEdgeColor','none'), h1.Children)
        hold on 
        h2=gscatterEH([SCORE_1(idx_sel1,1),SCORE_1(idx_sel1,2)],idx_clust_1(idx_sel1),ColorTag);
        arrayfun(@(x) set(x,'linewidth',1), h2.Children);
        ax.XTick=[];
        ax.YTick=[];
        %text(ax.XLim(1),ax.YLim(2)+.15,'C','HorizontalAlignment','right','FontSize',20,'FontWeight','bold')
        % 
        ax=axes('position',[0.47,0.53,.15,.15*aspect_ration]);
        h1=gscatterEH([SCORE_2(:,1),SCORE_2(:,2)],idx_clust_2,ColorTag);
        arrayfun(@(x) set(x,'MarkerFaceAlpha',.3), h1.Children)
        arrayfun(@(x) set(x,'MarkerEdgeColor','none'), h1.Children)
        hold on 
        h2=gscatterEH([SCORE_2(idx_sel2,1),SCORE_2(idx_sel2,2)],idx_clust_2(idx_sel2),ColorTag);
        arrayfun(@(x) set(x,'linewidth',1), h2.Children);
        ax.XTick=[];
        ax.YTick=[];
        %text(ax.XLim(1),ax.YLim(2)+.15,'B','HorizontalAlignment','right','FontSize',20,'FontWeight','bold')
        print(f, '-djpeg', strcat(analysis_path,'/sentence_selection_presentation_1.jpeg'));
        ar = annotation('arrow');
        ar.Position=[.6,.3,.042,.1];
        ar.LineWidth=2
        ar = annotation('arrow');
        ar.Position=[.6,.7,.042,-.1];
        ar.LineWidth=2
        %%%%%%%%%%%%%%%%% correlation 
        
        
        
        vector_1_spec=[vectors_1(idx_sel1,:);vectors_1(idx_sel2,:)];
        vector_2_spec=[vectors_2(idx_sel1,:);vectors_2(idx_sel2,:)];
        similarity_spec_1=squareform(pdist(vector_1_spec,'correlation'));
        similarity_spec_2=squareform(pdist(vector_2_spec,'correlation'));
        x = 1:size(vector_1_spec,1);
        
        ax=axes('position',[0.66,0.25+.01*aspect_ration,.2,.2*aspect_ration]);
        a=imagesc(x,x,similarity_spec_2)
        colors = cbrewer('div', 'RdBu', 128);
        colors = flipud(colors); % puts red on top, blue at the bottom
        colormap(colors);
        set(gca, 'ydir', 'reverse','box','off');
        daspect([1,1,1]);
        set(gca, 'ydir', 'reverse','box','on');
        ax.XTick=[];
        ax.YTick=[];
        
        
        ax=axes('position',[0.65,0.25,.2,.2*aspect_ration]);
        a=imagesc(x,x,similarity_spec_1)
        colors = cbrewer('div', 'RdBu', 128);
        colors = flipud(colors); % puts red on top, blue at the bottom
        colormap(colors);
        set(gca, 'ydir', 'reverse','box','on');
       ax.XTick=[];
       ax.YTick=[];
        
        daspect([1,1,1]);
        %ax.YAxis.Visible='off';
        %ax.XAxis.Visible='off';
        %text(ax.XLim(1),ax.YLim(1)+.15,'B','HorizontalAlignment','right','FontSize',20,'FontWeight','bold')
        %ar = annotation('arrow');
        %ar.Position=[.40,.25,.046,0];
        %ar = annotation('arrow');
        %ar.Position=[.413,.4,.033,0];

        % 
       %text(ax.XLim(1)-.5,ax.YLim(1)-5,'D','HorizontalAlignment','right','FontSize',20,'FontWeight','bold')
    end 
    
end 
%
corr_=sort([corr_result{:,3}]);
ax=axes('position',[0.9,0.15,.07,.7]);
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

threshold=0.01;
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
vectors_1=all_vectors_w2v(randperm(100,100),:);
vectors_2=all_vectors_w2v(randperm(100,100),:);
vectors_3=all_vectors_w2v(randperm(100,100),:);
similarity_1=squareform(pdist(vectors_1,'correlation'));
similarity_2=squareform(pdist(vectors_2,'correlation'));
similarity_3=squareform(pdist(vectors_3,'correlation'));
%% 
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
hold on 
bl.LineWidth=2
bl1=bar(2,[B(1,2)],'DisplayName','F_2');
bl1.LineWidth=2
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
close all 
f=figure;
aspect_ration=9.32./4.13;
y=600;
set(f,'position',[591 455 aspect_ration*y y]);
%
ax=axes('position',[0.05,0.5,.3,.45]);
voxels=150;
sentence=50;
voxel_pattern=rand(sentence,voxels);
a=imagesc(sentence,voxels,voxel_pattern);
colors = cbrewer('div', 'RdBu', 128);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);
set(gca, 'ydir', 'reverse','box','off');
daspect([1,1,1]);
set(gca, 'ydir', 'reverse','box','on');
ax.XTick=[];
ax.YLabel.String='Sentences';
ax.XLabel.FontWeight='bold'
ax.XLabel.String='Voxel pattern';
ax.YLabel.FontWeight='bold'
ax.YTick=[];
% 
a=annotation(f,'textbox',[.36 .7 .1 .1],'String','\approx','FontSize',25);
a.LineStyle='none';
% 
ax=axes('position',[0.45,0.4,.15,.5]);
a=text(0,0,'Response Pattern 1','FontSize',12,'FontWeight','bold');
a.Rotation=90;
a.Color=[1,0,0];
a.VerticalAlignment='top'
% 
a=text(.2,0,'Response Pattern 2','FontSize',12,'FontWeight','bold');
a.Rotation=90;
a.Color=[0,0,1];
a.VerticalAlignment='top'
a=text(.4,.5,'...','FontSize',12,'FontWeight','bold');
% 
a=text(.6,0,'Response Pattern N','FontSize',12,'FontWeight','bold');
a.Rotation=90;
a.Color=[0,.8,0];
a.VerticalAlignment='top'
ax.XAxis.Color=[1,1,1]
ax.XTick=[];
ax.YAxis.Color=[1,1,1]
ax.YTick=[];
%ax.XAxis.Visible='off'
%ax.YAxis.Visible='off'
ax.XLim=[-0.05,.9]
ax.YLim=[-0.05,.9]
ax.XLabel.String='Component';
ax.XLabel.FontWeight='bold'
ax.XLabel.Color=[0,0,0]
ax.YLabel.String='Sentences';
ax.YLabel.Color=[0,0,0]
ax.YLabel.FontWeight='bold'

a=annotation(f,'textbox',[.6 .7 .1 .1],'String','\times','FontSize',25);
a.LineStyle='none';
% 

ax=axes('position',[0.7,0.6,.25,.3]);
a=text(0,0,'Voxel Weight Pattern 1','FontSize',12,'FontWeight','bold');
a.Rotation=0;
a.Color=[1,0,0];
a.VerticalAlignment='top'
% 
a=text(0,.2,'Voxel Weight Pattern 2','FontSize',12,'FontWeight','bold');
a.Rotation=0;
a.Color=[0,0,1];
a.VerticalAlignment='top'
a=text(.4,.55,'...','FontSize',12,'FontWeight','bold');
a.Rotation=90;
a.VerticalAlignment='middle'
% 
a=text(.0,.6,'Voxel Weight Pattern N','FontSize',12,'FontWeight','bold');
a.Rotation=0;
a.Color=[0,.8,0];
a.VerticalAlignment='top'
ax.XAxis.Color=[1,1,1]
ax.XTick=[];
ax.YAxis.Color=[1,1,1]
ax.YTick=[];
%ax.XAxis.Visible='off'
%ax.YAxis.Visible='off'
ax.XLim=[-0.05,.9]
ax.YLim=[-0.05,.9]
ax.YLabel.String='Component';
ax.XLabel.FontWeight='bold'
ax.XLabel.Color=[0,0,0]
ax.XLabel.String='Voxel Weight';
ax.YLabel.Color=[0,0,0]
ax.YLabel.FontWeight='bold'
ax.YDir='reverse'


ax=axes('position',[0.35,0.15,.25,.1]);
x = 1:20;


Z=rand(size(x,1),size(x,2));
a=imagesc(x,x*0+1,Z);
grid on
ax.GridAlpha=1;
ax.LineWidth=1
ax.YTick=[]
ax.XTick=x+.5
ax.XTickLabel=[]
ax.XAxis.LineWidth=.001
ax.YAxis.LineWidth=.001
daspect([1.5,1,1]);
ax.TickLength=[0.0001,0.0001];%%
ax=axes('position',[0.35,0.05,.25,.1]);
x = 1:20;
y = 1:2;
[X,Y] = meshgrid(x,y);
C=ones([size(X,1),size(X,2),3]);

[R,TIEADJ] = tiedrank(Z);
a=imagesc(x,x*0+1,R)
grid on
ax.GridAlpha=1
ax.YTick=[]
ax.XTick=x+.5
ax.XTickLabel=[]
ax.XAxis.LineWidth=.001
ax.YAxis.LineWidth=.001
daspect([1.5,1,1]);
ax.TickLength=[0.0001,0.0001];
% 
ax=axes('position',[0.05,0.15,.25,.1]);
x = 1:20;
colors = cbrewer('seq', 'Greens', 128);
%colors = flipud(colors); % puts red on top, blue at the bottom


shg
Z=rand(size(x,1),size(x,2));
a=imagesc(x,x*0+1,Z);
colormap(ax,colors);
grid on
ax.GridAlpha=1;
ax.LineWidth=1
ax.YTick=[]
ax.XTick=x+.5
ax.XTickLabel=[]
ax.XAxis.LineWidth=.001
ax.YAxis.LineWidth=.001
daspect([1.5,1,1]);
ax.TickLength=[0.0001,0.0001];%%
ax=axes('position',[0.05,0.05,.25,.1]);
x = 1:20;


[R,TIEADJ] = tiedrank(Z);
a=imagesc(x,x*0+1,R)
colormap(ax,colors);
grid on
ax.GridAlpha=1
ax.YTick=[]
ax.XTick=x+.5
ax.XTickLabel=[]
ax.XAxis.LineWidth=.001
ax.YAxis.LineWidth=.001
daspect([1.5,1,1]);
ax.TickLength=[0.0001,0.0001];
%
ax=axes('position',[0.05,0.3,.25,.1]);
x = 1:20;
colors = cbrewer('seq', 'Reds', 128);
%colors = flipud(colors); % puts red on top, blue at the bottom


shg
Z=rand(size(x,1),size(x,2));
a=imagesc(x,x*0+1,Z);
colormap(ax,colors);
grid on
ax.GridAlpha=1;
ax.LineWidth=1
ax.YTick=[]
ax.XTick=x+.5
ax.XTickLabel=[]
ax.XAxis.LineWidth=.001
ax.YAxis.LineWidth=.001
daspect([1.5,1,1]);
ax.TickLength=[0.0001,0.0001];%%
ax=axes('position',[0.05,0.4,.25,.1]);
x = 1:20;


[R,TIEADJ] = tiedrank(Z);
a=imagesc(x,x*0+1,R)
colormap(ax,colors);
grid on
ax.GridAlpha=1
ax.YTick=[]
ax.XTick=x+.5
ax.XTickLabel=[]
ax.XAxis.LineWidth=.001
ax.YAxis.LineWidth=.001
daspect([1.5,1,1]);
ax.TickLength=[0.0001,0.0001];
%
% 
ax=axes('position',[0.7,0.1,.25,.3]);
AMC026_elec_fsaverage=load(['/Users/eghbalhosseiniasl1/MyData/ecog-sentence/subjects_raw/AMC026/IMAGING/','AMC026_elec_fsavg.mat']);
AMC026_elec_fsaverage=AMC026_elec_fsaverage.elec_fsavg_frs;
fspial_lh = ft_read_headshape('/Applications/freesurfer/subjects/fsaverage/surf/lh.pial');
fspial_lh.coordsys = 'fsaverage';
h=ft_plot_mesh(fspial_lh);
view([-100 10]);
material dull;
lighting gouraud;
camlight;
print(f, '-djpeg', strcat(analysis_path,'/ICA_schematic.jpeg'));
%%
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
hold on 
bl.LineWidth=2
bl1=bar(2,[B(1,2)],'DisplayName','F_2');
bl1.LineWidth=2
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
close all 
f=figure;
ax=axes('position',[0.15,0.15,.5,.7]);
x = 1:20;
y = 1:30;
[X,Y] = meshgrid(x,y);
C=ones([size(X,1),size(X,2),3]);
a=surf(X,Y,rand(size(X,1),size(X,2)));
a.FaceAlpha=1;
hold on
set(ax,'view',[0 90])
set(ax,'xlim',[0,20])
daspect([1.5,1.0,1]);
grid off
ax.YAxis.Visible='off';
ax.XAxis.Visible='off';
ax.ZAxis.Visible='off';
shg
colors = cbrewer('div', 'RdBu', 128);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);


print(f, '-djpeg', strcat(analysis_path,'/sentence_embedding.jpeg'));
    
%% 
close all 
f=figure;
ax=axes('position',[0.15,0.15,.5,.7]);
x = 1:2;
y = 1:20;
[X,Y] = meshgrid(x,y);
C=ones([size(X,1),size(X,2),3]);
a=surf(X,Y,0*rand(size(X,1),size(X,2)),C);
a.FaceAlpha=1;
hold on
set(ax,'view',[0 90])
set(ax,'xlim',[0,20])
daspect([1.5,1,1]);
grid off
ax.YAxis.Visible='off';
ax.XAxis.Visible='off';
ax.ZAxis.Visible='off';
shg
colors = cbrewer('div', 'RdBu', 128);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);


print(f, '-djpeg', strcat(analysis_path,'/sentence_embedding.jpeg'));
    
%% 
figure
AMC026_elec_fsaverage=load(['/Users/eghbalhosseiniasl1/MyData/ecog-sentence/subjects_raw/AMC026/IMAGING/','AMC026_elec_fsavg.mat']);
AMC026_elec_fsaverage=AMC026_elec_fsaverage.elec_fsavg_frs;
fspial_lh = ft_read_headshape('/Applications/freesurfer/subjects/fsaverage/surf/lh.pial');
fspial_lh.coordsys = 'fsaverage';
h=ft_plot_mesh(fspial_lh);
h.CData=[0*h.XData+1]';
view([-100 10]);
material dull;
lighting gouraud;
camlight;
print(f, '-djpeg', strcat(analysis_path,'/ICA_schematic.jpeg'));
