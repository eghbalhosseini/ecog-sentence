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
%
random_set=randperm(10000);
all_vectors_w2v=cell2mat(arrayfun(@(x) word2vec(emb,emb.Vocabulary(x)),...
    random_set,'UniformOutput',false)');

% 
vectors_1=all_vectors_w2v(randperm(10000,1000),:);
vectors_2=all_vectors_w2v(randperm(10000,1000),:);

%% 
corr_result={};
k = 100;
ColorTag = cbrewer('qual', 'Paired', k);
f=figure;
set(f,'position',[405 1252 900 900]);
    
vectors=vectors_1;
[idx_clust_1,~,~] = spectralcluster(vectors,k);

vectors=vectors_2;
[idx_clust_2,~,~] = spectralcluster(vectors,k);
    

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
        ax=axes('position',[0.06,0.19,.15,.29]);
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
        text(ax.XLim(1)-.5,ax.YLim(2)-4,'A','HorizontalAlignment','right','FontSize',20,'FontWeight','bold')
        %ar = annotation('arrow');
        %ar.Position=[.12,.25,.1,0.045];
        %ar = annotation('arrow');
        %ar.Position=[.125,.48,.1,-0.024];
        
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
        %text(ax.XLim(1),ax.YLim(1)+.15,'B','HorizontalAlignment','right','FontSize',20,'FontWeight','bold')
        ar = annotation('arrow');
        ar.Position=[.40,.25,.046,0];
        ar = annotation('arrow');
        ar.Position=[.413,.4,.033,0];
        
        % 
        vectors=vectors_1;
        [idx_clust,~,~] = spectralcluster(vectors,k);
        [~, SCORE] = pca(vectors);
        half_clust=(unique(idx_clust));
        half_clust=half_clust(randperm(k,floor(k/2)));
        b=arrayfun(@(x) find(idx_clust==x), half_clust, 'UniformOutput', false);
        idx_sel=cell2mat(cellfun(@(x) x(randi(length(x))),b,'UniformOutput',false));
        idx_sel1=idx_sel;
        %
        vectors=vectors_2;
        [idx_clust,~,~] = spectralcluster(vectors,k);
        [~, SCORE] = pca(vectors);
        half_clust=(unique(idx_clust));
        half_clust=half_clust(randperm(k,floor(k/2)));
        b=arrayfun(@(x) find(idx_clust==x), half_clust, 'UniformOutput', false);
        idx_sel=cell2mat(cellfun(@(x) x(randi(length(x))),b,'UniformOutput',false));
        idx_sel2=idx_sel;
        
        %%%%%%%%%%%%%%spectral clustering  
        ax=axes('position',[0.45,0.2,.1,.1]);
        h1=gscatterEH([SCORE(:,1),SCORE(:,2)],idx_clust,ColorTag);
        arrayfun(@(x) set(x,'MarkerFaceAlpha',.3), h1.Children)
        arrayfun(@(x) set(x,'MarkerEdgeColor','none'), h1.Children)
        hold on 
        h2=gscatterEH([SCORE(idx_sel1,1),SCORE(idx_sel1,2)],idx_clust(idx_sel1),ColorTag);
        arrayfun(@(x) set(x,'linewidth',1), h2.Children);
        ax.XTick=[];
        ax.YTick=[];
        text(ax.XLim(1),ax.YLim(2)+.15,'C','HorizontalAlignment','right','FontSize',20,'FontWeight','bold')
        ar = annotation('arrow');
        ar.Position=[.56,.25,.02,0];
        % 
        ax=axes('position',[0.45,0.35,.1,.1]);
        h1=gscatterEH([SCORE(:,1),SCORE(:,2)],idx_clust,ColorTag);
        arrayfun(@(x) set(x,'MarkerFaceAlpha',.3), h1.Children)
        arrayfun(@(x) set(x,'MarkerEdgeColor','none'), h1.Children)
        hold on 
        h2=gscatterEH([SCORE(idx_sel2,1),SCORE(idx_sel2,2)],idx_clust(idx_sel2),ColorTag);
        arrayfun(@(x) set(x,'linewidth',1), h2.Children);
        ax.XTick=[];
        ax.YTick=[];
        text(ax.XLim(1),ax.YLim(2)+.15,'B','HorizontalAlignment','right','FontSize',20,'FontWeight','bold')
        ar = annotation('arrow');
        ar.Position=[.56,.4,.02,0];
        %%%%%%%%%%%%%%%%% correlation 
        
        ax=axes('position',[0.57,0.24,.2,.2]);
        vector_1_spec=[vectors_1(idx_sel1,:);vectors_1(idx_sel2,:)];
        vector_2_spec=[vectors_2(idx_sel1,:);vectors_2(idx_sel2,:)];
        similarity_spec_1=squareform(pdist(vector_1_spec,'correlation'));
        similarity_spec_2=squareform(pdist(vector_2_spec,'correlation'));
        x = 1:size(vector_1_spec,1);
        [X,Y] = meshgrid(x,x);
        plot3([0,100,100,0,0],[0,0,100,100,0],[2,2,2,2,2],'Color','k');
        hold on ;
        plot3([0,100,100,0,0]+5,[0,0,100,100,0]-5,-1+[0,0,0,0,0],'Color','k');
        hold on ;
        a=surf(X,Y,similarity_spec_1);
        colors = cbrewer('div', 'RdBu', 128);
        colors = flipud(colors); % puts red on top, blue at the bottom
        colormap(colors);
        
        a.EdgeColor='none';
        set(gca, 'ydir', 'reverse','box','off');
        set(ax,'view',[0,90]);
        hold on ;
        a=surf(X+5,Y-5,similarity_spec_2);
        a.EdgeColor='none';
        daspect([1,1,1]);
        ax.YAxis.Visible='off';
        ax.XAxis.Visible='off';
        %text(ax.XLim(1),ax.YLim(1)+.15,'B','HorizontalAlignment','right','FontSize',20,'FontWeight','bold')
        ar = annotation('arrow');
        ar.Position=[.40,.25,.046,0];
        ar = annotation('arrow');
        ar.Position=[.413,.4,.033,0];

        % 
       text(ax.XLim(1)-.5,ax.YLim(1)-5,'D','HorizontalAlignment','right','FontSize',20,'FontWeight','bold')
    end 
    
end 
%
corr_=sort([corr_result{:,3}]);
ax=axes('position',[0.76,0.20,.05,.23]);


a=scatter(corr_,1:length(corr_));
xLims=ax.XLim;
ax.XLim=[-1,1]*1.05*max(xLims);
ax.YLim=[-5,105];
a.SizeData=15;
a.MarkerEdgeColor=[.3,.3,.3];
set(gca, 'ydir', 'normal','box','off');

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
ax.Title.String='Corr';
text(ax.XLim(1),ax.YLim(2)+10,'E','HorizontalAlignment','right','FontSize',20,'FontWeight','bold')

        if ~exist(strcat(analysis_path))
            mkdir(strcat(analysis_path))
        end 

print(f, '-djpeg', strcat(analysis_path,'/sentence_selection.jpeg'));

%%
