analysis_path='/Users/eghbalhosseiniasl1/MyData/ECoG-sentence/analysis/swjn_schematic_pca_analysis/';
print_opt=0;
%%
x=[repmat(1:8,20,1)]'./8;
x=x(:);
x_ramp=x;
x=[x';.5*rand(10,size(x,1))];
w=rand(70,size(x,1));

output=.3*(w*x)+.1*rand(size(w*x));
output=output-mean(output,2)*ones(1,160);
[coeff,score,latent,tsquared,explained,mu]= pca(output);




%
close all 
f=figure;
aspect_ration=9.32./4.13;
y=500;
set(f,'position',[591 455 aspect_ration*y y]);
a=axes('position',[.07,.88,.25,.1]);
colors = cbrewer('div', 'RdYlBu', 8);
colors=[repmat(colors,20,1)];
trials=1:length(x_ramp);
s=scatter(trials,x_ramp,15,'filled');
a.YLim=[-.2,1.2];
s.CData=colors;
a.XAxis.Visible='on'
a.YAxis.Visible='on'
a.XTick=[]
a.YTick=[0.1250,1];
a.YTickLabel={'1','8'};
a.XAxis.Color=[1,1,1]
a.YLabel.String={'word','position'}
a.YLabel.HorizontalAlignment='right '
a.YLabel.VerticalAlignment='middle '
a.YLabel.Rotation=0;
a.XLabel.String='Sentences'
a.XLabel.Color=[0,0,0]

hold on 
guides=arrayfun(@(x) plot([x,x],a.YLim,'-','color',[.5,.5,.5]),trials(~mod(trials,8)))
arrayfun(@(x) uistack(x,'bottom'),guides)
a.FontSize=12;
a.FontWeight='bold '
%
a=axes('position',[.07,.5,.25,.3]);
colors = cbrewer('div', 'RdBu', 128);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);
imagesc(output)
set(gca, 'ydir', 'reverse','box','off');
a.XTick=[];
a.XLabel.String='Sentences';

a.YTick=[];
a.YLabel.String={'Electrodes'}
a.YLabel.HorizontalAlignment='right '
a.YLabel.VerticalAlignment='middle '
a.YLabel.Rotation=0;
%text(a.XLim(1)-5,a.XLim(1)-10,'A','HorizontalAlignment','right','FontSize',12,'FontWeight','bold')
a.FontSize=12
a.FontWeight='bold '
%
a=axes('position',[.07,.1,.25,.3])
colors = cbrewer('div', 'RdBu', 128);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);
imagesc(score)
set(gca, 'ydir', 'reverse','box','off');
a.XTick=[];
a.XLabel.String ='Components';
a.YTick=[];
a.YLabel.String='Weights';
a.YLabel.HorizontalAlignment='right '
a.YLabel.Rotation=0;
a.YTickLabelRotation=90;
%text(a.XLim(1)-1.5,a.XLim(1)-2,'B','HorizontalAlignment','right','FontSize',12,'FontWeight','bold');
a.FontSize=12;
a.FontWeight='bold '
%
%
a=axes('position',[.36,.1,.15,aspect_ration*.31])
colors = cbrewer('div', 'RdBu', 128);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);
imagesc(coeff)
set(gca, 'ydir', 'reverse','box','off');
xlims=a.XLim;
a.XTick=[];

a.XLabel.String='Components';

a.YTick=[];
a.YLabel.String='Sentences';
a.YLabel.HorizontalAlignment='right '
a.YLabel.Rotation=90;
a.FontSize=12;
a.FontWeight='bold '
a=axes('position',[.36,.81,.15,.08])
bl=bar(explained,'Facecolor','flat');
set(gca,'box','off');
bl.EdgeColor='none'
x_lim=get(gca,'xlim');
a.XAxis.Visible='off';
a.XLim=xlims;
a.YAxis.Color=[1,1,1]
a.YLabel.String='%var'
a.YLabel.Color=[0,0,0]
a.YLabel.Rotation=0;
a.FontWeight='bold '
a.YLabel.VerticalAlignment='middle'
%text(a.XLim(1)-5,a.YLim(2)+7,'C','HorizontalAlignment','right','FontSize',12,'FontWeight','bold')
a.FontSize=12;
% if ~exist(strcat(analysis_path))
%     mkdir(strcat(analysis_path));
% end
% if ~exist(strcat(analysis_path))
%     mkdir(strcat(analysis_path));
% end
% if print_opt==1
%     set(gcf,'PaperPosition',[.25 .25 8 6])
%     print(f,'-painters', '-djpeg', strcat(analysis_path,'/','swjn_schematic_pca_1.jpeg'));
% end
%
%
R_width=.25;
R_bar=.08;
num_rows=3;
plot_dist=.05;
R_height=(.9-.1-plot_dist*num_rows)/num_rows;
R_start=R_height*(0:num_rows-1)+plot_dist*(1:num_rows)+.05;
num_columns=1;
total_plots=length(R_start);
p=0;
x=[repmat(1:8,20,1)]';
x=x(:);
word_pos=x';
sort_idx={};
unique_word_pos=unique(word_pos);
for i=1:length(unique(word_pos))
    sort_idx=[sort_idx,find(word_pos==unique_word_pos(i))];
end
ax_min=floor(min(coeff(:)));
ax_max=ceil(max(coeff(:)));
for i=1:3
    R_indx=i-num_rows*fix((i-1)/num_rows);
    %f=figure(fix((i-1)/num_rows)+2);
    %set(gcf,'position',[29,12,852,1270])
    colors = cbrewer('div', 'RdYlBu', 8);
    ax=axes('position',[.58,R_start(R_indx),R_width,R_height]);
    A=cellfun(@(x) transpose(coeff(x,i)),sort_idx,'UniformOutput',false);
    R_inferred_sort=cellfun(@(x) sort(x),A,'UniformOutput',false);
    R_sort_ave=cellfun(@mean,A);
    bl=bar(cell2mat(R_inferred_sort),'Facecolor','flat');
    bar_color=colors(word_pos(cell2mat(sort_idx)),:);
    set(bl,'Linestyle','none');
    
    bl.CData=bar_color;
    set(gca,'box','off')
    
    %ax.YLim=[ax_min,ax_max];
    ax.XAxis.Visible='off';
    ax.YAxis.Color=[1,1,1 ];
    ax.YTick=[];
    if i==1
        ax.XAxis.Visible='on';
        ax.XTick=[];
    ax.XLabel.String={'Components arranged based ','on sorted word position'};
   
    end
    ax.YLabel.String=[sprintf('Component %d ',i)];
    ax.YLabel.Color=[0,0,0]
    if i==3
       % text(ax.XLim(1),ax.YLim(2)+.2,'D','HorizontalAlignment','right','FontSize',12,'FontWeight','bold')
    end 
    ax.FontSize=12;
    ax.FontWeight='bold '
    %
    ax1=axes('position',[.85,R_start(R_indx),R_bar,R_height]);
    hold on;
    bl1=arrayfun(@(x,y) bar(x,y,'Facecolor','flat'),unique_word_pos,R_sort_ave);
    set(bl1,'Linestyle','none');
    arrayfun(@(x,y) set(x,'CData',colors(y,:)),bl1,unique_word_pos);
    arrayfun(@(x,y) set(x,'Displayname',num2str(y)),bl1,unique_word_pos);
    
    set(gca,'box','off')
    ax1.YAxis.Visible='off'
    ax1.XAxis.Visible='off'
    if i==3
        %text(ax1.XLim(1),ax1.YLim(2)+.02,'E','HorizontalAlignment','right','FontSize',12,'FontWeight','bold')
    end 
    if i==1
        text(4,1.05*min(ax1.YLim),{'mean over','word position'},'HorizontalAlignment','center','FontSize',12,'VerticalAlignment','top');
        [h,icons,plots,legend_text]=legend('show','position',[[.85+1.1*R_bar,R_start(R_indx)+.1,.05,R_height/5]])
        arrayfun(@(x) set(set(icons(x).Children,'XData',icons(x).Children.XData./[.2;.2;1;1;.3])),[9:16])
        arrayfun(@(x) set(icons(x),'FontSize',10,'fontweight','bold'),[1:8])
        h.Box='off';
    end
        ax.FontSize=12;
    set(gca,'box','off')
    
            if ~mod(i,total_plots) | i==3
            %legend('show','Location','northeastoutside')
            p=p+1;
           
           
    end
end

a=axes('position',[.58,.88,.25,.1]);
colors = cbrewer('div', 'RdYlBu', 8);
colors=[repmat(colors,20,1)];
trials=1:length(x_ramp);
[x_ramp_sort,idx]=sort(x_ramp);
s=scatter(trials,x_ramp_sort,15,'filled');
hold on 
xpos=[1:20:160]+10;
ypos=[0:.125:1]-.1;
arrayfun(@(z) text(xpos(z),ypos(z),num2str(z),'fontsize',12,'FontWeight','bold'),1:8)
a.YLim=[-.2,1.2];
s.CData=colors(idx,:);
a.XAxis.Visible='off'
a.YAxis.Visible='on'
a.YTick=[];

a.YLabel.String={'Sorted word','position'}
a.YLabel.HorizontalAlignment='right '
a.YLabel.VerticalAlignment='middle'
a.YLabel.Rotation=0;

hold on 
a.FontSize=12;
a.FontWeight='bold '

% 
 if ~exist(strcat(analysis_path))
                mkdir(strcat(analysis_path))
            end
    
            %set(gcf,'PaperPosition',[.25 .25 8 6])
            print(f, '-djpeg', strcat(analysis_path,'/','ECoG_pca_sort_comp_graphics_presentation.jpeg'));