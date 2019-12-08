analysis_path='/Users/eghbalhosseiniasl1/MyData/ECoG-sentence/analysis/swjn_schematic_pca_analysis/';
print_opt=1;
%%
x=[repmat(1:8,20,1)]'./8;
x=x(:);
x=[x';.5*rand(10,size(x,1))];
w=rand(70,size(x,1));

output=.3*(w*x)+.1*rand(size(w*x));
output=output-mean(output,2)*ones(1,160);
[coeff,score,latent,tsquared,explained,mu]= pca(output);




%%
f=figure;
set(gcf,'position',[1376 353 988 992]);
a=axes('position',[.1,.25,.2,.1]);
colors = cbrewer('div', 'RdBu', 128);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);
imagesc(output)
set(gca, 'ydir', 'reverse','box','off');
a.XTick=80;
a.XTickLabel='Words';
a.YTick=35;
a.YTickLabel='Electrode';
a.YTickLabelRotation=90;
text(a.XLim(1)-5,a.XLim(1)-10,'A','HorizontalAlignment','right','FontSize',20,'FontWeight','bold')
a.FontSize=16
%
a=axes('position',[.1,.1,.2,.1])
colors = cbrewer('div', 'RdBu', 128);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);
imagesc(score)
set(gca, 'ydir', 'reverse','box','off');
a.XTick=35;
a.XTickLabel='Component';
a.YTick=35;
a.YTickLabel='Weight';
a.YTickLabelRotation=90;
text(a.XLim(1)-1.5,a.XLim(1)-2,'B','HorizontalAlignment','right','FontSize',20,'FontWeight','bold');
a.FontSize=16;
%
a=axes('position',[.35,.1,.15,.2])
colors = cbrewer('div', 'RdBu', 128);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);
imagesc(coeff)
set(gca, 'ydir', 'reverse','box','off');
xlims=a.XLim;
a.XTick=35;
a.XTickLabel='Component';
a.YTick=70;
a.YTickLabel='Words';
a.YTickLabelRotation=90;
a.FontSize=16;
a=axes('position',[.35,.3,.15,.05])
bl=bar(explained,'Facecolor','flat');
set(gca,'box','off');
bl.EdgeColor='none'
x_lim=get(gca,'xlim');
a.XAxis.Visible='off';
a.XLim=xlims
a.YAxis.Visible='off'
text(a.XLim(1)-5,a.YLim(2)+7,'C','HorizontalAlignment','right','FontSize',20,'FontWeight','bold')
title('%var');
a.FontSize=16;
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

R_width=.13;
R_bar=.05;
num_rows=5;
plot_dist=.03;
R_height=(.55-.1-plot_dist*num_rows)/num_rows;
R_start=R_height*(0:num_rows-1)+plot_dist*(1:num_rows)+.07;
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
    ax=axes('position',[.51,R_start(R_indx),R_width,R_height]);
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
    ax.YAxis.Visible='off';
    if i==1
        ax.XAxis.Visible='on';
        ax.XTick=80;
    ax.XTickLabel='Sorted';
   
    end
    ax.Title.String=[sprintf('component %d ',i)];
    if i==3
        text(ax.XLim(1),ax.YLim(2)+.2,'D','HorizontalAlignment','right','FontSize',20,'FontWeight','bold')
    end 
    ax.FontSize=16;
    %
    ax1=axes('position',[.67,R_start(R_indx),R_bar,R_height]);
    hold on;
    bl1=arrayfun(@(x,y) bar(x,y,'Facecolor','flat'),unique_word_pos,R_sort_ave);
    set(bl1,'Linestyle','none');
    arrayfun(@(x,y) set(x,'CData',colors(y,:)),bl1,unique_word_pos);
    arrayfun(@(x,y) set(x,'Displayname',num2str(y)),bl1,unique_word_pos);
    
    set(gca,'box','off')
    ax1.YAxis.Visible='off'
    ax1.XAxis.Visible='off'
    if i==3
        text(ax1.XLim(1),ax1.YLim(2)+.02,'E','HorizontalAlignment','right','FontSize',20,'FontWeight','bold')
    end 
    if i==1
        text(4,1.05*min(ax1.YLim),'mean','HorizontalAlignment','center','FontSize',16,'VerticalAlignment','top');
        legend('show','position',[[.62+R_width,R_start(R_indx)+.1,.005,R_height/5]])
    end
        ax.FontSize=16;
    set(gca,'box','off')
    
            if ~mod(i,total_plots) | i==3
            %legend('show','Location','northeastoutside')
            p=p+1;
           
            if ~exist(strcat(analysis_path))
                mkdir(strcat(analysis_path))
            end
    
            %set(gcf,'PaperPosition',[.25 .25 8 6])
            print(f, '-painters','-djpeg', strcat(analysis_path,'/','ECoG_pca_sort_comp_graphics.jpeg'));
    end
end
