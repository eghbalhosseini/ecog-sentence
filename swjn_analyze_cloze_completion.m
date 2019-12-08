%% STEP 0: prepare the workspace
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
data_path='/Users/eghbalhosseiniasl1/MyData/ecog-sentence/cloze_data';

d_cum_cloze= dir([data_path,strcat('/*.csv')]);
fprintf(' %d cloze files were found \n', length(d_cum_cloze));
d_paragraph= dir([data_path,strcat('/*.txt')]);
fprintf(' %d paragraph files were found \n', length(d_paragraph));



save_path='/Users/eghbalhosseiniasl1/MyData/ECoG-sentence/cloze_data/';
analysis_path='/Users/eghbalhosseiniasl1/MyData/ECoG-sentence/analysis/swjn_analyze_cloze_completion';
%
if 1
    fprintf('adding basic ecog tools to path \n');
    addpath('~/MyCodes/basic-ecog-tools/');
    addpath(genpath('~/MyCodes/ecog-sentence/'));
    addpath(genpath('~/MyCodes/basic-ecog-tools/activeBrain'));
    addpath(genpath('~/MyCodes/basic-ecog-tools/ecog-filters'));
    addpath(genpath('~/MyCodes/basic-ecog-tools/mex'));
end


find_index= @(input_cell) find(cell2mat(cellfun(@(x) ~isempty(x), input_cell,'UniformOutput',false)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  STEP 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cloze_data=readtable(strcat(d_cum_cloze(1).folder,'/',d_cum_cloze(1).name))
%% 
all_paras={};
fileID = fopen(strcat(d_paragraph(1).folder,'/',d_paragraph(1).name),'r','n','UTF-8');
tline = 'Paragraph'
while ischar(tline)
    
    tline = fgetl(fileID);
    disp(tline)
    if isempty(strfind(tline,'Paragraph')) && ischar(tline)
    all_paras=[all_paras;(tline)];
    end 
end
fclose(fileID);
%%  
E=strjoin(cloze_data{:,2});
E=[cloze_data{:,2}];
all_idx={};
all_firsts={};
all_first_cloze_prob={};
for k=1:size(all_paras,1)
   B = convertCharsToStrings(all_paras{k});
   C=strsplit(B,'.')';
   D=strsplit(C(1),' ');
   first_sent=strjoin(D(2:end));
   first_sent=erase(first_sent,',');
   first_sent=erase(first_sent,"?");
   first_sent=erase(first_sent,"-");
   first_sent=strsplit(first_sent)';
   % 
   E=cloze_data{cloze_data{:,1}==k,2};
   probs=cloze_data{cloze_data{:,1}==k,3};
   all_firsts=[all_firsts;{first_sent'}];
   all_first_cloze_prob=[all_first_cloze_prob;{probs(1:length(first_sent))}];
   
end
% 
[sent_len,idx]=sort(cell2mat(cellfun(@(x) size(x,1),all_first_cloze_prob,'UniformOutput',false)))

longest=max(max(cell2mat(cellfun(@size,all_first_cloze_prob,'UniformOutput',false))));
all_cloze_mat=cell2mat(cellfun(@(x)[x;nan*ones(longest-size(x,1),1)]',all_first_cloze_prob(idx),'UniformOutput',false));

%% 
close all 
f=figure;
set(f,'position',[[134 278 1426 1060]])
ax=axes('position',[.1,.05,.5,.6])
x=1+(1:size(all_cloze_mat,2));
y=1:size(all_cloze_mat,1);
colors = cbrewer('div', 'RdBu', 128);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);
im=imagesc(all_cloze_mat);
hold on 
imagesc([40,41],1:40,nan)
%    handles = colorbar;
%    handles.TickDirection = 'out';
%    handles.Box = 'off';
%    handle.Ylabel.String='probability'
%    drawnow;
set(gca, 'ydir', 'normal','box','off');
hold on 
arrayfun(@(x,y) text(x+1,y,strcat(':  ',num2str(x+.5),'  sentence length'),'HorizontalAlignment','left','Color',[1,1,1],'fontsize',14,'FontWeight','bold'),sent_len(1)-.5,1);
arrayfun(@(x,y) text(x+1,y,strcat(': ',num2str(x+.5)),'HorizontalAlignment','left','Color',[1,1,1],'fontsize',14,'FontWeight','bold'),sent_len(2:end)-.5,[2:40]');
ax.XLim=[0,42.5]
ax.XAxis.Visible='off'
ax.YLabel.String='first sentence in paragraph'
ax.XLabel.String='sentence length'


%
ax.FontSize=14;
% 
ax=axes('position',[.1,.7,.5,.25])
y1=double(nanmean(all_cloze_mat,1));
x=1:length(y1);

e1=double(nanstd(all_cloze_mat,1)./sqrt(sum(~isnan(all_cloze_mat),1)));
[l,p] = boundedline(x, y1, e1, '-b');
hAnnotation=arrayfun(@(x) get(x,'Annotation'),p);
hLegendEntry = arrayfun(@(x) get(x,'LegendInformation'),hAnnotation);
arrayfun(@(x) set(x,'IconDisplayStyle','off'),hLegendEntry);

l.LineWidth=2;
hold on 
ax.YLim=[-.05,.8]
ax.XLim=[0,42.5]
ax.YTickLabel
ax.YLabel.String='probably over all sentences';
ax.Title.String='Cloze probability for first sentence in each paragraph'
ax.FontSize=14;
%
if ~exist(analysis_path)
    mkdir(analysis_path)
end
print(f,'-djpeg', strcat(analysis_path,'/','cloze_results'));

