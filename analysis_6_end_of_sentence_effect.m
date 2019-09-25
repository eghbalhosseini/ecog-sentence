clear all 
close all 
home
%% 
if 1
    fprintf('adding basic ecog tools to path \n');
    addpath('~/MyCodes/basic-ecog-tools/');
    addpath('~/MyCodes/ecog-sentence/');
    addpath(genpath('~/MyCodes/basic-ecog-tools/activeBrain'));
    addpath(genpath('~/MyCodes/basic-ecog-tools/ecog-filters'));
    addpath(genpath('~/MyCodes/basic-ecog-tools/mex'));
end 
save_path='/Users/eghbalhosseiniasl1/MyData/ecog-sentence/';
%% 
data_path='/Users/eghbalhosseiniasl1/MyData/ecog-sentence';
analysis_path='/Users/eghbalhosseiniasl1/MyData/ecog-sentence/analysis/analysis_3_visualize_ch_response/';
subject_id='AMC026';
d= dir([data_path,strcat('/**/',subject_id,'*_crunched.mat')]);
fprintf(' %d .mat files were found \n', length(d));

% find the node information
d_nmerge= dir([data_path,'/**/*nmerge.txt']);
fprintf(' %d nmerge files were found \n', length(d_nmerge));
i=1 
% load node info
[node_cell,node_table]=generate_node_table_eh(strcat(d_nmerge(i).folder,'/',d_nmerge(i).name));
% 
plot_width=.9/length(d);
plot_length=.87;
%% there are 31 unique patterns, make a table of each sentence and for each word location what pattern it is member of 
open_nodes_cell=cellfun(@(x) x(1:8), node_table.open_nodes,'UniformOutput',false);
unique_open_node_patterns=unique(cell2mat(open_nodes_cell),'rows');
all_pattern_locations={};
for p=1:size(open_nodes_cell,1)
    open_node_pattern=open_nodes_cell{p,:};
    
    pattern_locations={};
    for i=1:8
       
        [~,pattern_memebership_id]=ismember(unique_open_node_patterns(:,1:i),open_node_pattern(1:i),'rows');
        pattern_locations=[pattern_locations,find(pattern_memebership_id)];
    end
    all_pattern_locations=[all_pattern_locations;[{open_node_pattern},pattern_locations]];
end
%% 
find_index= @(input_cell) find(cell2mat(cellfun(@(x) ~isempty(x), input_cell,'UniformOutput',false)));
all_example_sentence_locations=[];
all_sentence_pattern=[];
sentence_electrode_with_langauge_accross_sessions=[];
words_electrode_with_langauge_accross_sessions=[];

for i=1:length(d)
    example_sentence_locations=[];
    fprintf('adding %s from %s \n',d(i).name, strcat(d(i).folder,'/',d(i).name));
    subj=load(strcat(d(i).folder,'/',d(i).name));
    subj_id=fieldnames(subj);
    subj=subj.(subj_id{1});
    data=subj.data;
    info=subj.info;
    language_electrode=info.ramp_electrodes;
     
    language_electrode_num=find(language_electrode);
    % step 1: extract electrodes with siginificant language response
    sentence_trial_index=~cellfun(@isempty,cellfun(@ (x) strfind(x,'S'),info.word_type,'UniformOutput',false));
    % sentence
    sentences=[data{sentence_trial_index}];
    sentence_gamma_band_envelope=cellfun(@(x) x(1:10),{sentences.signal_hilbert_zs_downsample_parsed},'UniformOutput',false);
    % creat a cell with wordposition(row)*time in trial(column) structure
    sentence_gamma_band_envelope=[sentence_gamma_band_envelope{:,:}];
    sentence_gamma_band_envelope=cellfun(@transpose,sentence_gamma_band_envelope,'UniformOutput',false);
    % make a words positions*channel* trial tensor 
    sentence_gamma_band_ave_envelope_tensor=cell2mat(permute(sentence_gamma_band_envelope,[1,3,2]));
    %append to langauge_channel*trial*words positions
    sentence_gamma_band_ave_envelope_tensor=permute(sentence_gamma_band_ave_envelope_tensor,[2,3,1]);
    electrodes_with_language_response=sentence_gamma_band_ave_envelope_tensor(find(language_electrode),:,:);
    sentence_electrode_with_langauge_accross_sessions=cat(2,sentence_electrode_with_langauge_accross_sessions,electrodes_with_language_response);
    % 
    
    % find sentence locations 
    example_sentence={sentences(:).trial_string};
    example_sentence=cellfun(@(x) x(2:end),example_sentence,'UniformOutput',false);
    % find trial location in node cell 
    example_sentence_locations=cellfun(@(x) (regexpi([node_cell{:,1}],x)),example_sentence,'UniformOutput',false);
    example_sentence_locations=cell2mat(cellfun(@(x) find_index(x), example_sentence_locations, 'UniformOutput',false));
 
    all_sentence_pattern=[all_sentence_pattern;
        cell2mat(cellfun(@(x) x(1:8),node_table.open_nodes(example_sentence_locations),'UniformOutput',false))];
    %% 
    words_trial_index=~cellfun(@isempty,cellfun(@ (x) strfind(x,'W'),info.word_type,'UniformOutput',false));
    % words
    words=[data{words_trial_index}];
    words_gamma_band_envelope=cellfun(@(x) x(1:10),{words.signal_hilbert_zs_downsample_parsed},'UniformOutput',false);
    % creat a cell with wordposition(row)*time in trial(column) structure
    words_gamma_band_envelope=[words_gamma_band_envelope{:,:}];
    words_gamma_band_envelope=cellfun(@transpose,words_gamma_band_envelope,'UniformOutput',false);
    % make a words positions*channel* trial tensor
    words_gamma_band_ave_envelope_tensor=cell2mat(permute(words_gamma_band_envelope,[1,3,2]));
    %append to langauge_channel*trial*words positions
    words_gamma_band_ave_envelope_tensor=permute(words_gamma_band_ave_envelope_tensor,[2,3,1]);
    electrodes_with_language_response=words_gamma_band_ave_envelope_tensor(find(language_electrode),:,:);
    words_electrode_with_langauge_accross_sessions=cat(2,words_electrode_with_langauge_accross_sessions,electrodes_with_language_response);
    
    
    fprintf('added %s from %s \n',d(i).name, strcat(d(i).folder,'/',d(i).name));
    
end


%% plot example electrode for all conditions close all;


%% 
temp=[mean(words_electrode_with_langauge_accross_sessions,3),...
    mean(sentence_electrode_with_langauge_accross_sessions,3)];
min_ax=min(temp(:));
max_ax=max(temp(:));

a=mean(sentence_electrode_with_langauge_accross_sessions,2);
b=mean(a,3);
[~,idx]=sort(b);

word_samples=135;
post_trial_samples=75;
probe_samples=420;
x=[-8*word_samples+1:post_trial_samples+probe_samples];
y=1:size(sentence_electrode_with_langauge_accross_sessions,1);
% 
figure 
plot_width=.18;
plot_length=.8;
set(gcf,'OuterPosition',[-1615 441 694 1337]);
a=mean(words_electrode_with_langauge_accross_sessions,2);
axes('position',[.08,.05,plot_length,plot_width]);
colors = cbrewer('div', 'RdBu', 128);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);
imagesc(x,y,a(idx,:),[min_ax,max_ax]);
set(gca, 'ydir', 'normal','box','off','xtick',[-8*135:135:0,post_trial_samples,probe_samples],'xticklabel',[-8:0,1,2]);
title('words');
handles = colorbar;
handles.TickDirection = 'out';
handles.Box = 'off';
handles.Position=[0.9 0.05 .02 plot_width];
handles.Label.String = 'z-score';
drawnow;
xlabel('word position relative to last word'); ylabel('Electrode');
 

% 
axes('position',[.08,.1+plot_width,plot_length,plot_width]);
colors = cbrewer('div', 'RdBu', 128);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);
a=mean(sentence_electrode_with_langauge_accross_sessions,2);
imagesc(x,y,a(idx,:),[min_ax,max_ax]);
set(gca, 'ydir', 'normal','box','off','xtick',[-8*135:135:0,post_trial_samples,probe_samples],'xticklabel',[-8:0,1,2]);
title('sentences');
handles = colorbar;
handles.TickDirection = 'out';
handles.Box = 'off';
handles.Position=[0.9 .1+plot_width .02 plot_width];
handles.Label.String = 'z-score';
drawnow;
xlabel('word position relative to last word'); ylabel('Electrode');

axes('position',[.08,.15+2*plot_width,plot_length,plot_width]);

S=squeeze(double(mean(sentence_electrode_with_langauge_accross_sessions,2)));
W=squeeze(double(mean(words_electrode_with_langauge_accross_sessions,2)));
x=[-8*word_samples+1:post_trial_samples+probe_samples];
mean_electrode_pattern_response=nanmean(S,1);
std_electrode_pattern_response=nanstd(S,1)./sqrt(sum(~isnan(S),1));
bl = boundedline(x, mean_electrode_pattern_response, std_electrode_pattern_response, 'cmap', [1,0,0],'alpha');
        % make the sudo line for open node pregression
mean_electrode_pattern_response=nanmean(W,1);
std_electrode_pattern_response=nanstd(W,1)./sqrt(sum(~isnan(W),1));
        
bl = boundedline(x, mean_electrode_pattern_response, std_electrode_pattern_response, 'cmap', [.3,.3,1],'alpha');
set(gca,'ylim',[min_ax,max_ax],'xlim',[min(x),max(x)]);
set(gca, 'ydir', 'normal','box','off','xtick',[-8*135:135:0,post_trial_samples,probe_samples],'xticklabel',[-8:0,1,2]);

if ~exist(strcat('/Users/eghbalhosseiniasl1/MyData/ecog-sentence/end_of_sentence/',info.subject))
        mkdir(strcat('/Users/eghbalhosseiniasl1/MyData/ecog-sentence/end_of_sentence/',info.subject))
end 
print(gcf, '-depsc', strcat('/Users/eghbalhosseiniasl1/MyData/ecog-sentence/end_of_sentence/',info.subject,'/',info.subject,'_channel_',num2str(language_electrode_num(i)),'_response_end_of_sentence','.eps')); 

%% 

S=squeeze(double(mean(sentence_electrode_with_langauge_accross_sessions,2)));
W=squeeze(double(mean(words_electrode_with_langauge_accross_sessions,2)));
cval={'S' ,'W'};
cind=[ones(size(S,1),1);2*ones(size(W,1),1)];
c=cval(cind);
x=[-8*word_samples+1:post_trial_samples+probe_samples];
y=[S;W];

clear g
g=gramm('x',x,'y',y,'color',c);
g.set_color_options('map','d3_20');
g.geom_hline();
g.set_order_options('color',[1,2]);
g.stat_summary();
g.set_title('mean activity across all channels ');
g.axe_property('Xtick',[-8*135:135:0,post_trial_samples,probe_samples]);
g.axe_property('Xticklabel',[-8:0,1,2]);
g.set_names('x','word position','y','z-score','color','stimuli');
figure('Position',[100 100 800 550]);
g.draw();


g.export('file_name',strcat(info.subject,'end_of_sentence_condition_ave'),'export_path',strcat(save_path),'file_type','pdf');
%print(gcf, '-depsc', strcat(info.subject,'_condition_comparison_ave','.eps'));
%% sort trials based on value of last merge 

[last_val,idx]=sort(all_sentence_pattern(:,end));

S=squeeze(double(mean(sentence_electrode_with_langauge_accross_sessions,1)));
W=squeeze(double(mean(words_electrode_with_langauge_accross_sessions,1)));
S_sorted=S(idx,:);
W_sorted=W(idx,:);
mat2cell((unique(last_val)),length((unique(last_val))),1);
%cval={strsplit(transpose((num2str(unique(last_val)))),'')};
c=last_val;
x=[-8*word_samples+1:post_trial_samples+probe_samples];
y=S_sorted;
y1=W_sorted;

clear g
g(1,1)=gramm('x',x,'y',y,'color',c);
g(1,1).set_color_options('map','d3_20');
g(1,1).geom_hline();
%g.set_order_options('color',[1,2]);
g(1,1).stat_summary();
g(1,1).set_title('Sentences :mean activity across all channels ');
g(1,1).axe_property('Xtick',[-8*135:135:0,post_trial_samples,probe_samples]);
g(1,1).axe_property('Xticklabel',[-8:0,1,2]);
g(1,1).set_names('x','word position','y','z-score','color','num of open nodes');

g(2,1)=copy(g(1));
g(2,1)=gramm('x',x,'y',y1);
g(2,1).geom_hline();
%g.set_order_options('color',[1,2]);
g(2,1).stat_summary();
g(2,1).set_title('Words: mean activity across all channels ');
g(2,1).axe_property('Xtick',[-8*135:135:0,post_trial_samples,probe_samples]);
g(2,1).axe_property('Xticklabel',[-8:0,1,2]);
g(2,1).set_names('x','word position','y','z-score','color','num of open nodes');



figure('Position',[100 100 800 550]);
g.draw();
%g.export('file_name',strcat(info.subject,'_condition_comparison_ave'),'file_type','pdf');

g.export('file_name',strcat(info.subject,'end_of_sentence_condition_ave'),'export_path',strcat(save_path,'end_of_sentence/',info.subject,'/'),'file_type','pdf');


