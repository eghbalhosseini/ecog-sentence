%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 0: prepare the workspace 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all 
close all 
data_path='/Users/eghbalhosseiniasl1/MyData/ecog-sentence';
subject_id='AMC026';
d= dir([data_path,strcat('/**/',subject_id,'*_crunched.mat')]);
fprintf(' %d .mat files were found \n', length(d));
gamma_band_index=4;
save_path='/Users/eghbalhosseiniasl1/MyData/ECoG-sentence/crunched/';
analysis_path='/Users/eghbalhosseiniasl1/MyData/ECoG-sentence/analysis/swjn_analysis_1_2_relate_gamma_to_syn_word_pos_open_node_effect/';
%
if 1
    fprintf('adding basic ecog tools to path \n');
    addpath('~/MyCodes/basic-ecog-tools/');
    addpath('~/MyCodes/ecog-sentence/');
    addpath(genpath('~/MyCodes/basic-ecog-tools/activeBrain'));
    addpath(genpath('~/MyCodes/basic-ecog-tools/ecog-filters'));
    addpath(genpath('~/MyCodes/basic-ecog-tools/mex'));
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find the node information 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d_nmerge= dir([data_path,'/**/*nmerge.txt']);
fprintf(' %d nmerge files were found \n', length(d_nmerge));
i=1;
[node_cell,node_table]=generate_node_table_eh(strcat(d_nmerge(i).folder,'/',d_nmerge(i).name));

sentences_with_8_node_operations=cellfun(@(x) length(x)==8, node_table.open_nodes,'UniformOutput',false);
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

unique_pattern_locations={};
for p=1:size(unique_open_node_patterns,1)
    unique_open_node_pattern=unique_open_node_patterns(p,:);
    pattern_locations={};
    for i=1:8
        a=unique_open_node_pattern(1:i);
        pattern_locations=[pattern_locations,
        find(~cell2mat(cellfun(@isempty,cellfun(@(x) strfind(num2str(x),num2str(a)),open_nodes_cell,'UniformOutput',false),'UniformOutput',false)))];
    end
    unique_pattern_locations=[unique_pattern_locations;[unique_open_node_pattern,transpose(pattern_locations)]];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  STEP 1 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
find_index= @(input_cell) find(cell2mat(cellfun(@(x) ~isempty(x), input_cell,'UniformOutput',false)));
all_sentence_pattern=[];
all_closing_pattern=[];
all_sentence_pattern_id=[];
sentence_electrode_with_langauge_accross_sessions=[];
words_electrode_with_langauge_accross_sessions=[];
session_sentence_hilbert_band_envelope_lang_elec_tensor=[];
session_sentence_hilbert_zs_band_envelope_lang_elec_tensor=[];
session_nonwords_hilbert_band_envelope_lang_elec_tensor=[];
session_words_hilbert_band_envelope_lang_elec_tensor=[];
session_jabberwocky_hilbert_band_envelope_lang_elec_tensor=[];

for i=1:length(d)
    subj=load(strcat(d(i).folder,'/',d(i).name));
    subj_id=fieldnames(subj);
    subj=subj.(subj_id{1});
    data=subj.data;
    info=subj.info;
    language_electrode=info.language_responsive_electrodes_hilbert_odd;
    language_electrode_num=find(language_electrode);
    % step 1: extract electrodes with siginificant language response
    sentence_trial_index=~cellfun(@isempty,cellfun(@ (x) strfind(x,'S'),info.word_type,'UniformOutput',false));
    %%%%%%%%%%%%%%% sentence
    sentences=[data{sentence_trial_index}];% 
    word_length=sentences(1).signal_range_downsample(1,2)-sentences(1).signal_range_downsample(1,1)+1;
    % creat a cell with wordposition(row)*time in trial(column) structure
    %  
    hilbert_band_envelope=cellfun(@(x) x(1:8),{sentences.signal_hilbert_downsample_parsed},'UniformOutput',false);
    hilbert_band_envelope=[hilbert_band_envelope{:,:}];
    hilbert_band_envelope=cellfun(@transpose,hilbert_band_envelope,'UniformOutput',false);
    % make a words positions*channel* trial tensor 
    hilbert_band_envelope_tensor=cell2mat(permute(hilbert_band_envelope,[1,3,2]));
    %append to langauge_channel*trial*words positions
    hilbert_band_envelope_tensor=permute(hilbert_band_envelope_tensor,[2,3,1]);
    hilbert_band_envelope_lang_elec_tensor=hilbert_band_envelope_tensor(find(language_electrode),:,:);

    session_sentence_hilbert_band_envelope_lang_elec_tensor=cat(2,session_sentence_hilbert_band_envelope_lang_elec_tensor,hilbert_band_envelope_lang_elec_tensor);
    % find sentence locations 
    example_sentence=cellfun(@(x) x(2:end),{sentences(:).trial_string},'UniformOutput',false);
    example_sentence_locations=cellfun(@(x) (regexpi([node_cell{:,1}],x)),example_sentence,'UniformOutput',false);
    example_sentence_locations=cell2mat(cellfun(@(x) find_index(x), example_sentence_locations, 'UniformOutput',false));
    all_sentence_pattern=[all_sentence_pattern;
        cell2mat(cellfun(@(x) x(1:8),node_table.open_nodes(example_sentence_locations),'UniformOutput',false))];
    all_sentence_pattern_id=[all_sentence_pattern_id;cell2mat(cellfun(@(x) find(ismember(unique_open_node_patterns,x,'rows')),...
        cellfun(@(x) x(1:8),node_table.open_nodes(example_sentence_locations),'UniformOutput',false),'UniformOutput',false))];
    all_closing_pattern=[all_closing_pattern;
        cell2mat(cellfun(@(x) x(1:8),node_table.closed_nodes(example_sentence_locations),'UniformOutput',false))];
    %%%%%%%%%%%%%%%%% words
    words_trial_index=~cellfun(@isempty,cellfun(@ (x) strfind(x,'W'),info.word_type,'UniformOutput',false));
    % words
    words=[data{words_trial_index}];
    hilbert_band_envelope=cellfun(@(x) x(1:8),{words.signal_hilbert_downsample_parsed},'UniformOutput',false);
    % creat a cell with wordposition(row)*time in trial(column) structure
    hilbert_band_envelope=[hilbert_band_envelope{:,:}];
    hilbert_band_envelope=cellfun(@transpose,hilbert_band_envelope,'UniformOutput',false);
    % make a words positions*channel* trial tensor
    hilbert_band_envelope_tensor=cell2mat(permute(hilbert_band_envelope,[1,3,2]));
    %append to langauge_channel*trial*words positions
    hilbert_band_envelope_tensor=permute(hilbert_band_envelope_tensor,[2,3,1]);
    hilbert_band_envelope_lang_elec_tensor=hilbert_band_envelope_tensor(find(language_electrode),:,:);
    session_words_hilbert_band_envelope_lang_elec_tensor=cat(2,session_words_hilbert_band_envelope_lang_elec_tensor,hilbert_band_envelope_lang_elec_tensor);
    %%%%%%%%%%%%%%%%% Nonwords
    nonwords_trial_index=~cellfun(@isempty,cellfun(@ (x) strfind(x,'N'),info.word_type,'UniformOutput',false));
    % nonwords
    nonwords=[data{nonwords_trial_index}];

    
    hilbert_band_ave_envelope_pre_trial=cellfun(@(x) reshape(x,size(x,1),1,size(x,2)),{nonwords.signal_ave_pre_trial_hilbert_downsample},'UniformOutput',false);
    hilbert_band_ave_envelope_pre_trial=[hilbert_band_ave_envelope_pre_trial{:,:}];
    hilbert_band_envelope=cellfun(@(x) x(1:8),{nonwords.signal_hilbert_downsample_parsed},'UniformOutput',false);
    % creat a cell with wordposition(row)*time in trial(column) structure
    hilbert_band_envelope=[hilbert_band_envelope{:,:}];
    hilbert_band_envelope=cellfun(@transpose,hilbert_band_envelope,'UniformOutput',false);
    % make a nonwords positions*channel* trial tensor
    hilbert_band_envelope_tensor=cell2mat(permute(hilbert_band_envelope,[1,3,2]));
    %append to langauge_channel*trial*nonwords positions
    hilbert_band_envelope_tensor=permute(hilbert_band_envelope_tensor,[2,3,1]);
    hilbert_band_envelope_lang_elec_tensor=hilbert_band_envelope_tensor(find(language_electrode),:,:);
    session_nonwords_hilbert_band_envelope_lang_elec_tensor=cat(2,session_nonwords_hilbert_band_envelope_lang_elec_tensor,hilbert_band_envelope_lang_elec_tensor);
    %%%%%%%%%%%%%%%%%Jabberwocky
    jabberwocky_trial_index=~cellfun(@isempty,cellfun(@ (x) strfind(x,'J'),info.word_type,'UniformOutput',false));
    jabberwocky=[data{jabberwocky_trial_index}];

    hilbert_band_ave_envelope_pre_trial=cellfun(@(x) reshape(x,size(x,1),1,size(x,2)),{jabberwocky.signal_ave_pre_trial_hilbert_downsample},'UniformOutput',false);
    hilbert_band_ave_envelope_pre_trial=[hilbert_band_ave_envelope_pre_trial{:,:}];
    hilbert_band_envelope=cellfun(@(x) x(1:8),{jabberwocky.signal_hilbert_downsample_parsed},'UniformOutput',false);
    % creat a cell with wordposition(row)*time in trial(column) structure
    hilbert_band_envelope=[hilbert_band_envelope{:,:}];
    hilbert_band_envelope=cellfun(@transpose,hilbert_band_envelope,'UniformOutput',false);
    % make a nonwords positions*channel* trial tensor
    hilbert_band_envelope_tensor=cell2mat(permute(hilbert_band_envelope,[1,3,2]));
    %append to langauge_channel*trial*nonwords positions
    hilbert_band_envelope_tensor=permute(hilbert_band_envelope_tensor,[2,3,1]);
    hilbert_band_envelope_lang_elec_tensor=hilbert_band_envelope_tensor(find(language_electrode),:,:);
    session_jabberwocky_hilbert_band_envelope_lang_elec_tensor=cat(2,session_jabberwocky_hilbert_band_envelope_lang_elec_tensor,hilbert_band_envelope_lang_elec_tensor);
    
    fprintf('added %s from %s \n',d(i).name, strcat(d(i).folder,'/',d(i).name));
    clear data subj
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  compute correlation of actvitiy with the number of open nodes 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all 

num_rows=4;
num_columns=2;
total_plots=num_rows*num_columns;
look_out_window=100; %200 ms
offset=30; % 100 ms
p_val_threshold=0.05;
word_range=[3:8];
j=0;
[row,column]=ind2sub(size(ones(2,3)),find(ones(2,3)));
for i=1:size(session_sentence_hilbert_band_envelope_lang_elec_tensor,1)
    electrode_response=squeeze(session_sentence_hilbert_band_envelope_lang_elec_tensor(i,:,:));    
  
    electrode_response_sentence_cell=mat2cell(electrode_response,ones(1,size(electrode_response,1)),word_length*ones(1,8));
    electrode_post_word_sentence_mean=double(cell2mat(cellfun(@(x) nanmean(x(offset+(1:look_out_window))),electrode_response_sentence_cell,'UniformOutput',false)));
    word_position=repmat(1:size(electrode_post_word_sentence_mean,2),[size(electrode_post_word_sentence_mean,1),1]);
    % 
    electrode_post_word_sentence_mean_cell=cellfun(@transpose,mat2cell(electrode_post_word_sentence_mean,ones(1,size(electrode_post_word_sentence_mean,1)),size(electrode_post_word_sentence_mean,2)),'UniformOutput',false);
    all_sentence_pattern_cell=cellfun(@transpose,mat2cell(all_sentence_pattern,ones(1,size(all_sentence_pattern,1)),size(all_sentence_pattern,2)),'uniformoutput',false);
    word_position_cell=cellfun(@transpose,mat2cell(word_position,ones(1,size(word_position,1)),size(word_position,2)),'uniformoutput',false);
    % select a subset of words 
    electrode_post_word_sentence_mean_cell=cellfun(@(x) x(word_range),electrode_post_word_sentence_mean_cell,'UniformOutput',false);
    all_sentence_pattern_cell=cellfun(@(x) x(word_range),all_sentence_pattern_cell,'UniformOutput',false);
    word_position_cell=cellfun(@(x) x(word_range),word_position_cell,'UniformOutput',false);
    
    [rho,p_val]=cellfun(@(x,y) corr(x,y,'type','Spearman'),electrode_post_word_sentence_mean_cell,all_sentence_pattern_cell,'UniformOutput',false);
    elec_response_to_open_node_corr=[rho{:}];
    elec_response_to_open_node_corr_p_val=[p_val{:}];
    elec_response_to_open_node_corr_significant=elec_response_to_open_node_corr_p_val<p_val_threshold;
    
    % 
    [rho,p_val]=cellfun(@(x,y) corr(x,y,'type','Spearman'),electrode_post_word_sentence_mean_cell,word_position_cell,'UniformOutput',false);
    elec_response_to_word_position_corr=[rho{:}];
    elec_response_to_word_position_corr_p_val=[p_val{:}];
    elec_response_to_word_position_corr_significant=elec_response_to_word_position_corr_p_val<p_val_threshold;
    %
    [rho,p_val]=cellfun(@(x,y) corr(x,y,'type','Spearman'),all_sentence_pattern_cell,word_position_cell,'UniformOutput',false);
    corr_range=([rho{:}]-min([rho{:}]))./(max([rho{:}])-min([rho{:}]));
    color_range=transpose(([1,1,1]')*corr_range);
    % 
    significance_overlap= 1 |elec_response_to_open_node_corr_significant | elec_response_to_word_position_corr_significant;
    positive_rhos= elec_response_to_word_position_corr>0 & elec_response_to_open_node_corr>0;
    elec_response_to_open_node_corr=elec_response_to_open_node_corr(positive_rhos);
    elec_response_to_word_position_corr=elec_response_to_word_position_corr(positive_rhos);
    color_range=color_range(positive_rhos,:);
    % 
    figure(fix((i-1)/total_plots)+1);
    set(gcf,'position',[1451,-285,957,1342]);
    sub_title=sprintf('subj: %s,\n electrodes: %d ', info.subject ,language_electrode_num(i));
    ax=subplot(num_rows,num_columns,i-total_plots*fix((i-1)/total_plots));
    figures=scatter_corner_hist(elec_response_to_word_position_corr,elec_response_to_open_node_corr,'word corr','open node corr',color_range);
    sub_title=sprintf('subj: %s,\n electrodes: %d ', info.subject ,language_electrode_num(i));
    figures(1).Title=title(sub_title);
    figures(2).FontSize=10;
    ~mod(i,total_plots)
    if ~mod(i,total_plots) | i==size(session_sentence_hilbert_band_envelope_lang_elec_tensor,1)
    if ~exist(strcat(analysis_path,info.subject))
        mkdir(strcat(analysis_path,info.subject))
    end
    print(gcf, '-depsc', strcat(analysis_path,info.subject,'/',info.subject,'_ch_',num2str(language_electrode_num(i)),'_corr_open_node_vs_word_pos.eps')); 
    close(gcf)
    end 
end 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  compute correlation based on unique patterns by taking the mean of values   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all 
unique_sentence_pattern=mat2cell(unique(all_sentence_pattern,'row'),ones(1,size(unique(all_sentence_pattern,'row'),1)),size(all_sentence_pattern,2));

uniq_sent_patt_loc_in_all_sentences=cellfun(@(x) find(ismember(all_sentence_pattern,x,'row')),unique_sentence_pattern,'UniformOutput',false);
look_out_window=100; %200 ms
offset=15; % 100 ms
word_range=[1:8];
num_rows=4;
num_columns=2;
Beta_mat_unique=[];
[row,column]=ind2sub(size(ones(2,3)),find(ones(2,3)));
p_val_threshold=0.05;
total_plots=num_rows*num_columns;

j=0;
for i=1:length(language_electrode_num)
    elec_response_to_open_node_pattern=[];
    electrode_post_word_sentence_mean_cell=[];
    elec_response_to_open_node_corr=[];
    elec_response_to_word_position_corr=[];
    electrode_response=squeeze(session_sentence_hilbert_band_envelope_lang_elec_tensor(i,:,:));
    electrode_response_sentence_cell=mat2cell(electrode_response,ones(1,size(electrode_response,1)),word_length*ones(1,8));
    electrode_post_word_sentence_mean=cell2mat(cellfun(@(x) nanmean(x(offset+(1:look_out_window))),electrode_response_sentence_cell,'UniformOutput',false));
    for k=1:size(unique_sentence_pattern,1)
        pattern_response=mean(electrode_post_word_sentence_mean(uniq_sent_patt_loc_in_all_sentences{k},:),1);
        elec_response_to_open_node_pattern=[elec_response_to_open_node_pattern;pattern_response];
    end
    electrode_post_word_sentence_mean_cell=cellfun(@transpose,mat2cell(elec_response_to_open_node_pattern,ones(1,size(elec_response_to_open_node_pattern,1)),size(elec_response_to_open_node_pattern,2)),'UniformOutput',false);
    unique_sentence_pattern_cell=cellfun(@transpose,unique_sentence_pattern,'uniformoutput',false);
    word_position=repmat(1:size(elec_response_to_open_node_pattern,2),[size(elec_response_to_open_node_pattern,1),1]);
    word_position_cell=cellfun(@transpose,mat2cell(word_position,ones(1,size(word_position,1)),size(word_position,2)),'uniformoutput',false);
    % select a subset of words
    electrode_post_word_sentence_mean_cell=cellfun(@(x) x(word_range),electrode_post_word_sentence_mean_cell,'UniformOutput',false);
    unique_sentence_pattern_cell=cellfun(@(x) x(word_range),unique_sentence_pattern_cell,'UniformOutput',false);
    word_position_cell=cellfun(@(x) x(word_range),word_position_cell,'UniformOutput',false);
    % compute correlation for open node
    [rho,p_val]=cellfun(@(x,y) corr(x,y,'type','Spearman'),electrode_post_word_sentence_mean_cell,unique_sentence_pattern_cell,'UniformOutput',false);
    elec_response_to_open_node_corr=[rho{:}];
    elec_response_to_open_node_corr_p_val=[p_val{:}];
    elec_response_to_open_node_corr_significant=elec_response_to_open_node_corr_p_val<p_val_threshold;
    % compute correlation for word position 
    [rho,p_val]=cellfun(@(x,y) corr(x,y,'type','Spearman'),electrode_post_word_sentence_mean_cell,word_position_cell,'UniformOutput',false);
    elec_response_to_word_position_corr=[rho{:}];
    elec_response_to_word_position_corr_p_val=[p_val{:}];
    elec_response_to_word_position_corr_significant=elec_response_to_word_position_corr_p_val<p_val_threshold;
    % choose  color based on the correlation between open node and 
    [rho,p_val]=cellfun(@(x,y) corr(x,y,'type','Spearman'),unique_sentence_pattern_cell,word_position_cell,'UniformOutput',false);
    corr_range=([rho{:}]-min([rho{:}]))./(max([rho{:}])-min([rho{:}]));
    color_range=transpose(([1,1,1]')*corr_range);
    %
    significance_overlap= 1 |elec_response_to_open_node_corr_significant | elec_response_to_word_position_corr_significant;
    positive_rhos= elec_response_to_word_position_corr>0 & elec_response_to_open_node_corr>0;
    elec_response_to_open_node_corr=elec_response_to_open_node_corr(positive_rhos);
    elec_response_to_word_position_corr=elec_response_to_word_position_corr(positive_rhos);
    color_range=color_range(positive_rhos,:);
  % do the do the correlation;
    figure(fix((i-1)/total_plots)+1);
    set(gcf,'position',[1451,-285,957,1342]);
    sub_title=sprintf('subj: %s,\n electrodes: %d ', info.subject ,language_electrode_num(i));
    ax=subplot(num_rows,num_columns,i-total_plots*fix((i-1)/total_plots));
    figures=scatter_corner_hist(elec_response_to_word_position_corr,elec_response_to_open_node_corr,'word corr','open node corr',color_range);
    sub_title=sprintf('subj: %s,\n electrodes: %d ', info.subject ,language_electrode_num(i));
    figures(1).Title=title(sub_title);
    figures(2).FontSize=10;
    if ~mod(i,total_plots) | i==size(session_sentence_hilbert_band_envelope_lang_elec_tensor,1)
    if ~exist(strcat(analysis_path,info.subject))
        mkdir(strcat(analysis_path,info.subject))
    end
    print(gcf, '-depsc', strcat(analysis_path,info.subject,'/',info.subject,'_ch_',num2str(language_electrode_num(i)),'_corr_unique_open_node_vs_word_pos_mean.eps')); 
    close(gcf)
    end   
end

% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  compute correlation based on unique patterns by taking the repetition of values   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all 
unique_sentence_pattern=mat2cell(unique(all_sentence_pattern,'row'),ones(1,size(unique(all_sentence_pattern,'row'),1)),size(all_sentence_pattern,2));
uniq_sent_patt_loc_in_all_sentences=cellfun(@(x) find(ismember(all_sentence_pattern,x,'row')),unique_sentence_pattern,'UniformOutput',false);
num_rows=4;
num_columns=2;
total_plots=num_rows*num_columns;
look_out_window=100; %200 ms
offset=30; % 100 ms
p_val_threshold=0.05;
word_range=[3:8];
% 
for i=1:length(language_electrode_num)
    elec_response_to_open_node_pattern=[];
    electrode_post_word_sentence_resp_cell=[];
    elec_response_to_open_node_corr=[];
    elec_response_to_word_position_corr=[];
    unique_sentence_pattern_cell=[];
    word_position_cell=[];
    electrode_response=double(squeeze(session_sentence_hilbert_band_envelope_lang_elec_tensor(i,:,:)));
    electrode_response_sentence_cell=mat2cell(electrode_response,ones(1,size(electrode_response,1)),word_length*ones(1,8));
    electrode_post_word_sentence_mean=cell2mat(cellfun(@(x) nanmean(x(offset+(1:look_out_window))),electrode_response_sentence_cell,'UniformOutput',false));
    for k=1:size(unique_sentence_pattern,1)
        pattern_responses=transpose(electrode_post_word_sentence_mean(uniq_sent_patt_loc_in_all_sentences{k},:));
        pattern_responses=pattern_responses(word_range,:);
        %
        sentence_open_node_pattern=repmat(transpose(unique_sentence_pattern{k}),1,length(uniq_sent_patt_loc_in_all_sentences{k}));
        sentence_open_node_pattern=sentence_open_node_pattern(word_range,:);
        %
        sentence_word_position_pattern=transpose(repmat(1:size(transpose(electrode_post_word_sentence_mean(uniq_sent_patt_loc_in_all_sentences{k},:)),1),...
            size(transpose(electrode_post_word_sentence_mean(uniq_sent_patt_loc_in_all_sentences{k},:)),2),1));
        sentence_word_position_pattern=sentence_word_position_pattern(word_range,:);
        %
        electrode_post_word_sentence_resp_cell{k,1}=transpose(reshape(pattern_responses,1,[]));
        unique_sentence_pattern_cell{k,1}=transpose(reshape(sentence_open_node_pattern,1,[]));
        word_position_cell{k,1}=transpose(reshape(sentence_word_position_pattern,1,[]));
    end
    % compute correlation for open node
    [rho,p_val]=cellfun(@(x,y) corr(x,y,'type','Spearman'),electrode_post_word_sentence_resp_cell,unique_sentence_pattern_cell,'UniformOutput',false);
    elec_response_to_open_node_corr=[rho{:}];
    elec_response_to_open_node_corr_p_val=[p_val{:}];
    elec_response_to_open_node_corr_significant=elec_response_to_open_node_corr_p_val<p_val_threshold;
    % compute correlation for wregord position 
    [rho,p_val]=cellfun(@(x,y) corr(x,y,'type','Spearman'),electrode_post_word_sentence_resp_cell,word_position_cell,'UniformOutput',false);
    elec_response_to_word_position_corr=[rho{:}];
    elec_response_to_word_position_corr_p_val=[p_val{:}];
    elec_response_to_word_position_corr_significant=elec_response_to_word_position_corr_p_val<p_val_threshold;
    % choose  color based on the correlation between open node and 
    [rho,p_val]=cellfun(@(x,y) corr(x,y,'type','Spearman'),unique_sentence_pattern_cell,word_position_cell,'UniformOutput',false);
    corr_range=([rho{:}]-min([rho{:}]))./(max([rho{:}])-min([rho{:}]));
    color_range=transpose(([1,1,1]')*corr_range);
    %
    significance_overlap= 1 |elec_response_to_open_node_corr_significant | elec_response_to_word_position_corr_significant;
    positive_rhos= elec_response_to_word_position_corr>0 & elec_response_to_open_node_corr>0;
    elec_response_to_open_node_corr=elec_response_to_open_node_corr(positive_rhos);
    elec_response_to_word_position_corr=elec_response_to_word_position_corr(positive_rhos);
    color_range=color_range(positive_rhos,:);
  % do the do the correlation;
    figure(fix((i-1)/total_plots)+1);
    set(gcf,'position',[1451,-285,957,1342]);
    sub_title=sprintf('subj: %s,\n electrodes: %d ', info.subject ,language_electrode_num(i));
    ax=subplot(num_rows,num_columns,i-total_plots*fix((i-1)/total_plots));
    figures=scatter_corner_hist(elec_response_to_word_position_corr,elec_response_to_open_node_corr,'word corr','open node corr',color_range);
    sub_title=sprintf('subj: %s,\n electrodes: %d ', info.subject ,language_electrode_num(i));
    figures(1).Title=title(sub_title);
    figures(2).FontSize=10;

    if ~mod(i,total_plots) | i==size(session_sentence_hilbert_band_envelope_lang_elec_tensor,1)
    if ~exist(strcat(analysis_path,info.subject))
        mkdir(strcat(analysis_path,info.subject))
    end
    print(gcf, '-depsc', strcat(analysis_path,info.subject,'/',info.subject,'_ch_',num2str(language_electrode_num(i)),'_corr_unique_open_node_vs_word_pos_repetition.eps')); 
    close(gcf)
    end   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  compute regression unique patterns by taking the repetition of values   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all 
unique_sentence_pattern=mat2cell(unique(all_sentence_pattern,'row'),ones(1,size(unique(all_sentence_pattern,'row'),1)),size(all_sentence_pattern,2));
uniq_sent_patt_loc_in_all_sentences=cellfun(@(x) find(ismember(all_sentence_pattern,x,'row')),unique_sentence_pattern,'UniformOutput',false);
num_rows=4;
num_columns=2;
total_plots=num_rows*num_columns;
look_out_window=100; %200 ms
offset=30; % 100 ms
p_val_threshold=0.05;
word_range=[3:8];
% 
for i=1:length(language_electrode_num)
    elec_response_to_open_node_pattern=[];
    electrode_post_word_sentence_resp_cell=[];
    elec_response_to_open_node_corr=[];
    elec_response_to_word_position_corr=[];
    unique_sentence_pattern_cell=[];
    word_position_cell=[];
    electrode_response=double(squeeze(session_sentence_hilbert_band_envelope_lang_elec_tensor(i,:,:)));
    electrode_response_sentence_cell=mat2cell(electrode_response,ones(1,size(electrode_response,1)),word_length*ones(1,8));
    electrode_post_word_sentence_mean=cell2mat(cellfun(@(x) nanmean(x(offset+(1:look_out_window))),electrode_response_sentence_cell,'UniformOutput',false));
    for k=1:size(unique_sentence_pattern,1)
        pattern_responses=transpose(electrode_post_word_sentence_mean(uniq_sent_patt_loc_in_all_sentences{k},:));
        pattern_responses=pattern_responses(word_range,:);
        %
        sentence_open_node_pattern=repmat(transpose(unique_sentence_pattern{k}),1,length(uniq_sent_patt_loc_in_all_sentences{k}));
        sentence_open_node_pattern=sentence_open_node_pattern(word_range,:);
        %
        sentence_word_position_pattern=transpose(repmat(1:size(transpose(electrode_post_word_sentence_mean(uniq_sent_patt_loc_in_all_sentences{k},:)),1),...
            size(transpose(electrode_post_word_sentence_mean(uniq_sent_patt_loc_in_all_sentences{k},:)),2),1));
        sentence_word_position_pattern=sentence_word_position_pattern(word_range,:);
        %
        electrode_post_word_sentence_resp_cell{k,1}=transpose(reshape(pattern_responses,1,[]));
        unique_sentence_pattern_cell{k,1}=transpose(reshape(sentence_open_node_pattern,1,[]));
        word_position_cell{k,1}=transpose(reshape(sentence_word_position_pattern,1,[]));
    end
    % compute correlation for open node
    [beta]=cellfun(@(x,y,z) ridge(x,[y,z],0.3,0),electrode_post_word_sentence_resp_cell,unique_sentence_pattern_cell,word_position_cell,'UniformOutput',false);
    electrode_beta_open_node_pattern=cell2mat(cellfun(@(x) x(2),beta,'UniformOutput',false));
    electrode_beta_word_pos_pattern=cell2mat(cellfun(@(x) x(3),beta,'UniformOutput',false));
    % 
    [rho,p_val]=cellfun(@(x,y) corr(x,y,'type','Spearman'),unique_sentence_pattern_cell,word_position_cell,'UniformOutput',false);
    corr_range=([rho{:}]-min([rho{:}]))./(max([rho{:}])-min([rho{:}]));
    color_range=transpose(([1,1,1]')*corr_range);
    %
  % do the do the correlation;
    figure(fix((i-1)/total_plots)+1);
    set(gcf,'position',[1451,-285,957,1342]);
    sub_title=sprintf('subj: %s,\n electrodes: %d ', info.subject ,language_electrode_num(i));
    ax=subplot(num_rows,num_columns,i-total_plots*fix((i-1)/total_plots));
    figures=scatter_corner_hist(electrode_beta_word_pos_pattern,electrode_beta_open_node_pattern,'word beta','open node beta',color_range);
    sub_title=sprintf('subj: %s,\n electrodes: %d ', info.subject ,language_electrode_num(i));
    figures(1).Title=title(sub_title);
    figures(2).FontSize=10;

    if ~mod(i,total_plots) | i==size(session_sentence_hilbert_band_envelope_lang_elec_tensor,1)
    if ~exist(strcat(analysis_path,info.subject))
        mkdir(strcat(analysis_path,info.subject))
    end
    print(gcf, '-depsc', strcat(analysis_path,info.subject,'/',info.subject,'_ch_',num2str(language_electrode_num(i)),'_corr_unique_open_node_vs_word_pos_repetition.eps')); 
    close(gcf)
    end   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  compute regression independently for unique patterns 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all 
unique_sentence_pattern=mat2cell(unique(all_sentence_pattern,'row'),ones(1,size(unique(all_sentence_pattern,'row'),1)),size(all_sentence_pattern,2));
uniq_sent_patt_loc_in_all_sentences=cellfun(@(x) find(ismember(all_sentence_pattern,x,'row')),unique_sentence_pattern,'UniformOutput',false);
num_rows=4;
num_columns=2;
total_plots=num_rows*num_columns;
look_out_window=100; %200 ms
offset=30; % 100 ms
p_val_threshold=0.05;
word_range=[3:8];
% 
for i=1:length(language_electrode_num)
    elec_response_to_open_node_pattern=[];
    electrode_post_word_sentence_resp_cell=[];
    elec_response_to_open_node_corr=[];
    elec_response_to_word_position_corr=[];
    unique_sentence_pattern_cell=[];
    word_position_cell=[];
    electrode_response=double(squeeze(session_sentence_hilbert_band_envelope_lang_elec_tensor(i,:,:)));
    electrode_response_sentence_cell=mat2cell(electrode_response,ones(1,size(electrode_response,1)),word_length*ones(1,8));
    electrode_post_word_sentence_mean=cell2mat(cellfun(@(x) nanmean(x(offset+(1:look_out_window))),electrode_response_sentence_cell,'UniformOutput',false));
    for k=1:size(unique_sentence_pattern,1)
        pattern_responses=transpose(electrode_post_word_sentence_mean(uniq_sent_patt_loc_in_all_sentences{k},:));
        pattern_responses=pattern_responses(word_range,:);
        %
        sentence_open_node_pattern=repmat(transpose(unique_sentence_pattern{k}),1,length(uniq_sent_patt_loc_in_all_sentences{k}));
        sentence_open_node_pattern=sentence_open_node_pattern(word_range,:);
        %
        sentence_word_position_pattern=transpose(repmat(1:size(transpose(electrode_post_word_sentence_mean(uniq_sent_patt_loc_in_all_sentences{k},:)),1),...
            size(transpose(electrode_post_word_sentence_mean(uniq_sent_patt_loc_in_all_sentences{k},:)),2),1));
        sentence_word_position_pattern=sentence_word_position_pattern(word_range,:);
        %
        electrode_post_word_sentence_resp_cell{k,1}=transpose(reshape(pattern_responses,1,[]));
        unique_sentence_pattern_cell{k,1}=transpose(reshape(sentence_open_node_pattern,1,[]));
        word_position_cell{k,1}=transpose(reshape(sentence_word_position_pattern,1,[]));
    end
    
    % compute correlation for open node
    [beta_open_pos,bint_open_pos,r_open_pos,rint_open_pos,stat_open_pos]=cellfun(@(x,y) regress(y,[ones(size(x)),x]),unique_sentence_pattern_cell,electrode_post_word_sentence_resp_cell,'UniformOutput',false);
    [beta_word_pos,bint_word_pos,r_word_pos,rint_word_pos,stat_word_pos]=cellfun(@(x,y) regress(y,[ones(size(x)),x]),word_position_cell,electrode_post_word_sentence_resp_cell,'UniformOutput',false);
    electrode_r_squared_open_node_pattern=cell2mat(cellfun(@(x) x(1),stat_open_pos,'UniformOutput',false));
    electrode_r_squared_word_pos_pattern=cell2mat(cellfun(@(x) x(1),stat_word_pos,'UniformOutput',false));
    % 
    [rho,p_val]=cellfun(@(x,y) corr(x,y,'type','Spearman'),unique_sentence_pattern_cell,word_position_cell,'UniformOutput',false);
    corr_range=([rho{:}]-min([rho{:}]))./(max([rho{:}])-min([rho{:}]));
    color_range=transpose(([1,1,1]')*corr_range);
    %
  % do the do the correlation;
    figure(fix((i-1)/total_plots)+1);
    set(gcf,'position',[1442 -160 957 1342]);
    sub_title=sprintf('subj: %s,\n electrodes: %d ', info.subject ,language_electrode_num(i));
    ax=subplot(num_rows,num_columns,i-total_plots*fix((i-1)/total_plots));
    figures=scatter_corner_hist(electrode_r_squared_word_pos_pattern,electrode_r_squared_open_node_pattern,'word pos r^2','open node r^2',color_range);
    sub_title=sprintf('subj: %s,\n electrodes: %d ', info.subject ,language_electrode_num(i));
    figures(1).Title=title(sub_title);
    figures(2).FontSize=10;

    if ~mod(i,total_plots) | i==size(session_sentence_hilbert_band_envelope_lang_elec_tensor,1)
    if ~exist(strcat(analysis_path,info.subject))
        mkdir(strcat(analysis_path,info.subject))
    end
%    print(gcf, '-depsc', strcat(analysis_path,info.subject,'/',info.subject,'_ch_',num2str(language_electrode_num(i)),'_corr_unique_open_node_vs_word_pos_repetition.eps')); 
   % close(gcf)
    end   
end



