%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 0: prepare the workspace 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all 
close all 
data_path='/Users/eghbalhosseiniasl1/MyData/ecog-sentence';
subject_id='AMC038';
d= dir([data_path,strcat('/**/',subject_id,'*_crunched.mat')]);
fprintf(' %d .mat files were found \n', length(d));
gamma_band_index=4;
save_path='/Users/eghbalhosseiniasl1/MyData/ECoG-sentence/crunched/';
analysis_path='/Users/eghbalhosseiniasl1/MyData/ECoG-sentence/analysis/swjn_analysis_3_relate_gamma_to_structural_probe_BERT/';
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
all_example_sentences=[];
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
    all_example_sentences=[all_example_sentences;example_sentence'];
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
%%  STEP 2: get the structural probe data    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
probe_path='/Users/eghbalhosseiniasl1/MyData/ecog-sentence/structural_probe/';
BERT_representation=load(strcat(probe_path,'representation_cells.mat'));
structural_distance_prediction=load(strcat(probe_path,'distance_predictions_cells.mat'));
structural_depth_prediction=load(strcat(probe_path,'depth_predictions_cells.mat'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find probe data 
structure_sentences=BERT_representation.representation_arr(:,1);
example_BERT_representation=[];
example_structural_distance_prediction=[];
example_structural_depth_prediction=[];
for i=1:length(all_example_sentences)
    indx=cell2mat(cellfun(@(x) contains(x,all_example_sentences{i},'IgnoreCase',true),structure_sentences,'UniformOutput',false));
    example_BERT_representation=[example_BERT_representation;BERT_representation.representation_arr(find(indx),2)];
    example_structural_distance_prediction=[example_structural_distance_prediction;...
        structural_distance_prediction.distance_predictions_arr(find(indx),2)];
    example_structural_depth_prediction=[example_structural_depth_prediction;...
        structural_depth_prediction.depth_predictions_arr(find(indx),2)];    
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
word_range=[1:8];
j=0;
[row,column]=ind2sub(size(ones(2,3)),find(ones(2,3)));
for i=1:size(session_sentence_hilbert_band_envelope_lang_elec_tensor,1)
    electrode_response=squeeze(session_sentence_hilbert_band_envelope_lang_elec_tensor(i,:,:));      
    electrode_response_sentence_cell=mat2cell(electrode_response,ones(1,size(electrode_response,1)),word_length*ones(1,8));
    electrode_post_word_sentence_mean=double(cell2mat(cellfun(@(x) nanmean(x(offset:end)),electrode_response_sentence_cell,'UniformOutput',false)));
    word_position=repmat(1:size(electrode_post_word_sentence_mean,2),[size(electrode_post_word_sentence_mean,1),1]);
    % 
    electrode_post_word_sentence_mean_cell=cellfun(@transpose,mat2cell(electrode_post_word_sentence_mean,ones(1,size(electrode_post_word_sentence_mean,1)),size(electrode_post_word_sentence_mean,2)),'UniformOutput',false);
    all_sentence_pattern_cell=cellfun(@transpose,mat2cell(all_sentence_pattern,ones(1,size(all_sentence_pattern,1)),size(all_sentence_pattern,2)),'uniformoutput',false);
    word_position_cell=cellfun(@transpose,mat2cell(word_position,ones(1,size(word_position,1)),size(word_position,2)),'uniformoutput',false);
    
    % select a subset of words 
    electrode_post_word_sentence_mean_cell=cellfun(@(x) x(word_range),electrode_post_word_sentence_mean_cell,'UniformOutput',false);
    all_structural_depth_cell=cellfun(@(x) double(transpose(x(word_range))),example_structural_depth_prediction,'UniformOutput',false);
    word_position_cell=cellfun(@(x) (x(word_range)),word_position_cell,'UniformOutput',false);
    
    [rho,p_val]=cellfun(@(x,y) corr(x,y,'type','Spearman'),electrode_post_word_sentence_mean_cell,all_structural_depth_cell,'UniformOutput',false);
    elec_response_to_structural_depth_corr=[rho{:}];
    elec_response_to_structural_depth_corr_p_val=[p_val{:}];
    elec_response_to_structural_depth_corr_significant=elec_response_to_structural_depth_corr_p_val<p_val_threshold;
    
    % 
    [rho,p_val]=cellfun(@(x,y) corr(x,y,'type','Spearman'),electrode_post_word_sentence_mean_cell,word_position_cell,'UniformOutput',false);
    elec_response_to_word_position_corr=[rho{:}];
    elec_response_to_word_position_corr_p_val=[p_val{:}];
    elec_response_to_word_position_corr_significant=elec_response_to_word_position_corr_p_val<p_val_threshold;
    % 
    positive_rhos= elec_response_to_word_position_corr>0 & elec_response_to_structural_depth_corr>0;
    elec_response_to_structural_depth_corr=elec_response_to_structural_depth_corr(positive_rhos);
    elec_response_to_word_position_corr=elec_response_to_word_position_corr(positive_rhos);
    if ~ishandle((fix((i-1)/total_plots)+1))
    figure(fix((i-1)/total_plots)+1);
    set(gcf,'position',[1451,-285,957,1342]);
    end 
    
    sub_title=sprintf('subj: %s,\n electrodes: %d ', info.subject ,language_electrode_num(i));
    ax=subplot(num_rows,num_columns,i-total_plots*fix((i-1)/total_plots));
    figures=scatter_corner_hist(elec_response_to_word_position_corr,elec_response_to_structural_depth_corr,'word corr','strucrual depth corr',[1,0,0]);
    sub_title=sprintf('subj: %s,\n electrodes: %d ', info.subject ,language_electrode_num(i));
    figures(2).FontSize=10;
    figures(1).Title=title(sub_title);
    if ~mod(i,total_plots) | i==size(session_sentence_hilbert_band_envelope_lang_elec_tensor,1)
    if ~exist(strcat(analysis_path,info.subject))
        mkdir(strcat(analysis_path,info.subject))
    end
    set(gcf,'PaperPosition',[.25 .25 8 6])
    pause(1)
    print(gcf, '-fillpage','-dpdf', strcat(analysis_path,info.subject,'/',info.subject,'_ch_',num2str(language_electrode_num(i)),'_corr_structural_depth_vs_word_pos.pdf')); 
    close(gcf)
    end 
end 



