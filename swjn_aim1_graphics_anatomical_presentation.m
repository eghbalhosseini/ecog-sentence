clear all 
close all 
analysis_path='/Users/eghbalhosseiniasl1/MyData/ECoG-sentence/analysis/swjn_Aim_1_graphics';
data_path='/Users/eghbalhosseiniasl1/MyData/ecog-sentence';
if 1
    fprintf('adding basic ecog tools to path \n');
    addpath('~/MyCodes/basic-ecog-tools/');
    addpath(genpath('~/MyCodes/ecog-sentence/'));
    addpath(genpath('~/MyCodes/basic-ecog-tools/activeBrain'));
    addpath(genpath('~/MyCodes/basic-ecog-tools/ecog-filters'));
    addpath(genpath('~/MyCodes/basic-ecog-tools/mex'));
end
d_pmi= dir([data_path,'/**/*pmi_sentences.csv']);
fprintf(' %d pmi files were found \n', length(d_pmi));
i=1;
[pmi_cell]=generate_pmi_table_eh(strcat(d_pmi(i).folder,'/',d_pmi(i).name));
subject_id='AMC026';
d= dir([data_path,strcat('/**/',subject_id,'*_crunched.mat')]);
fprintf(' %d .mat files were found \n', length(d));

d_nmerge= dir([data_path,'/**/*nmerge.txt']);
fprintf(' %d nmerge files were found \n', length(d_nmerge));
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

%%  
AMC026_elec_fsaverage=load(['/Users/eghbalhosseiniasl1/MyData/ecog-sentence/subjects_raw/AMC026/IMAGING/','AMC026_elec_fsavg.mat']);
AMC026_elec_fsaverage=AMC026_elec_fsaverage.elec_fsavg_frs;
AMC029_elec_fsaverage=load(['/Users/eghbalhosseiniasl1/MyData/ecog-sentence/subjects_raw/AMC029/IMAGING/','AMC029_elec_fsavg.mat']);
AMC029_elec_fsaverage=AMC029_elec_fsaverage.elec_fsavg_frs;

%% get 

find_index= @(input_cell) find(cell2mat(cellfun(@(x) ~isempty(x), input_cell,'UniformOutput',false)));
sentence_patterns=[];
wordlist_patterns=[];
nonwords_patterns=[];
jabber_patterns=[];

sentence_patterns={};
wordlist_patterns={};

session_sentence_hilbert_band_envelope_lang_elec_tensor=[];
session_sentence_hilbert_zs_band_envelope_lang_elec_tensor=[];
session_nonwords_hilbert_band_envelope_lang_elec_tensor=[];
session_words_hilbert_band_envelope_lang_elec_tensor=[];
session_jabberwocky_hilbert_band_envelope_lang_elec_tensor=[];

session_sentence_hilbert_band_envelope_pre_trial_tensor=[];
session_sentence_hilbert_zs_band_envelope_pre_trial_tensor=[];
session_sentence_hilbert_band_ave_envelope_pre_trial_tensor=[];
session_sentence_hilbert_band_envelope_tensor=[];
session_sentence_hilbert_zs_band_envelope_tensor=[];
session_sentence_hilbert_band_ave_envelope_tensor=[];

session_sentence_hilbert_pca_zs_band_envelope_tensor=[];
session_sentence_hilbert_pca_zs_band_envelope_pre_trial_tensor=[];
% 
session_words_hilbert_band_envelope_pre_trial_tensor=[];
session_words_hilbert_zs_band_envelope_pre_trial_tensor=[];
session_words_hilbert_pca_zs_band_envelope_pre_trial_tensor=[];
session_words_hilbert_band_ave_envelope_pre_trial_tensor=[];
session_words_hilbert_band_envelope_tensor=[];
session_words_hilbert_zs_band_envelope_tensor=[];
session_words_hilbert_pca_zs_band_envelope_tensor=[];
session_words_hilbert_band_ave_envelope_tensor=[];
% 
session_nonwords_hilbert_band_envelope_pre_trial_tensor=[];
session_nonwords_hilbert_zs_band_envelope_pre_trial_tensor=[];
session_nonwords_hilbert_band_ave_envelope_pre_trial_tensor=[];
session_nonwords_hilbert_band_envelope_tensor=[];
session_nonwords_hilbert_zs_band_envelope_tensor=[];
session_nonwords_hilbert_band_ave_envelope_tensor=[];
session_nonwords_hilbert_pca_zs_band_envelope_tensor=[];
session_nonwords_hilbert_pca_zs_band_envelope_pre_trial_tensor=[];

% 
session_jabberwocky_hilbert_band_envelope_pre_trial_tensor=[];
session_jabberwocky_hilbert_zs_band_envelope_pre_trial_tensor=[];
session_jabberwocky_hilbert_band_ave_envelope_pre_trial_tensor=[];
session_jabberwocky_hilbert_band_envelope_tensor=[];
session_jabberwocky_hilbert_zs_band_envelope_tensor=[];
session_jabberwocky_hilbert_band_ave_envelope_tensor=[];
session_jabberwocky_hilbert_pca_zs_band_envelope_tensor=[];
session_jabber_hilbert_pca_zs_band_envelope_pre_trial_tensor=[];



for i=1:length(d)
    subj=load(strcat(d(i).folder,'/',d(i).name));
    subj_id=fieldnames(subj);
    subj=subj.(subj_id{1});
    data=subj.data;
    info=subj.info;
    language_electrode=info.ramp_electrodes_hilbert_odd;
    language_electrode_num=find(language_electrode);
    
    % step 1: extract electrodes with siginificant language response
    sentence_trial_index=~cellfun(@isempty,cellfun(@ (x) strfind(x,'S'),info.word_type,'UniformOutput',false));
    sentences=[data{sentence_trial_index}];
    %
    hilbert_band_ave_envelope=cellfun(@(x) x(1:8),{sentences.signal_ave_hilbert_downsample_parsed},'UniformOutput',false);
    hilbert_band_ave_envelope=[hilbert_band_ave_envelope{:,:}];
    hilbert_band_ave_envelope_pre_trial=cellfun(@(x) reshape(x,size(x,1),1,size(x,2)),{sentences.signal_ave_pre_trial_hilbert_downsample},'UniformOutput',false);
    hilbert_band_ave_envelope_pre_trial=[hilbert_band_ave_envelope_pre_trial{:,:}];
    hilbert_band_envelope=cellfun(@(x) x(1:8),{sentences.signal_hilbert_downsample_parsed},'UniformOutput',false);
    hilbert_band_envelope=[hilbert_band_envelope{:,:}];
    hilbert_band_envelope_pre_trial=cellfun(@(x) reshape(x,size(x,1),1,size(x,2)),{sentences.signal_pre_trial_hilbert_downsample},'UniformOutput',false);
    hilbert_band_envelope_pre_trial=[hilbert_band_envelope_pre_trial{:,:}];
    hilbert_zs_band_envelope=cellfun(@(x) x(1:8),{sentences.signal_hilbert_zs_downsample_parsed},'UniformOutput',false);
    hilbert_zs_band_envelope=[hilbert_zs_band_envelope{:,:}];
    hilbert_zs_band_envelope_pre_trial=cellfun(@(x) reshape(x,size(x,1),1,size(x,2)),{sentences.signal_pre_trial_hilbert_zs_downsample},'UniformOutput',false);
    hilbert_zs_band_envelope_pre_trial=[hilbert_zs_band_envelope_pre_trial{:,:}];
    try 
        hilbert_pca_zs_band_envelope=cellfun(@(x) x(1:8),{sentences.signal_hilbert_pca_zs_downsample_parsed},'UniformOutput',false);
        hilbert_pca_zs_band_envelope=[hilbert_pca_zs_band_envelope{:,:}];
        hilbert_pca_zs_band_envelope_pre_trial=cellfun(@(x) reshape(x,size(x,1),1,size(x,2)),{sentences.signal_pre_trial_hilbert_pca_zs_downsample},'UniformOutput',false);
        hilbert_pca_zs_band_envelope_pre_trial=[hilbert_pca_zs_band_envelope_pre_trial{:,:}];
    end 
    %
    sentence_hilbert_band_ave_envelope_pre_trial_tensor=permute(hilbert_band_ave_envelope_pre_trial,[1,3,2]);
    sentence_hilbert_band_envelope_pre_trial_tensor=permute(hilbert_band_envelope_pre_trial,[1,3,2]);
    sentence_hilbert_band_ave_envelope_tensor=cell2mat(permute(hilbert_band_ave_envelope,[3,1,2]));
    sentence_hilbert_band_envelope_tensor=cell2mat(permute(hilbert_band_envelope,[3,1,2]));
    sentence_hilbert_zs_band_envelope_pre_trial_tensor=permute(hilbert_zs_band_envelope_pre_trial,[1,3,2]);
    sentence_hilbert_zs_band_envelope_tensor=cell2mat(permute(hilbert_zs_band_envelope,[3,1,2]));
    try 
         sentence_hilbert_pca_zs_band_envelope_pre_trial_tensor=permute(hilbert_pca_zs_band_envelope_pre_trial,[1,3,2]);
    sentence_hilbert_pca_zs_band_envelope_tensor=cell2mat(permute(hilbert_pca_zs_band_envelope,[3,1,2]));
    end 
    %
    session_sentence_hilbert_band_envelope_pre_trial_tensor=cat(3,session_sentence_hilbert_band_envelope_pre_trial_tensor,sentence_hilbert_band_envelope_pre_trial_tensor);
    session_sentence_hilbert_band_ave_envelope_pre_trial_tensor=cat(3,session_sentence_hilbert_band_ave_envelope_pre_trial_tensor,sentence_hilbert_band_ave_envelope_pre_trial_tensor);
    session_sentence_hilbert_band_envelope_tensor=cat(3,session_sentence_hilbert_band_envelope_tensor,sentence_hilbert_band_envelope_tensor);
    session_sentence_hilbert_band_ave_envelope_tensor=cat(3,session_sentence_hilbert_band_ave_envelope_tensor,sentence_hilbert_band_ave_envelope_tensor);
    session_sentence_hilbert_zs_band_envelope_pre_trial_tensor=cat(3,session_sentence_hilbert_zs_band_envelope_pre_trial_tensor,sentence_hilbert_zs_band_envelope_pre_trial_tensor);
    session_sentence_hilbert_zs_band_envelope_tensor=cat(3,session_sentence_hilbert_zs_band_envelope_tensor,sentence_hilbert_zs_band_envelope_tensor);
    try 
        session_sentence_hilbert_pca_zs_band_envelope_pre_trial_tensor=cat(3,session_sentence_hilbert_pca_zs_band_envelope_pre_trial_tensor,sentence_hilbert_pca_zs_band_envelope_pre_trial_tensor);
    session_sentence_hilbert_pca_zs_band_envelope_tensor=cat(3,session_sentence_hilbert_pca_zs_band_envelope_tensor,sentence_hilbert_pca_zs_band_envelope_tensor);
    end
    
    % find sentence locations
    example_sentence=cellfun(@(x) x(2:end),{sentences(:).trial_string},'UniformOutput',false);
    sentence_patterns=[sentence_patterns,example_sentence];
    
    %%%%%%%%%%%%%%%%% words
    %%%%%%%%%%%%%%%%% words
    words_trial_index=~cellfun(@isempty,cellfun(@ (x) strfind(x,'W'),info.word_type,'UniformOutput',false));
    words=[data{words_trial_index}];
    %
    hilbert_band_ave_envelope=cellfun(@(x) x(1:8),{words.signal_ave_hilbert_downsample_parsed},'UniformOutput',false);
    hilbert_band_ave_envelope=[hilbert_band_ave_envelope{:,:}];
    hilbert_band_ave_envelope_pre_trial=cellfun(@(x) reshape(x,size(x,1),1,size(x,2)),{words.signal_ave_pre_trial_hilbert_downsample},'UniformOutput',false);
    hilbert_band_ave_envelope_pre_trial=[hilbert_band_ave_envelope_pre_trial{:,:}];
    hilbert_band_envelope=cellfun(@(x) x(1:8),{words.signal_hilbert_downsample_parsed},'UniformOutput',false);
    hilbert_band_envelope=[hilbert_band_envelope{:,:}];
    hilbert_band_envelope_pre_trial=cellfun(@(x) reshape(x,size(x,1),1,size(x,2)),{words.signal_pre_trial_hilbert_downsample},'UniformOutput',false);
    hilbert_band_envelope_pre_trial=[hilbert_band_envelope_pre_trial{:,:}];
    hilbert_zs_band_envelope=cellfun(@(x) x(1:8),{words.signal_hilbert_zs_downsample_parsed},'UniformOutput',false);
    hilbert_zs_band_envelope=[hilbert_zs_band_envelope{:,:}];
    hilbert_zs_band_envelope_pre_trial=cellfun(@(x) reshape(x,size(x,1),1,size(x,2)),{words.signal_pre_trial_hilbert_zs_downsample},'UniformOutput',false);
    hilbert_zs_band_envelope_pre_trial=[hilbert_zs_band_envelope_pre_trial{:,:}];
    try
        hilbert_pca_zs_band_envelope=cellfun(@(x) x(1:8),{words.signal_hilbert_pca_zs_downsample_parsed},'UniformOutput',false);
        hilbert_pca_zs_band_envelope=[hilbert_pca_zs_band_envelope{:,:}];
        hilbert_pca_zs_band_envelope_pre_trial=cellfun(@(x) reshape(x,size(x,1),1,size(x,2)),{words.signal_pre_trial_hilbert_pca_zs_downsample},'UniformOutput',false);
        hilbert_pca_zs_band_envelope_pre_trial=[hilbert_pca_zs_band_envelope_pre_trial{:,:}];
    end
 
    %
    words_hilbert_band_ave_envelope_pre_trial_tensor=permute(hilbert_band_ave_envelope_pre_trial,[1,3,2]);
    words_hilbert_band_envelope_pre_trial_tensor=permute(hilbert_band_envelope_pre_trial,[1,3,2]);
    words_hilbert_band_ave_envelope_tensor=cell2mat(permute(hilbert_band_ave_envelope,[3,1,2]));
    words_hilbert_band_envelope_tensor=cell2mat(permute(hilbert_band_envelope,[3,1,2]));
    words_hilbert_zs_band_envelope_pre_trial_tensor=permute(hilbert_zs_band_envelope_pre_trial,[1,3,2]);
    words_hilbert_zs_band_envelope_tensor=cell2mat(permute(hilbert_zs_band_envelope,[3,1,2]));
    try
        sentence_hilbert_pca_zs_band_envelope_pre_trial_tensor=permute(hilbert_pca_zs_band_envelope_pre_trial,[1,3,2]);
        sentence_hilbert_pca_zs_band_envelope_tensor=cell2mat(permute(hilbert_pca_zs_band_envelope,[3,1,2]));
    end
    
    %
    session_words_hilbert_band_envelope_pre_trial_tensor=cat(3,session_words_hilbert_band_envelope_pre_trial_tensor,words_hilbert_band_envelope_pre_trial_tensor);
    session_words_hilbert_band_ave_envelope_pre_trial_tensor=cat(3,session_words_hilbert_band_ave_envelope_pre_trial_tensor,words_hilbert_band_ave_envelope_pre_trial_tensor);
    session_words_hilbert_band_envelope_tensor=cat(3,session_words_hilbert_band_envelope_tensor,words_hilbert_band_envelope_tensor);
    session_words_hilbert_band_ave_envelope_tensor=cat(3,session_words_hilbert_band_ave_envelope_tensor,words_hilbert_band_ave_envelope_tensor);
    session_words_hilbert_zs_band_envelope_pre_trial_tensor=cat(3,session_words_hilbert_zs_band_envelope_pre_trial_tensor,words_hilbert_zs_band_envelope_pre_trial_tensor);
    session_words_hilbert_zs_band_envelope_tensor=cat(3,session_words_hilbert_zs_band_envelope_tensor,words_hilbert_zs_band_envelope_tensor);
        try 
        session_words_hilbert_pca_zs_band_envelope_pre_trial_tensor=cat(3,session_words_hilbert_pca_zs_band_envelope_pre_trial_tensor,sentence_hilbert_pca_zs_band_envelope_pre_trial_tensor);
    session_words_hilbert_pca_zs_band_envelope_tensor=cat(3,session_words_hilbert_pca_zs_band_envelope_tensor,sentence_hilbert_pca_zs_band_envelope_tensor);
        end
     %
    example_words=cellfun(@(x) x(2:end),{words(:).trial_string},'UniformOutput',false);
    wordlist_patterns=[wordlist_patterns,example_words];
    % 
    %%%%%%%%%%%%%%%% nonwords
    nonwords_trial_index=~cellfun(@isempty,cellfun(@ (x) strfind(x,'N'),info.word_type,'UniformOutput',false));
    nonwords=[data{nonwords_trial_index}];
    %
    hilbert_band_ave_envelope=cellfun(@(x) x(1:8),{nonwords.signal_ave_hilbert_downsample_parsed},'UniformOutput',false);
    hilbert_band_ave_envelope=[hilbert_band_ave_envelope{:,:}];
    hilbert_band_ave_envelope_pre_trial=cellfun(@(x) reshape(x,size(x,1),1,size(x,2)),{nonwords.signal_ave_pre_trial_hilbert_downsample},'UniformOutput',false);
    hilbert_band_ave_envelope_pre_trial=[hilbert_band_ave_envelope_pre_trial{:,:}];
    hilbert_band_envelope=cellfun(@(x) x(1:8),{nonwords.signal_hilbert_downsample_parsed},'UniformOutput',false);
    hilbert_band_envelope=[hilbert_band_envelope{:,:}];
    hilbert_band_envelope_pre_trial=cellfun(@(x) reshape(x,size(x,1),1,size(x,2)),{nonwords.signal_pre_trial_hilbert_downsample},'UniformOutput',false);
    hilbert_band_envelope_pre_trial=[hilbert_band_envelope_pre_trial{:,:}];
    hilbert_zs_band_envelope=cellfun(@(x) x(1:8),{nonwords.signal_hilbert_zs_downsample_parsed},'UniformOutput',false);
    hilbert_zs_band_envelope=[hilbert_zs_band_envelope{:,:}];
    hilbert_zs_band_envelope_pre_trial=cellfun(@(x) reshape(x,size(x,1),1,size(x,2)),{nonwords.signal_pre_trial_hilbert_zs_downsample},'UniformOutput',false);
    hilbert_zs_band_envelope_pre_trial=[hilbert_zs_band_envelope_pre_trial{:,:}];
    try
        hilbert_pca_zs_band_envelope=cellfun(@(x) x(1:8),{nonwords.signal_hilbert_pca_zs_downsample_parsed},'UniformOutput',false);
        hilbert_pca_zs_band_envelope=[hilbert_pca_zs_band_envelope{:,:}];
        hilbert_pca_zs_band_envelope_pre_trial=cellfun(@(x) reshape(x,size(x,1),1,size(x,2)),{nonwords.signal_pre_trial_hilbert_pca_zs_downsample},'UniformOutput',false);
        hilbert_pca_zs_band_envelope_pre_trial=[hilbert_pca_zs_band_envelope_pre_trial{:,:}];
    end
 
    %
    nonwords_hilbert_band_ave_envelope_pre_trial_tensor=permute(hilbert_band_ave_envelope_pre_trial,[1,3,2]);
    nonwords_hilbert_band_envelope_pre_trial_tensor=permute(hilbert_band_envelope_pre_trial,[1,3,2]);
    nonwords_hilbert_band_ave_envelope_tensor=cell2mat(permute(hilbert_band_ave_envelope,[3,1,2]));
    nonwords_hilbert_band_envelope_tensor=cell2mat(permute(hilbert_band_envelope,[3,1,2]));
    nonwords_hilbert_zs_band_envelope_pre_trial_tensor=permute(hilbert_zs_band_envelope_pre_trial,[1,3,2]);
    nonwords_hilbert_zs_band_envelope_tensor=cell2mat(permute(hilbert_zs_band_envelope,[3,1,2]));
      try
        nonwords_hilbert_pca_zs_band_envelope_pre_trial_tensor=permute(hilbert_pca_zs_band_envelope_pre_trial,[1,3,2]);
        nonwords_hilbert_pca_zs_band_envelope_tensor=cell2mat(permute(hilbert_pca_zs_band_envelope,[3,1,2]));
    end
    %
    session_nonwords_hilbert_band_envelope_pre_trial_tensor=cat(3,session_nonwords_hilbert_band_envelope_pre_trial_tensor,nonwords_hilbert_band_envelope_pre_trial_tensor);
    session_nonwords_hilbert_band_ave_envelope_pre_trial_tensor=cat(3,session_nonwords_hilbert_band_ave_envelope_pre_trial_tensor,nonwords_hilbert_band_ave_envelope_pre_trial_tensor);
    session_nonwords_hilbert_band_envelope_tensor=cat(3,session_nonwords_hilbert_band_envelope_tensor,nonwords_hilbert_band_envelope_tensor);
    session_nonwords_hilbert_band_ave_envelope_tensor=cat(3,session_nonwords_hilbert_band_ave_envelope_tensor,nonwords_hilbert_band_ave_envelope_tensor);
    
    session_nonwords_hilbert_zs_band_envelope_pre_trial_tensor=cat(3,session_nonwords_hilbert_zs_band_envelope_pre_trial_tensor,nonwords_hilbert_zs_band_envelope_pre_trial_tensor);
    session_nonwords_hilbert_zs_band_envelope_tensor=cat(3,session_nonwords_hilbert_zs_band_envelope_tensor,nonwords_hilbert_zs_band_envelope_tensor);
    try
        session_nonwords_hilbert_pca_zs_band_envelope_pre_trial_tensor=cat(3,session_nonwords_hilbert_pca_zs_band_envelope_pre_trial_tensor,sentence_hilbert_pca_zs_band_envelope_pre_trial_tensor);
        session_nonwords_hilbert_pca_zs_band_envelope_tensor=cat(3,session_nonwords_hilbert_pca_zs_band_envelope_tensor,sentence_hilbert_pca_zs_band_envelope_tensor);
    end
    example_nonwords=cellfun(@(x) x(2:end),{nonwords(:).trial_string},'UniformOutput',false);
    nonwords_patterns=[nonwords_patterns,example_nonwords];
    
    %%%%%%%%%%%%%%% jabberwocky
    jabberwocky_trial_index=~cellfun(@isempty,cellfun(@ (x) strfind(x,'J'),info.word_type,'UniformOutput',false));
    jabberwocky=[data{jabberwocky_trial_index}];
    %
    hilbert_band_ave_envelope=cellfun(@(x) x(1:8),{jabberwocky.signal_ave_hilbert_downsample_parsed},'UniformOutput',false);
    hilbert_band_ave_envelope=[hilbert_band_ave_envelope{:,:}];
    hilbert_band_ave_envelope_pre_trial=cellfun(@(x) reshape(x,size(x,1),1,size(x,2)),{jabberwocky.signal_ave_pre_trial_hilbert_downsample},'UniformOutput',false);
    hilbert_band_ave_envelope_pre_trial=[hilbert_band_ave_envelope_pre_trial{:,:}];
    hilbert_band_envelope=cellfun(@(x) x(1:8),{jabberwocky.signal_hilbert_downsample_parsed},'UniformOutput',false);
    hilbert_band_envelope=[hilbert_band_envelope{:,:}];
    hilbert_band_envelope_pre_trial=cellfun(@(x) reshape(x,size(x,1),1,size(x,2)),{jabberwocky.signal_pre_trial_hilbert_downsample},'UniformOutput',false);
    hilbert_band_envelope_pre_trial=[hilbert_band_envelope_pre_trial{:,:}];
    hilbert_zs_band_envelope=cellfun(@(x) x(1:8),{jabberwocky.signal_hilbert_zs_downsample_parsed},'UniformOutput',false);
    hilbert_zs_band_envelope=[hilbert_zs_band_envelope{:,:}];
    hilbert_zs_band_envelope_pre_trial=cellfun(@(x) reshape(x,size(x,1),1,size(x,2)),{jabberwocky.signal_pre_trial_hilbert_zs_downsample},'UniformOutput',false);
    hilbert_zs_band_envelope_pre_trial=[hilbert_zs_band_envelope_pre_trial{:,:}];
    try
        hilbert_pca_zs_band_envelope=cellfun(@(x) x(1:8),{jabberwocky.signal_hilbert_pca_zs_downsample_parsed},'UniformOutput',false);
        hilbert_pca_zs_band_envelope=[hilbert_pca_zs_band_envelope{:,:}];
        hilbert_pca_zs_band_envelope_pre_trial=cellfun(@(x) reshape(x,size(x,1),1,size(x,2)),{jabberwocky.signal_pre_trial_hilbert_pca_zs_downsample},'UniformOutput',false);
        hilbert_pca_zs_band_envelope_pre_trial=[hilbert_pca_zs_band_envelope_pre_trial{:,:}];
    end
 
    %
    jabberwocky_hilbert_band_ave_envelope_pre_trial_tensor=permute(hilbert_band_ave_envelope_pre_trial,[1,3,2]);
    jabberwocky_hilbert_band_envelope_pre_trial_tensor=permute(hilbert_band_envelope_pre_trial,[1,3,2]);
    jabberwocky_hilbert_band_ave_envelope_tensor=cell2mat(permute(hilbert_band_ave_envelope,[3,1,2]));
    jabberwocky_hilbert_band_envelope_tensor=cell2mat(permute(hilbert_band_envelope,[3,1,2]));
    jabberwocky_hilbert_zs_band_envelope_pre_trial_tensor=permute(hilbert_zs_band_envelope_pre_trial,[1,3,2]);
    jabberwocky_hilbert_zs_band_envelope_tensor=cell2mat(permute(hilbert_zs_band_envelope,[3,1,2]));
    try
        jabberwocky_hilbert_pca_zs_band_envelope_pre_trial_tensor=permute(hilbert_pca_zs_band_envelope_pre_trial,[1,3,2]);
        jabberwocky_hilbert_pca_zs_band_envelope_tensor=cell2mat(permute(hilbert_pca_zs_band_envelope,[3,1,2]));
    end
    %
    session_jabberwocky_hilbert_band_envelope_pre_trial_tensor=cat(3,session_jabberwocky_hilbert_band_envelope_pre_trial_tensor,jabberwocky_hilbert_band_envelope_pre_trial_tensor);
    session_jabberwocky_hilbert_band_ave_envelope_pre_trial_tensor=cat(3,session_jabberwocky_hilbert_band_ave_envelope_pre_trial_tensor,jabberwocky_hilbert_band_ave_envelope_pre_trial_tensor);
    session_jabberwocky_hilbert_band_envelope_tensor=cat(3,session_jabberwocky_hilbert_band_envelope_tensor,jabberwocky_hilbert_band_envelope_tensor);
    session_jabberwocky_hilbert_band_ave_envelope_tensor=cat(3,session_jabberwocky_hilbert_band_ave_envelope_tensor,jabberwocky_hilbert_band_ave_envelope_tensor);
    session_jabberwocky_hilbert_zs_band_envelope_pre_trial_tensor=cat(3,session_jabberwocky_hilbert_zs_band_envelope_pre_trial_tensor,jabberwocky_hilbert_zs_band_envelope_pre_trial_tensor);
    session_jabberwocky_hilbert_zs_band_envelope_tensor=cat(3,session_jabberwocky_hilbert_zs_band_envelope_tensor,jabberwocky_hilbert_zs_band_envelope_tensor);
    
    try
        session_jabber_hilbert_pca_zs_band_envelope_pre_trial_tensor=cat(3,session_jabber_hilbert_pca_zs_band_envelope_pre_trial_tensor,sentence_hilbert_pca_zs_band_envelope_pre_trial_tensor);
        session_jabberwocky_hilbert_pca_zs_band_envelope_tensor=cat(3,session_jabberwocky_hilbert_pca_zs_band_envelope_tensor,sentence_hilbert_pca_zs_band_envelope_tensor);
    end
    example_jabber=cellfun(@(x) x(2:end),{jabberwocky(:).trial_string},'UniformOutput',false);
    jabber_patterns=[jabber_patterns,example_jabber];
    
    
    fprintf('added %s from %s \n',d(i).name, strcat(d(i).folder,'/',d(i).name));
    clear data subj
end
%% 
idx=17;
electrode_id=language_electrode_num(idx);

x=squeeze(session_sentence_hilbert_band_envelope_tensor(electrode_id,:,:));
x_pre=squeeze(session_sentence_hilbert_band_envelope_pre_trial_tensor(electrode_id,:,:));
elec_sentence_resp=double([x_pre;x]);
% 
x=squeeze(session_words_hilbert_band_envelope_tensor(electrode_id,:,:));
x_pre=squeeze(session_words_hilbert_band_envelope_pre_trial_tensor(electrode_id,:,:));
elec_words_resp=double([x_pre;x]);
% 
x=squeeze(session_nonwords_hilbert_band_envelope_tensor(electrode_id,:,:));
x_pre=squeeze(session_nonwords_hilbert_band_envelope_pre_trial_tensor(electrode_id,:,:));
elec_nonwords_resp=double([x_pre;x]);
%
x=squeeze(session_jabberwocky_hilbert_band_envelope_tensor(electrode_id,:,:));
x_pre=squeeze(session_jabberwocky_hilbert_band_envelope_pre_trial_tensor(electrode_id,:,:));
elec_jabber_resp=double([x_pre;x]);


%% for presentation 

colors=cbrewer('qual','Paired',10);
time_color=cbrewer('qual','Set1',10);
close all 
f=figure;
aspect_ration=9.32./4.13;
y=600;
set(f,'position',[591 455 aspect_ration*y y]);
f.Units = 'normalized';
ax=axes('position',[.05,.4,.25,.25*aspect_ration]);
fspial_lh = ft_read_headshape('/Applications/freesurfer/subjects/fsaverage/surf/lh.pial');
fspial_lh.coordsys = 'fsaverage';
h=ft_plot_mesh(fspial_lh);
h.DisplayName='Average';
h_1=ft_plot_sens(AMC026_elec_fsaverage,'elecshape','point','elecsize',15,'facecolor',colors(5,:),'edgecolor',[0,0,0]);
h_1.DisplayName='AMC026';
h_2=ft_plot_sens(AMC029_elec_fsaverage,'elecshape','point','elecsize',15,'facecolor',colors(6,:),'edgecolor',[0,0,0]);
h_2.DisplayName='AMC029';
view([-100 10]);
material dull;
lighting gouraud;
camlight; 
legend('position',[.035,.45,.04,0.05],'fontsize',12,'fontweight','bold')

%a=annotation(f,'textbox',[.1 .74 .1 .1],'String','A','FontSize',20,'FontWeight','bold');
%a.LineStyle='none';
%
k1=70;y1=double(elec_sentence_resp(:,k1));
k2=45;y2=double(elec_sentence_resp(:,k2));
% 
k3=70;y3=double(elec_words_resp(:,k3));
k4=72;y4=double(elec_words_resp(:,k4));
% 
k5=70;y5=double(elec_nonwords_resp(:,k5));
k6=72;y6=double(elec_nonwords_resp(:,k6));
% 
k7=70;y7=double(elec_jabber_resp(:,k7));
k8=72;y8=double(elec_jabber_resp(:,k8));

% 
ax_max=max([y1;y2;y3;y4;y5;y6;y7;y8]);
ax_min=min([y1;y2;y3;y4;y5;y6;y7;y8]);

%
ax=axes('position',[.35,.8,.23,0.1]);

x=[1:length(y1)]-120;
l=plot(x,y1);
hold on
l.Color=time_color(1,:)
l.LineWidth=2;
ax.YLim=[ax_min,1.2*ax_max];
ax.XLim=[-120,1080];
ax.XTick=[-120,1:135:1080];
arrayfun(@(x)plot([x,x],ax.YLim,'k--'),ax.XTick(2:end));
ticks=strsplit(sentence_patterns{k1},' ');
ticks=[{' '},ticks];
str=arrayfun(@(x) text(ax.XTick(x),ax.YLim(2),ticks(x),...
    'VerticalAlignment','top','HorizontalAlignment','left','FontSize',12,'FontWeight','bold','Rotation',45),1:9);
ax.YAxis.Visible='on';
ax.XAxis.Visible='off';
box off
ax=axes('position',[.35,.6,.23,0.1]);
l=plot(x,y2);
hold on;
l.LineWidth=2;
ax.XLim=[-120,1080];
ax.YLim=[ax_min,1.2*ax_max];
ax.XTick=[1:135:1080];
arrayfun(@(x)plot([x,x],ax.YLim,'k--'),ax.XTick);
l.Color=time_color(1,:)
ticks=strsplit(sentence_patterns{k2},' ');
str=arrayfun(@(x) text(ax.XTick(x),ax.YLim(2),ticks(x),...
    'VerticalAlignment','top','HorizontalAlignment','left','FontSize',12,'FontWeight','bold','Rotation',45),1:8)
ax.YAxis.Visible='off';
ax.XAxis.Visible='off';
box off

% words 
ax=axes('position',[.68,.8,.23,0.1]);

l=plot(x,y3);
hold on
l.LineWidth=2;
l.Color=time_color(2,:)
ax.XLim=[-120,1080];
ax.YLim=[ax_min,1.2*ax_max];
ax.XLim=[-120,1080];
ax.XTick=[1:135:1080];
arrayfun(@(x)plot([x,x],ax.YLim,'k--'),ax.XTick)
ticks=strsplit(wordlist_patterns{k3},' ')

str=arrayfun(@(x) text(ax.XTick(x),ax.YLim(2),ticks(x),'VerticalAlignment','top','HorizontalAlignment','left','FontSize',12,'FontWeight','bold','Rotation',45),1:8)
ax.YAxis.Visible='off';
ax.XAxis.Visible='off';
box off

ax=axes('position',[.68,.6,.23,0.1]);
l=plot(x,y4);
hold on
l.LineWidth=2;
l.Color=time_color(2,:)
ax.XLim=[-120,1080];
ax.YLim=[ax_min,1.2*ax_max];
ax.XLim=[-120,1080];
ax.XTick=[1:135:1080];
arrayfun(@(x)plot([x,x],ax.YLim,'k--'),ax.XTick)

ticks=strsplit(wordlist_patterns{k4},' ')
str=arrayfun(@(x) text(ax.XTick(x),ax.YLim(2),ticks(x),'VerticalAlignment','top','HorizontalAlignment','left','FontSize',12,'FontWeight','bold','Rotation',45),1:8)
ax.YAxis.Visible='off';
ax.XAxis.Visible='off';
box off

ax=axes('position',[.35,.35,.23,0.1]);

x=[1:length(y1)]-120;
l=plot(x,y5);
hold on
l.LineWidth=2;
l.Color=time_color(3,:)
ax.YLim=[ax_min,1.2*ax_max];
ax.XTick=[-120,1:135:1080];
ax.XLim=[-120,1080];
arrayfun(@(x)plot([x,x],ax.YLim,'k--'),ax.XTick(2:end));
ticks=strsplit(nonwords_patterns{k5},' ');
ticks=[{' '},ticks];
str=arrayfun(@(x) text(ax.XTick(x),ax.YLim(2),ticks(x),...
    'VerticalAlignment','top','HorizontalAlignment','left','FontSize',12,'FontWeight','bold','Rotation',45),1:9);
ax.YAxis.Visible='off';
ax.XAxis.Visible='off';

ax=axes('position',[.35,.15,.23,0.1]);

x=[1:length(y1)]-120;
l=plot(x,y6);
hold on
l.LineWidth=2;
l.Color=time_color(3,:)
ax.YLim=[ax_min,1.2*ax_max];
ax.XTick=[-120,1:135:1080];
ax.XLim=[-120,1080];
arrayfun(@(x)plot([x,x],ax.YLim,'k--'),ax.XTick(2:end));
ticks=strsplit(nonwords_patterns{k6},' ');
ticks=[{' '},ticks];
str=arrayfun(@(x) text(ax.XTick(x),ax.YLim(2),ticks(x),...
    'VerticalAlignment','top','HorizontalAlignment','left','FontSize',12,'FontWeight','bold','Rotation',45),1:9);
ax.YAxis.Visible='off';
ax.XAxis.Visible='off';


ax=axes('position',[.68,.35,.23,0.1]);

x=[1:length(y1)]-120;
l=plot(x,y7);
hold on
l.LineWidth=2;
l.Color=time_color(4,:)
ax.YLim=[ax_min,1.2*ax_max];
ax.XTick=[-120,1:135:1080];
ax.XLim=[-120,1080];
arrayfun(@(x)plot([x,x],ax.YLim,'k--'),ax.XTick(2:end));
ticks=strsplit(jabber_patterns{k7},' ');
ticks=[{' '},ticks];
str=arrayfun(@(x) text(ax.XTick(x),ax.YLim(2),ticks(x),...
    'VerticalAlignment','top','HorizontalAlignment','left','FontSize',12,'FontWeight','bold','Rotation',45),1:9);
ax.YAxis.Visible='off';
ax.XAxis.Visible='off';
uistack(ax.Children(7),'top');
ax=axes('position',[.68,.15,.23,0.1]);

x=[1:length(y1)]-120;
l=plot(x,y8);
hold on
l.LineWidth=2;
l.Color=time_color(4,:)
ax.YLim=[ax_min,1.2*ax_max];
ax.XTick=[-120,1:135:1080];
ax.XLim=[-120,1080];
arrayfun(@(x)plot([x,x],ax.YLim,'k--'),ax.XTick(2:end));
ticks=strsplit(nonwords_patterns{k8},' ');
ticks=[{' '},ticks];
str=arrayfun(@(x) text(ax.XTick(x),ax.YLim(2),ticks(x),...
    'VerticalAlignment','top','HorizontalAlignment','left','FontSize',12,'FontWeight','bold','Rotation',45),1:9);
ax.YAxis.Visible='off';
ax.XAxis.Visible='off';


%
ax=axes('position',[.05,.15,.23,.2]);
y1=double(mean(elec_sentence_resp,2));
e1=nanstd(elec_sentence_resp,0,2)./sqrt(sum(~isnan(elec_sentence_resp),2));
y2=double(mean(elec_words_resp,2));
e2=double(nanstd(elec_words_resp,0,2)./sqrt(sum(~isnan(elec_words_resp),2)));
y3=double(mean(elec_nonwords_resp,2));
e3=double(nanstd(elec_nonwords_resp,0,2)./sqrt(sum(~isnan(elec_nonwords_resp),2)));
y4=double(mean(elec_jabber_resp,2));
e4=double(nanstd(elec_jabber_resp,0,2)./sqrt(sum(~isnan(elec_jabber_resp),2)));
x=[1:length(y1)]-120;
[l,p] = boundedline(x, y1, e1, 'cmap',time_color(1,:));
hAnnotation=arrayfun(@(x) get(x,'Annotation'),p);
hLegendEntry = arrayfun(@(x) get(x,'LegendInformation'),hAnnotation);
arrayfun(@(x) set(x,'IconDisplayStyle','off'),hLegendEntry);
l.DisplayName='S';
l.LineWidth=2;

hold on 
% 
[l,p] = boundedline(x, y2, e2, 'cmap',time_color(2,:));
l.LineWidth=2;
hAnnotation=arrayfun(@(x) get(x,'Annotation'),p);
hLegendEntry = arrayfun(@(x) get(x,'LegendInformation'),hAnnotation);
arrayfun(@(x) set(x,'IconDisplayStyle','off'),hLegendEntry);
l.DisplayName='W';
% 
[l,p] = boundedline(x, y3, e3, 'cmap',time_color(3,:));
l.LineWidth=2;
hAnnotation=arrayfun(@(x) get(x,'Annotation'),p);
hLegendEntry = arrayfun(@(x) get(x,'LegendInformation'),hAnnotation);
arrayfun(@(x) set(x,'IconDisplayStyle','off'),hLegendEntry);
l.DisplayName='N';
% 
[l,p] = boundedline(x, y4, e4, 'cmap',time_color(4,:));
l.LineWidth=2;
hAnnotation=arrayfun(@(x) get(x,'Annotation'),p);
hLegendEntry = arrayfun(@(x) get(x,'LegendInformation'),hAnnotation);
arrayfun(@(x) set(x,'IconDisplayStyle','off'),hLegendEntry);
l.DisplayName='J';


uistack(ax.Children(7),'top');
leg=legend('position',[.25,.4,.03,0.05],'fontsize',12);
ax.XLim=[-120,1080];
ax.XTick=[-120,1:135:1080];
ax.FontSize=12;
ref=arrayfun(@(x)plot([x,x],ax.YLim,'k--'),ax.XTick(2:end));
Annotation=arrayfun(@(x) get(x,'Annotation'),ref);
hLegendEntry = arrayfun(@(x) get(x,'LegendInformation'),Annotation);
arrayfun(@(x) set(x,'IconDisplayStyle','off'),hLegendEntry);

ax.XTickLabel={'Word:','1','2','3','4','5','6','7','8'};

ax.YAxis.Label.String='\gamma'

ax.FontWeight='bold'

%a=annotation(f,'textbox',[.1 .32 .1 .1],'String','B','FontSize',8,'FontWeight','bold');
%a.LineStyle='none';
%a=annotation(f,'textbox',[.1 .17 .55 .1],'String',...
%    'Figure.1 (A) Example of two subjects overlayed on an average brain, (B) Response of a sample electrode to S and N conditions ',...
%    'FontSize',25);
%a.LineStyle='none';
 
if ~exist(strcat(analysis_path))
    mkdir(strcat(analysis_path))
end

print(f, '-djpeg', strcat(analysis_path,'/average_anatomicals_presentation.jpeg'));
%print(f, '-painters', '-depsc', strcat(analysis_path,'/average_anatomicals_presentation.eps'));