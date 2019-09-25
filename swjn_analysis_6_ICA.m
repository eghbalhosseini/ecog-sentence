%% STEP 0: prepare the workspace 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all 
close all 
data_path='/Users/eghbalhosseiniasl1/MyData/ecog-sentence';
subject_id='AMC';
d= dir([data_path,strcat('/**/',subject_id,'*_crunched.mat')]);
fprintf(' %d .mat files were found \n', length(d));
gamma_band_index=4;


% 
subject_ids=table2cell(unique(cell2table(cellfun(@(x) x(regexpi(x,'AMC')+[0:5]), {d.name},'UniformOutput',false)')));
% 

save_path='/Users/eghbalhosseiniasl1/MyData/ECoG-sentence/crunched/';
analysis_path='/Users/eghbalhosseiniasl1/MyData/ECoG-sentence/analysis/swjn_analysis_5_construct_RDM_for_all_subjects/';
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
d_nmerge= dir([data_path,'/**/*nmerge.txt']);
fprintf(' %d nmerge files were found \n', length(d_nmerge));
i=1;
[node_cell,node_table]=generate_node_table_eh(strcat(d_nmerge(i).folder,'/',d_nmerge(i).name));
% 
%% get wordlist
d_wordlist=dir([data_path,'/**/*words.txt']);
i=1;

fid=fopen(strcat(d_wordlist(i).folder,'/',d_wordlist(i).name));
tline = fgetl(fid);
wordlist_cell = cell(0,1);
while ischar(tline)
    wordlist_cell{end+1,1} = tline;
    tline = fgetl(fid);
end
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  STEP 1 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
find_index= @(input_cell) find(cell2mat(cellfun(@(x) ~isempty(x), input_cell,'UniformOutput',false)));
all_subjects_data={};
all_wordlists={};
all_subjects_info.fields={'subject_id';'sentence';'words';'nonwords';'jabberwocky';'sentences';'sentence_ids'};
for k=1:length(subject_ids)
    all_sentence_pattern=[];
    all_sentence_pattern_id=[];
    all_example_sentences=[];
    all_sentence_locations=[];
    all_example_wordlist=[];
    all_wordlist_locations=[];
    
    sentence_electrode_with_langauge_accross_sessions=[];
    words_electrode_with_langauge_accross_sessions=[];
    session_sentence_hilbert_band_envelope_lang_elec_tensor=[];
    session_sentence_hilbert_zs_band_envelope_lang_elec_tensor=[];
    session_nonwords_hilbert_band_envelope_lang_elec_tensor=[];
    session_words_hilbert_band_envelope_lang_elec_tensor=[];
    session_jabberwocky_hilbert_band_envelope_lang_elec_tensor=[];
    
    d_sub=(find(~cellfun(@isempty,cellfun(@(x) strfind(x,subject_ids{k}),{d.name},'UniformOutput',false))));
    for i=d_sub
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
        hilbert_band_envelope=cellfun(@(x) x(1:8),{sentences.signal_hilbert_zs_downsample_parsed},'UniformOutput',false);
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
        all_sentence_locations=[all_sentence_locations;example_sentence_locations'];
        
        %%%%%%%%%%%%%%%%% words
        words_trial_index=~cellfun(@isempty,cellfun(@ (x) strfind(x,'W'),info.word_type,'UniformOutput',false));
        % words
        words=[data{words_trial_index}];
        hilbert_band_envelope=cellfun(@(x) x(1:8),{words.signal_hilbert_zs_downsample_parsed},'UniformOutput',false);
        % creat a cell with wordposition(row)*time in trial(column) structure
        hilbert_band_envelope=[hilbert_band_envelope{:,:}];
        hilbert_band_envelope=cellfun(@transpose,hilbert_band_envelope,'UniformOutput',false);
        % make a words positions*channel* trial tensor
        hilbert_band_envelope_tensor=cell2mat(permute(hilbert_band_envelope,[1,3,2]));
        %append to langauge_channel*trial*words positions
        hilbert_band_envelope_tensor=permute(hilbert_band_envelope_tensor,[2,3,1]);
        hilbert_band_envelope_lang_elec_tensor=hilbert_band_envelope_tensor(find(language_electrode),:,:);
        session_words_hilbert_band_envelope_lang_elec_tensor=cat(2,session_words_hilbert_band_envelope_lang_elec_tensor,hilbert_band_envelope_lang_elec_tensor);
        % 
        example_wordlist=cellfun(@(x) x(2:end),{words(:).trial_string},'UniformOutput',false);
        all_example_wordlist=[all_example_wordlist;example_wordlist'];
        all_wordlists=[all_wordlists;example_wordlist'];
        
        
        fprintf('added %s from %s \n',d(i).name, strcat(d(i).folder,'/',d(i).name));
        clear data subj
    end
    all_subjects_data=[all_subjects_data;
                        [subject_ids(k),{session_sentence_hilbert_band_envelope_lang_elec_tensor},{session_words_hilbert_band_envelope_lang_elec_tensor},...
                        {session_nonwords_hilbert_band_envelope_lang_elec_tensor},{session_jabberwocky_hilbert_band_envelope_lang_elec_tensor},...
                        {all_example_sentences},{all_sentence_locations},...
                        {all_example_wordlist},...
                        {word_length}]];
end
%% 
unique_wordlist=table2cell(unique(cell2table(all_wordlists)));
for k=1:size(all_subjects_data(:,8),1)
    example_wordlist_locations=cellfun(@(x) (regexpi(unique_wordlist,x)),[all_subjects_data{k,8}],'UniformOutput',false);
    example_wordlist_locations=cell2mat(cellfun(@(x) find_index(x), example_wordlist_locations, 'UniformOutput',false));
    all_subjects_data{k,10}=example_wordlist_locations;
end 
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  find the intersections of the sentences 
all_sentence_locations=all_subjects_data(:,7);
H=all_sentence_locations{1};
for i=2:length(all_sentence_locations)
    H=intersect(H,all_sentence_locations{i});
end 

%%
%%  find the intersections of the wordlists 
all_wordlist_locations=all_subjects_data(:,10);
H1=all_wordlist_locations{1};
for i=2:length(all_wordlist_locations)
    H1=intersect(H1,all_wordlist_locations{i});
end 
%% 
close all 
look_out_window=135; %200 ms
offset=0; % 100 ms
word_range=[1:8];
electrode_words_in_sentence_mat=[];
for k=1:length(subject_ids)
    num_electrodes=size(all_subjects_data{k,2},1);
    electrode_sentence_time_tensor=all_subjects_data{k,2};
    word_length=all_subjects_data{k,9};
    subject_sentences=all_subjects_data{k,6};
    sentence_locations=cellfun(@(x) find(all_subjects_data{k,7}==x),mat2cell(H,ones(1,size(H,1)),1),'UniformOutput',false);
    shared_sentences=cellfun(@(x) subject_sentences{x(1)},sentence_locations,'UniformOutput',false);
    a=[];b=[];
    for i=1:num_electrodes
        electrode_response=squeeze(electrode_sentence_time_tensor(i,:,:));  
        electrode_sentence_time_cell=mat2cell(electrode_response,ones(1,size(electrode_response,1)),word_length*ones(1,8));
        electrode_shared_sentence_time_cell=cellfun(@(x) electrode_sentence_time_cell(x,:),sentence_locations,'UniformOutput',false);
        a=cellfun(@(x) cellfun(@(y) double(nanmean(y(offset+[1:look_out_window]))),x),electrode_shared_sentence_time_cell,'UniformOutput',false);
        b=cell2mat(cellfun(@(x) nanmean(x,1),a,'UniformOutput',false));

        electrode_words_in_sentence_mat=[electrode_words_in_sentence_mat;reshape(b,1,[])];
    end 
end 


sentence_RDM=squareform(pdist(electrode_words_in_sentence_mat','correlation'));
save(strcat(analysis_path,'/','electrode_words_in_sentence_matrix.mat'),'electrode_words_in_sentence_mat');
%% 

% dimensionality of data and components
M = 98; % number of features (e.g. sounds)
N = 416; % number of measures (e.g. fMRI voxels)
K = 10; % number of components

% create the data matrix
% decomposition analysis
N_RANDOM_INITS = 20;
PLOT_FIGURES = 1;
RAND_SEED = 1;
[R_inferred, W_inferred] = nonparametric_ica(electrode_words_in_sentence_mat, K, N_RANDOM_INITS, PLOT_FIGURES, RAND_SEED);

figure; 
set(gcf,'position',[1376 353 988 992])
axes('position',[.1,.6,.5,.3])
colors = cbrewer('div', 'RdBu', 128);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);

imagesc(electrode_words_in_sentence_mat)
set(gca, 'ydir', 'reverse','box','off');
xlabel('Sentences')
ylabel('Electrodes')
title('Response')
% 
axes('position',[.1,.1,.5,.2])
colors = cbrewer('div', 'RdBu', 128);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);

imagesc(W_inferred)
set(gca, 'ydir', 'reverse','box','off');
title('Weights')

%
axes('position',[.65,.1,.3,.5])
colors = cbrewer('div', 'RdBu', 128);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);

imagesc(R_inferred)
set(gca, 'ydir', 'reverse','box','off');
title('Features')
