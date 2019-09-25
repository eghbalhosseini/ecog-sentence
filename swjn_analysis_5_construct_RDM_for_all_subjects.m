%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    addpath('~/MyCodes/ecog-sentence/');
    addpath(genpath('~/MyCodes/basic-ecog-tools/activeBrain'));
    addpath(genpath('~/MyCodes/basic-ecog-tools/ecog-filters'));
    addpath(genpath('~/MyCodes/basic-ecog-tools/mex'));
end 
%% 
d_nmerge= dir([data_path,'/**/*nmerge.txt']);
fprintf(' %d nmerge files were found \n', length(d_nmerge));
i=1;
[node_cell,node_table]=generate_node_table_eh(strcat(d_nmerge(i).folder,'/',d_nmerge(i).name));

%% get wordlist
d_wordlist=dir([data_path,'/**/*words.txt']);
i=1;

fid=fopen(strcat(d_wordlist(i).folder,'/',d_wordlist(i).name));
tline = fgetl(fid);
wordlist_cell = cell(0,1);
while ischar(tline)
%    try 
%         tline = strrep(tline," '","'");
%     end 
%     try 
%         tline = strrep(tline,'can not','cannot');  
%     end 
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
emb = fastTextWordEmbedding;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% construct electrode by word matrices for sentences 
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

shared_word_2_vec_sentences=cell2mat(cellfun(@(y) double(cell2mat(cellfun(@(x) word2vec(emb,lower(x)),strsplit(strrep(y,"'S",""),' '),'UniformOutput',false)')),shared_sentences,'UniformOutput',false));

shared_sentence_words=(cellfun(@(y) strsplit(strrep(y,"'S",""),' '),shared_sentences,'UniformOutput',false));
shared_sentence_words=[shared_sentence_words{:}]';
%
Y_sentences = pdist(shared_word_2_vec_sentences,'cosine');

sentence_RDM=squareform(pdist(electrode_words_in_sentence_mat','correlation'));
save(strcat(analysis_path,'/','electrode_words_in_sentence_matrix.mat'),'electrode_words_in_sentence_mat');

tria_ecog_sentence_RDM=triu(sentence_RDM,1);
tria_ones=triu(ones(size(tria_bert_sentence_RDM)),1);
tria_ones(tria_ones==0)=nan;

tria_ecog_sentence_RDM=tria_ecog_sentence_RDM.*tria_ones;
flatten_ecog_sentence_RDM=tria_ecog_sentence_RDM(:);
flatten_ecog_sentence_RDM(isnan(flatten_ecog_sentence_RDM))=[];


% wordlists 
electrode_words_in_wordlist_mat=[];
for k=1:length(subject_ids)
    num_electrodes=size(all_subjects_data{k,2},1);
    electrode_wordlist_time_tensor=all_subjects_data{k,3};
    word_length=all_subjects_data{k,9};
    subject_sentences=all_subjects_data{k,8};
    wordlist_locations=cellfun(@(x) find(all_subjects_data{k,10}==x),mat2cell(H1,ones(1,size(H1,1)),1),'UniformOutput',false);
    shared_wordlist=cellfun(@(x) subject_sentences{x(1)},sentence_locations,'UniformOutput',false);
    a=[];b=[];
    for i=1:num_electrodes
        electrode_response=squeeze(electrode_wordlist_time_tensor(i,:,:));
        electrode_wordlist_time_cell=mat2cell(electrode_response,ones(1,size(electrode_response,1)),word_length*ones(1,8));
        electrode_shared_wordlist_time_cell=cellfun(@(x) electrode_wordlist_time_cell(x,:),wordlist_locations,'UniformOutput',false);
        a=cellfun(@(x) cellfun(@(y) double(nanmean(y(offset+[1:look_out_window]))),x),electrode_shared_wordlist_time_cell,'UniformOutput',false);
        b=cell2mat(cellfun(@(x) nanmean(x,1),a,'UniformOutput',false));
        electrode_words_in_wordlist_mat=[electrode_words_in_wordlist_mat;reshape(b,1,[])];
    end 
end 

shared_word_2_vec_wordlist=cell2mat(cellfun(@(y) double(cell2mat(cellfun(@(x) word2vec(emb,lower(x)),strsplit(strrep(strrep(y,"TASHA","SASHA"),"'S",""),' '),'UniformOutput',false)')),shared_wordlist,'UniformOutput',false));

shared_wordlist_words=(cellfun(@(y) strsplit(strrep(strrep(y,"TASHA","SASHA"),"'S",""),' '),shared_wordlist,'UniformOutput',false));
shared_wordlist_words=[shared_wordlist_words{:}]';
%
Y_wordlist = pdist(shared_word_2_vec_wordlist,'cosine');
save(strcat(analysis_path,'/','electrode_words_in_wordlist_matrix.mat'),'electrode_words_in_wordlist_mat');
wordlist_RDM=squareform(pdist(electrode_words_in_wordlist_mat','correlation'));

tria_ecog_wordlist_RDM=triu(wordlist_RDM,1);
tria_ecog_wordlist_RDM=tria_ecog_wordlist_RDM.*tria_ones;
flatten_ecog_wordlist_RDM=tria_ecog_wordlist_RDM(:);
flatten_ecog_wordlist_RDM(isnan(flatten_ecog_wordlist_RDM))=[];

% plot both rdms 

min_ax=min([sentence_RDM(:);wordlist_RDM(:)]);
max_ax=max([sentence_RDM(:);wordlist_RDM(:)]);
 
figure;
set(gcf,'position',[808 19 752 1319])

subplot(2,2,1)
colors = cbrewer('div', 'RdBu', 128);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);
imagesc(sentence_RDM,[min_ax,max_ax]);
daspect([1,1,1])
set(gca, 'ydir', 'reverse','box','off');
title('Sentences');
handles = colorbar;
handles.TickDirection = 'out';
handles.Box = 'off';
drawnow;
% 
subplot(2,2,2)
colors = cbrewer('div', 'RdBu', 128);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);
imagesc(squareform(Y_sentences));
daspect([1,1,1])
set(gca, 'ydir', 'reverse','box','off');
title('Sentences:distance between words in word2vec');
handles = colorbar;
handles.TickDirection = 'out';
handles.Box = 'off';
drawnow;
%
subplot(2,2,3)
colors = cbrewer('div', 'RdBu', 128);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);
imagesc(wordlist_RDM,[min_ax,max_ax]);
daspect([1,1,1])
set(gca, 'ydir', 'reverse','box','off');
title('Words');
handles = colorbar;
handles.TickDirection = 'out';
handles.Box = 'off';
drawnow;

subplot(2,2,4)
colors = cbrewer('div', 'RdBu', 128);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);
imagesc(squareform(Y_wordlist));
daspect([1,1,1])
set(gca, 'ydir', 'reverse','box','off');
title('Wordlists:distance between words in word2vec');
handles = colorbar;
handles.TickDirection = 'out';
handles.Box = 'off';
drawnow;


if ~exist(analysis_path)
    mkdir(analysis_path)
end

print(gcf,'-fillpage', '-dpdf', strcat(analysis_path,'/','ECoG_RDM','_window_450.pdf'));


%% do the RDM based on sorted sentences and wordlist 
%sentence
shared_word_2_vec_sentences=cell2mat(cellfun(@(y) double(cell2mat(cellfun(@(x) word2vec(emb,lower(x)),strsplit(strrep(y,"'S",""),' '),'UniformOutput',false)')),shared_sentences,'UniformOutput',false));

shared_sentence_words=(cellfun(@(y) strsplit(strrep(y,"'S",""),' '),shared_sentences,'UniformOutput',false));
shared_sentence_words=[shared_sentence_words{:}]';
%
Y = pdist(shared_word_2_vec_sentences,'cosine');
Z = linkage(Y,'average');
leafOrder=optimalleaforder(Z,Y);
Y_order_sentence = pdist(shared_word_2_vec_sentences(leafOrder,:),'cosine');

sentence_RDM_ordered=squareform(pdist(electrode_words_in_sentence_mat(:,leafOrder)','correlation'));
% wordlist
shared_word_2_vec_wordlist=cell2mat(cellfun(@(y) double(cell2mat(cellfun(@(x) word2vec(emb,lower(x)),strsplit(strrep(strrep(y,"TASHA","SASHA"),"'S",""),' '),'UniformOutput',false)')),shared_wordlist,'UniformOutput',false));

shared_wordlist_words=(cellfun(@(y) strsplit(strrep(strrep(y,"TASHA","SASHA"),"'S",""),' '),shared_wordlist,'UniformOutput',false));
shared_wordlist_words=[shared_wordlist_words{:}]';
%
Y = pdist(shared_word_2_vec_wordlist,'cosine');
Z = linkage(Y,'average');
leafOrder=optimalleaforder(Z,Y);
Y_order_wordlist = pdist(shared_word_2_vec_wordlist(leafOrder,:),'cosine');

wordlist_RDM_ordered=squareform(pdist(electrode_words_in_wordlist_mat(:,leafOrder)','correlation'));

min_ax=min([wordlist_RDM_ordered(:);sentence_RDM_ordered(:)]);
max_ax=max([wordlist_RDM_ordered(:);sentence_RDM_ordered(:)]);
 
figure;
set(gcf,'position',[808 19 752 1319])

subplot(2,2,1)
colors = cbrewer('div', 'RdBu', 128);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);
imagesc(sentence_RDM_ordered,[min_ax,max_ax]);
daspect([1,1,1])
set(gca, 'ydir', 'reverse','box','off');
title('Sentences');
handles = colorbar;
handles.TickDirection = 'out';
handles.Box = 'off';
drawnow;
% 
subplot(2,2,2)
colors = cbrewer('div', 'RdBu', 128);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);
imagesc(squareform(Y_order_sentence));
daspect([1,1,1])
set(gca, 'ydir', 'reverse','box','off');
title('Sentences:distance between words in word2vec (grouped)');
handles = colorbar;
handles.TickDirection = 'out';
handles.Box = 'off';
drawnow;
%
subplot(2,2,3)
colors = cbrewer('div', 'RdBu', 128);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);
imagesc(wordlist_RDM_ordered,[min_ax,max_ax]);
daspect([1,1,1])
set(gca, 'ydir', 'reverse','box','off');
title('Words');
handles = colorbar;
handles.TickDirection = 'out';
handles.Box = 'off';
drawnow;

subplot(2,2,4)
colors = cbrewer('div', 'RdBu', 128);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);
imagesc(squareform(Y_order_wordlist));
daspect([1,1,1])
set(gca, 'ydir', 'reverse','box','off');
title('Wordlists:distance between words in word2vec (grouped)');
handles = colorbar;
handles.TickDirection = 'out';
handles.Box = 'off';
drawnow;

if ~exist(analysis_path)
    mkdir(analysis_path)
end

print(gcf,'-fillpage', '-dpdf', strcat(analysis_path,'/','ECoG_RDM_sorted','_window_450.pdf'));

%% compare against network 
% load the pickle data 
d_model= dir([data_path,'/models_output/*_model.mat']);
fprintf(' %d .pkl files were found \n', length(d_model));
% bert base 
i=1 ;
S_bert_cell=load([data_path,'/models_output/',d_model(i).name]);
S_NN_cell=S_bert_cell.data;
S_NN_cell_cut=cellfun(@(x) x(1:8,:),S_NN_cell,'UniformOutput',false);
NN_sentences=S_NN_cell_cut(H);
NN_sentences=cellfun(@double,NN_sentences,'UniformOutput',false);
NN_sentnce_mat=cell2mat(NN_sentences');
bert_sentence_RDM=squareform(pdist(NN_sentnce_mat,'correlation'));

i=2 ;
S_gpt2_cell=load([data_path,'/models_output/',d_model(i).name]);
S_NN_cell=S_gpt2_cell.data;
S_NN_cell_cut=cellfun(@(x) x(1:8,:),S_NN_cell,'UniformOutput',false);
NN_sentences=S_NN_cell_cut(H);
NN_sentences=cellfun(@double,NN_sentences,'UniformOutput',false);
NN_sentnce_mat=cell2mat(NN_sentences');
gpt2_sentence_RDM=squareform(pdist(NN_sentnce_mat,'correlation'));
i=3 ;
S_roberta_cell=load([data_path,'/models_output/',d_model(i).name]);
S_NN_cell=S_roberta_cell.data;
S_NN_cell_cut=cellfun(@(x) x(1:8,:),S_NN_cell,'UniformOutput',false);
NN_sentences=S_NN_cell_cut(H);
NN_sentences=cellfun(@double,NN_sentences,'UniformOutput',false);
NN_sentnce_mat=cell2mat(NN_sentences');
roberta_sentence_RDM=squareform(pdist(NN_sentnce_mat,'correlation'));
i=4 ;
S_transfo_cell=load([data_path,'/models_output/',d_model(i).name]);
S_NN_tensor=S_transfo_cell.data;
NN_sentences=S_NN_tensor(H,:,:);
NN_sentences=double(NN_sentences);
NN_sentences_mat=reshape(NN_sentences,size(NN_sentences,1)*size(NN_sentences,2),[]);
transfo_sentence_RDM=squareform(pdist(NN_sentences_mat,'correlation'));
% 
i=5 ;
W_bert_cell=load([data_path,'/models_output/',d_model(i).name]);
W_NN_cell=W_bert_cell.data;
W_NN_cell_cut=cellfun(@(x) x(1:8,:),W_NN_cell,'UniformOutput',false);
NN_wordlists=W_NN_cell_cut(H);
NN_wordlists=cellfun(@double,NN_wordlists,'UniformOutput',false);
NN_wordlist_mat=cell2mat(NN_wordlists');
bert_wordlist_RDM=squareform(pdist(NN_wordlist_mat,'correlation'));
i=6 ;
W_gpt2_cell=load([data_path,'/models_output/',d_model(i).name]);
W_NN_cell=W_gpt2_cell.data;
W_NN_cell_cut=cellfun(@(x) x(1:8,:),W_NN_cell,'UniformOutput',false);
NN_wordlists=W_NN_cell_cut(H);
NN_wordlists=cellfun(@double,NN_wordlists,'UniformOutput',false);
NN_wordlist_mat=cell2mat(NN_wordlists');
gpt2_wordlist_RDM=squareform(pdist(NN_wordlist_mat,'correlation'));
i=7 ;
W_roberta_cell=load([data_path,'/models_output/',d_model(i).name]);
W_NN_cell=W_roberta_cell.data;
W_NN_cell_cut=cellfun(@(x) x(1:8,:),W_NN_cell,'UniformOutput',false);
NN_wordlists=W_NN_cell_cut(H);
NN_wordlists=cellfun(@double,NN_wordlists,'UniformOutput',false);
NN_wordlist_mat=cell2mat(NN_wordlists');
roberta_wordlist_RDM=squareform(pdist(NN_wordlist_mat,'correlation'));
i=8 ;
W_transfo_cell=load([data_path,'/models_output/',d_model(i).name]);
W_NN_tensor=W_transfo_cell.data;
NN_wordlists=W_NN_tensor(H,:,:);
NN_wordlists=double(NN_wordlists);
NN_wordlist_mat=reshape(NN_wordlists,size(NN_wordlists,1)*size(NN_wordlists,2),[]);
transfo_wordlist_RDM=squareform(pdist(NN_wordlist_mat,'correlation'));


figure;
set(gcf,'position',[808 19 752 1319])

subplot(4,2,1)
colors = cbrewer('div', 'RdBu', 128);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);
imagesc(bert_sentence_RDM);
daspect([1,1,1])
set(gca, 'ydir', 'reverse','box','off');
title('bert-Sentences');
handles = colorbar;
handles.TickDirection = 'out';
handles.Box = 'off';
drawnow;

subplot(4,2,2)
colors = cbrewer('div', 'RdBu', 128);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);
imagesc(bert_wordlist_RDM);
daspect([1,1,1])
set(gca, 'ydir', 'reverse','box','off');
title('bert-words');
handles = colorbar;
handles.TickDirection = 'out';
handles.Box = 'off';
drawnow;


subplot(4,2,3)
colors = cbrewer('div', 'RdBu', 128);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);
imagesc(gpt2_sentence_RDM);
daspect([1,1,1])
set(gca, 'ydir', 'reverse','box','off');
title('gpt2-Sentences');
handles = colorbar;
handles.TickDirection = 'out';
handles.Box = 'off';
drawnow;

subplot(4,2,4)
colors = cbrewer('div', 'RdBu', 128);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);
imagesc(gpt2_wordlist_RDM);
daspect([1,1,1])
set(gca, 'ydir', 'reverse','box','off');
title('gpt2-words');
handles = colorbar;
handles.TickDirection = 'out';
handles.Box = 'off';
drawnow;


subplot(4,2,5)
colors = cbrewer('div', 'RdBu', 128);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);
imagesc(roberta_sentence_RDM);
daspect([1,1,1])
set(gca, 'ydir', 'reverse','box','off');
title('Roberta-Sentences');
handles = colorbar;
handles.TickDirection = 'out';
handles.Box = 'off';
drawnow;

subplot(4,2,6)
colors = cbrewer('div', 'RdBu', 128);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);
imagesc(roberta_wordlist_RDM);
daspect([1,1,1])
set(gca, 'ydir', 'reverse','box','off');
title('Roberta-words');
handles = colorbar;
handles.TickDirection = 'out';
handles.Box = 'off';
drawnow;


subplot(4,2,7)
colors = cbrewer('div', 'RdBu', 128);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);
imagesc(transfo_sentence_RDM);
daspect([1,1,1])
set(gca, 'ydir', 'reverse','box','off');
title('transfo-Sentences');
handles = colorbar;
handles.TickDirection = 'out';
handles.Box = 'off';
drawnow;


subplot(4,2,8)
colors = cbrewer('div', 'RdBu', 128);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);
imagesc(transfo_wordlist_RDM);
daspect([1,1,1])
set(gca, 'ydir', 'reverse','box','off');
title('transfo-words');
handles = colorbar;
handles.TickDirection = 'out';
handles.Box = 'off';
drawnow;

print(gcf,'-fillpage', '-dpdf', strcat(analysis_path,'/','NN_RDM','.pdf'));

%% do correlation between RDMS 

tria_bert_sentence_RDM=triu(bert_sentence_RDM,1);
tria_gpt2_sentence_RDM=triu(gpt2_sentence_RDM,1);
tria_roberta_sentence_RDM=triu(roberta_sentence_RDM,1);
tria_transfo_sentence_RDM=triu(transfo_sentence_RDM,1);

% 
tria_ones=triu(ones(size(tria_bert_sentence_RDM)),1);
tria_ones(tria_ones==0)=nan;
tria_bert_sentence_RDM=tria_bert_sentence_RDM.*tria_ones;
flatten_bert_sentence_RDM=tria_bert_sentence_RDM(:);
flatten_bert_sentence_RDM(isnan(flatten_bert_sentence_RDM))=[];

tria_gpt2_sentence_RDM=tria_gpt2_sentence_RDM.*tria_ones;
flatten_gpt2_sentence_RDM=tria_gpt2_sentence_RDM(:);
flatten_gpt2_sentence_RDM(isnan(flatten_gpt2_sentence_RDM))=[];

tria_roberta_sentence_RDM=tria_roberta_sentence_RDM.*tria_ones;
flatten_roberta_sentence_RDM=tria_roberta_sentence_RDM(:);
flatten_roberta_sentence_RDM(isnan(flatten_roberta_sentence_RDM))=[];

tria_transfo_sentence_RDM=tria_transfo_sentence_RDM.*tria_ones;
flatten_transfo_sentence_RDM=tria_transfo_sentence_RDM(:);
flatten_transfo_sentence_RDM(isnan(flatten_transfo_sentence_RDM))=[];

%% 
tria_bert_wordlist_RDM=triu(bert_wordlist_RDM,1);
tria_gpt2_wordlist_RDM=triu(gpt2_wordlist_RDM,1);
tria_roberta_wordlist_RDM=triu(roberta_wordlist_RDM,1);
tria_transfo_wordlist_RDM=triu(transfo_wordlist_RDM,1);

tria_bert_wordlist_RDM=tria_bert_wordlist_RDM.*tria_ones;
flatten_bert_wordlist_RDM=tria_bert_wordlist_RDM(:);
flatten_bert_wordlist_RDM(isnan(flatten_bert_wordlist_RDM))=[];

tria_gpt2_wordlist_RDM=tria_gpt2_wordlist_RDM.*tria_ones;
flatten_gpt2_wordlist_RDM=tria_gpt2_wordlist_RDM(:);
flatten_gpt2_wordlist_RDM(isnan(flatten_gpt2_wordlist_RDM))=[];

tria_roberta_wordlist_RDM=tria_roberta_wordlist_RDM.*tria_ones;
flatten_roberta_wordlist_RDM=tria_roberta_wordlist_RDM(:);
flatten_roberta_wordlist_RDM(isnan(flatten_roberta_wordlist_RDM))=[];

tria_transfo_wordlist_RDM=tria_transfo_wordlist_RDM.*tria_ones;
flatten_transfo_wordlist_RDM=tria_transfo_wordlist_RDM(:);
flatten_transfo_wordlist_RDM(isnan(flatten_transfo_wordlist_RDM))=[];

%% compute a correlation between ecog data and networks 
R_ecog_bert_sentence = corrcoef(flatten_ecog_sentence_RDM,flatten_bert_sentence_RDM);
R_ecog_gpt2_sentence = corrcoef(flatten_ecog_sentence_RDM,flatten_gpt2_sentence_RDM);
R_ecog_roberta_sentence = corrcoef(flatten_ecog_sentence_RDM,flatten_roberta_sentence_RDM);
R_ecog_transfo_sentence = corrcoef(flatten_ecog_sentence_RDM,flatten_transfo_sentence_RDM);
% 
R_ecog_bert_wordlist = corrcoef(flatten_ecog_wordlist_RDM,flatten_bert_wordlist_RDM);
R_ecog_gpt2_wordlist = corrcoef(flatten_ecog_wordlist_RDM,flatten_gpt2_wordlist_RDM);
R_ecog_roberta_wordlist = corrcoef(flatten_ecog_wordlist_RDM,flatten_roberta_wordlist_RDM);
R_ecog_transfo_wordlist = corrcoef(flatten_ecog_wordlist_RDM,flatten_transfo_wordlist_RDM);

y=[R_ecog_bert_sentence(1,2),R_ecog_bert_wordlist(1,2);...
   R_ecog_gpt2_sentence(1,2),R_ecog_gpt2_wordlist(1,2);...
   R_ecog_roberta_sentence(1,2),R_ecog_roberta_wordlist(1,2);...
   R_ecog_transfo_sentence(1,2),R_ecog_transfo_wordlist(1,2)];
% 
figure;
bl=bar(y,'Facecolor','flat');
set(bl(1),'Linestyle','none','BarWidth',.8);
set(bl(2),'Linestyle','none','BarWidth',.8);
set(bl(1),'Displayname','Sentences');
set(bl(2),'Displayname','Wordlist');
set(gca,'xticklabel',{'bert','gpt2','roberta','transfo'});
set(gca,'FontSize',14)
set(gca,'box','off')
legend('location','northeastoutside')
title('Correlation between ECoG RDM and Network RDMs');
ylabel('Correlation');

print(gcf,'-bestfit', '-dpdf', strcat(analysis_path,'/','ECoG_NN_RDM_correlation','.pdf'));
