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
    addpath(genpath('~/MyCodes/ecog-sentence/'));
    addpath(genpath('~/MyCodes/basic-ecog-tools/activeBrain'));
    addpath(genpath('~/MyCodes/basic-ecog-tools/ecog-filters'));
    addpath(genpath('~/MyCodes/basic-ecog-tools/mex'));
end


find_index= @(input_cell) find(cell2mat(cellfun(@(x) ~isempty(x), input_cell,'UniformOutput',false)));
electrode_group='language';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  STEP 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

all_sub_dat=[];

for k=1:length(subject_ids)-1
    session_sentence_examples={};
    session_wordlist_examples={};
    session_nonwords_examples={};
    session_jabberwocky_examples={};
    session_sentence_hilbert_band_envelope_tensor=[];
    session_words_hilbert_band_envelope_tensor=[];
    session_nonwords_hilbert_band_envelope_tensor=[];
    session_jabberwocky_hilbert_band_envelope_tensor=[];
    
    d_sub=(find(~cellfun(@isempty,cellfun(@(x) strfind(x,subject_ids{k}),{d.name},'UniformOutput',false))));
    for i=d_sub
        subj=load(strcat(d(i).folder,'/',d(i).name));
        subj_id=fieldnames(subj);
        subj=subj.(subj_id{1});
        data=subj.data;
        info=subj.info;
        if strcmp(electrode_group,'language');
            language_electrode=info.language_responsive_electrodes_hilbert_odd;
        elseif strcmp(electrode_group,'buildup');
            language_electrode=info.ramp_electrodes_hilbert_odd;
        end
        language_electrode_num=find(language_electrode);
        % step 1: extract electrodes with siginificant language response
        sentence_trial_index=~cellfun(@isempty,cellfun(@ (x) strfind(x,'S'),info.word_type,'UniformOutput',false));
        %%%%%%%%%%%%%%% sentence
        hilbert_band_envelope=[];
        sentences=[data{sentence_trial_index}];%
        word_length=sentences(1).signal_range_downsample(1,2)-sentences(1).signal_range_downsample(1,1)+1;
        % creat a cell with wordposition(row)*time in trial(column) structure
        % hilbert mean (changlab)
        hilbert_band_envelope=cellfun(@(x) x(1:8),{sentences.signal_hilbert_downsample_parsed},'UniformOutput',false);
        hilbert_band_envelope=[hilbert_band_envelope{:,:}];
        hilbert_band_envelope=cellfun(@transpose,hilbert_band_envelope,'UniformOutput',false);
        % make a words positions*channel* trial tensor
        hilbert_band_envelope_tensor=cell2mat(permute(hilbert_band_envelope,[1,3,2]));
        %append to langauge_channel*trial*words positions
        hilbert_band_envelope_tensor=permute(hilbert_band_envelope_tensor,[2,3,1]);
        hilbert_band_envelope_tensor=hilbert_band_envelope_tensor(find(language_electrode),:,:);
        
        session_sentence_hilbert_band_envelope_tensor=cat(2,session_sentence_hilbert_band_envelope_tensor,hilbert_band_envelope_tensor);
        % find sentence locations
        example_sentence=cellfun(@(x) x(2:end),{sentences(:).trial_string},'UniformOutput',false);
        session_sentence_examples=[session_sentence_examples;example_sentence'];
        %%%%%%%%%%%%%%%%% words
        hilbert_band_envelope=[];
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
        hilbert_band_envelope_tensor=hilbert_band_envelope_tensor(find(language_electrode),:,:);
        session_words_hilbert_band_envelope_tensor=cat(2,session_words_hilbert_band_envelope_tensor,hilbert_band_envelope_tensor);
        %
        example_words=cellfun(@(x) x(2:end),{words(:).trial_string},'UniformOutput',false);
        session_wordlist_examples=[session_wordlist_examples;example_words'];
        %%%%%%%%%%%%%%%%%%% nonwords
        hilbert_band_envelope=[];
        nonwords_trial_index=~cellfun(@isempty,cellfun(@ (x) strfind(x,'N'),info.word_type,'UniformOutput',false));
        % nonwords
        nonwords=[data{nonwords_trial_index}];
        hilbert_band_envelope=cellfun(@(x) x(1:8),{nonwords.signal_hilbert_downsample_parsed},'UniformOutput',false);
        % creat a cell with wordposition(row)*time in trial(column) structure
        hilbert_band_envelope=[hilbert_band_envelope{:,:}];
        hilbert_band_envelope=cellfun(@transpose,hilbert_band_envelope,'UniformOutput',false);
        % make a nonwords positions*channel* trial tensor
        hilbert_band_envelope_tensor=cell2mat(permute(hilbert_band_envelope,[1,3,2]));
        %append to langauge_channel*trial*nonwords positions
        hilbert_band_envelope_tensor=permute(hilbert_band_envelope_tensor,[2,3,1]);
        hilbert_band_envelope_tensor=hilbert_band_envelope_tensor(find(language_electrode),:,:);
        session_nonwords_hilbert_band_envelope_tensor=cat(2,session_nonwords_hilbert_band_envelope_tensor,hilbert_band_envelope_tensor);
        %
        example_nonwords=cellfun(@(x) x(2:end),{nonwords(:).trial_string},'UniformOutput',false);
        session_nonwords_examples=[session_nonwords_examples;example_nonwords'];
        %
        hilbert_band_envelope=[];
        jabberwocky_trial_index=~cellfun(@isempty,cellfun(@ (x) strfind(x,'J'),info.word_type,'UniformOutput',false));
        % jabberwocky
        jabberwocky=[data{jabberwocky_trial_index}];
        hilbert_band_envelope=cellfun(@(x) x(1:8),{jabberwocky.signal_hilbert_downsample_parsed},'UniformOutput',false);
        % creat a cell with wordposition(row)*time in trial(column) structure
        hilbert_band_envelope=[hilbert_band_envelope{:,:}];
        hilbert_band_envelope=cellfun(@transpose,hilbert_band_envelope,'UniformOutput',false);
        % make a jabberwocky positions*channel* trial tensor
        hilbert_band_envelope_tensor=cell2mat(permute(hilbert_band_envelope,[1,3,2]));
        %append to langauge_channel*trial*jabberwocky positions
        hilbert_band_envelope_tensor=permute(hilbert_band_envelope_tensor,[2,3,1]);
        hilbert_band_envelope_tensor=hilbert_band_envelope_tensor(find(language_electrode),:,:);
        session_jabberwocky_hilbert_band_envelope_tensor=cat(2,session_jabberwocky_hilbert_band_envelope_tensor,hilbert_band_envelope_tensor);
        %
        example_jabberwocky=cellfun(@(x) x(2:end),{jabberwocky(:).trial_string},'UniformOutput',false);
        session_jabberwocky_examples=[session_jabberwocky_examples;example_jabberwocky'];
        %
        fprintf('added %s from %s \n',d(i).name, strcat(d(i).folder,'/',d(i).name));
        clear data subj
    end
    sub.sub_id=subject_ids{k};
    % signal
    sub.sess_sentence_hilb_tensor=session_sentence_hilbert_band_envelope_tensor;
    sub.sess_wordlist_hilb_tensor=session_words_hilbert_band_envelope_tensor;
    sub.sess_nonwords_hilb_tensor=session_nonwords_hilbert_band_envelope_tensor;
    sub.sess_jabberwocky_hilb_tensor=session_jabberwocky_hilbert_band_envelope_tensor;
    % stimuli
    sub.sess_sentence_stim=session_sentence_examples;
    sub.sess_wordlist_stim=session_wordlist_examples;
    sub.sess_nonwords_stim=session_nonwords_examples;
    sub.sess_jabberwocky_stim=session_jabberwocky_examples;
    sub.sess_stim_length=word_length;
    all_sub_dat=[all_sub_dat;sub];
    
    
    
end
%% 
%% find shared stimuli 
stim_type={'sentence','wordlist','nonwords','jabberwocky'};
for k=1:size(stim_type,2)
    fprintf(['adding ',stim_type{k},'\n'])
    eval(strcat('A=all_sub_dat(1).sess_',stim_type{k},'_stim;')) ;
    for kk=1:size(all_sub_dat,1)
        eval(strcat('B=all_sub_dat(',num2str(kk),').sess_',stim_type{k},'_stim;')) ;
        A=intersect(A,B);
    end 
    for kk=1:size(all_sub_dat,1)
        eval(strcat('all_sub_dat(kk).sess_',stim_type{k},'_stim_shared=A;')) ;
    end
end
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
A=cell2table(all_sub_dat(1).sess_sentence_stim_shared)
writetable(A,'shared_sentence.txt')
%% 
emb = fastTextWordEmbedding;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
sub_rdm_dat=struct;
find_first_index= @(input_cell,val) find(cell2mat(cellfun(@(x) ~isempty(x), input_cell,'UniformOutput',false)),val,'first');
stim_type={'sentence','wordlist','nonwords','jabberwocky'};
for k=1:size(stim_type,2)
    fprintf(['adding ',stim_type{k},' electrode-stim mat\n'])
    elec_stim_mat=[];
    elec_stim_str=[];
    elec_stim_wordpos=[];
    for kk=1:size(all_sub_dat,1)
        eval(strcat('A=all_sub_dat(',num2str(kk),').sess_',stim_type{k},'_stim;'));% get_stim
        eval(strcat('B=all_sub_dat(',num2str(kk),').sess_',stim_type{k},'_stim_shared;'));% get_stim_shared
        idx=cellfun(@(x) (regexpi(A,x)),B,'UniformOutput',false);
        idx_1st=cellfun(@(x) find_first_index(x,1), idx, 'UniformOutput',false);
        idx=cellfun(@(x) find_index(x), idx, 'UniformOutput',false);
        % get hilbert
        eval(strcat('Hilb=all_sub_dat(',num2str(kk),').sess_',stim_type{k},'_hilb_tensor;'));% get_stim_shared
        eval(strcat('Stim_L=all_sub_dat(',num2str(kk),').sess_stim_length;'));% get_word_length
        El=cellfun(@(x) double(squeeze(Hilb(:,x,:))),idx_1st','UniformOutput',false);
        El_demeaned=cellfun(@(x) transpose(bsxfun(@minus, x', mean(x'))),El,'UniformOutput',false);
        El_cut_demeaned=[];
        El_cut=[];
        a=size(El_demeaned{1},1);
        for p=1:size(El_demeaned,2)
            El_cut=[El_cut,mat2cell(El{p},a,Stim_L*ones(1,8))];
            El_cut_demeaned=[El_cut_demeaned,mat2cell(El_demeaned{p},a,Stim_L*ones(1,8))];
        end
        El_cut_ave=cell2mat(cellfun(@(x) nanmean(x,2),El_cut,'UniformOutput',false));
        El_cut_demeaned_ave=cell2mat(cellfun(@(x) nanmean(x,2),El_cut_demeaned,'UniformOutput',false));
        El_string=cellfun(@(x) strsplit(A{x},' '),idx_1st','UniformOutput',false);
        El_wordpos=cellfun(@(x) 1:size(x,2),El_string,'UniformOutput',false);
        El_string=[El_string{:}];
        El_wordpos=[El_wordpos{:}];
        elec_stim_mat=[elec_stim_mat;El_cut_demeaned_ave];      
    end
     % make word 2 vec
     S=cellfun(@(x) strrep(x,"'S",""),El_string);
     S=cellfun(@(x) strrep(x,"TASHA","SASHA"),S);
     S=cellfun(@lower,S,'Uniformoutput',false);
     A=cellfun(@(x) double(word2vec(emb,x)),S,'UniformOutput',false);
        
    eval(strcat('sub_rdm_dat.sess_',stim_type{k},'_elec_stim_mat=elec_stim_mat;')) ;
    eval(strcat('sub_rdm_dat.sess_',stim_type{k},'_elec_stim_str=El_string;')) ;
    eval(strcat('sub_rdm_dat.sess_',stim_type{k},'_elec_stim_wordpos=El_wordpos;')) ;
    eval(strcat('sub_rdm_dat.sess_',stim_type{k},'_elec_stim_str_w2v=A;')) ;
end
%% 
conditions={'sentence','wordlist','jabberwocky','nonwords'};
all_RDMs=struct;
for i=1:length(conditions )
    A=fieldnames(sub_rdm_dat);
    idx=strcmp(A,strcat('sess_',conditions{i},'_elec_stim_mat'));
    B_cond=transpose(sub_rdm_dat.(A{idx}));
    RDM=squareform(pdist(B_cond,'correlation'));
    % 
    idx=strcmp(A,strcat('sess_',conditions{i},'_elec_stim_str_w2v'));
    B_w2v=transpose(sub_rdm_dat.(A{idx}));
    B_w2v=cell2mat(B_w2v);
    RDM_w2v=squareform(pdist(B_w2v,'correlation'));
    % 
    idx=strcmp(A,strcat('sess_',conditions{i},'_elec_stim_str'));
    STR=sub_rdm_dat.(A{idx});
    try
        Y= pdist(B_w2v,'correlation');
        Z = linkage(Y,'average');
        LO=optimalleaforder(Z,Y);
        RDM_w2v_order = squareform(pdist(B_w2v(LO,:),'cosine'));
        RDM_order=squareform(pdist(B_cond(LO,:),'correlation'));
        STR_order=STR(LO);
    catch err 
        RDM_order=nan*RDM;
        RDM_w2v_order=nan*RDM_w2v;
        STR_order=cellfun(@(x) nan,STR,'UniformOutput',false);
    end 
    % 
    all_RDMs(i).cond=conditions{i};
    all_RDMs(i).cond_rdm=RDM;
    all_RDMs(i).cond_str=STR;
    all_RDMs(i).cond_w2v_rdm=RDM_w2v;
    all_RDMs(i).cond_rdm_order=RDM_order;
    all_RDMs(i).cond_w2v_rdm_order=RDM_w2v_order;
    all_RDMs(i).cond_w2v_str_order=STR_order;
    fprintf(['added ',conditions{i},' rdm \n'])
end 

%% 

X=[all_RDMs(:).cond_rdm];
ax_lim=[min(X(:)),max(X(:))];
figure;
set(gcf,'position',[1000 152 1384 1186])
colors = cbrewer('div', 'RdBu', 128);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);

for i=1:size(all_RDMs,2)
    subplot(2,2,i)
    imagesc(all_RDMs(i).cond_rdm,[0,2]);
    daspect([1,1,1])
    title(all_RDMs(i).cond);
    set(gca, 'ydir', 'reverse','box','off');
    handles = colorbar;
    handles.TickDirection = 'out';
    handles.Box = 'off';
    drawnow;
end

if ~exist(analysis_path)
    mkdir(analysis_path)
end
print(gcf,'-fillpage', '-dpdf', strcat(analysis_path,'/','ECoG_SWJN_RDM.pdf'));

%% do the RDM based on sorted sentences and wordlist 
%sentence
conditions={'sentence','wordlist'};
X=cellfun(@(x) strcmp({all_RDMs(:).cond},x),conditions,'UniformOutput',false);
w2v_RDM=cellfun(@(x) all_RDMs(x).cond_rdm_order,X,'UniformOutput',false);
w2v_order=cellfun(@(x) all_RDMs(x).cond_w2v_rdm_order,X,'UniformOutput',false);
figure;
set(gcf,'position',[1000 152 1384 1186])
colors = cbrewer('div', 'RdBu', 128);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);
for i=1:size(w2v_RDM,2)
    subplot(2,2,2*i-1)
    imagesc(w2v_RDM{i},[0,2]);
    daspect([1,1,1])
    title({conditions{i},'word2vec order'});
    set(gca, 'ydir', 'reverse','box','off');
    handles = colorbar;
    handles.TickDirection = 'out';
    handles.Box = 'off';
    drawnow;
    subplot(2,2,2*i)
    imagesc(w2v_order{i},[0,2]);
    daspect([1,1,1])
    title({conditions{i},'word2vec strings ordered'});
    set(gca, 'ydir', 'reverse','box','off');
    handles = colorbar;
    handles.TickDirection = 'out';
    handles.Box = 'off';
    drawnow;
end
print(gcf,'-fillpage', '-dpdf', strcat(analysis_path,'/','ECoG_RDM_w2v_sorted.pdf'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compare against network 
%% get sentence and wordlist for model 
% wordlist 
model_stim=struct;
d_w=dir([data_path,'/models_output/*words.txt']);
i=1;
fid=fopen(strcat(d_w(i).folder,'/',d_w(i).name));
tline = fgetl(fid);
w_cell = cell(0,1);
while ischar(tline)
    w_cell{end+1,1} = strtrim(tline);
    tline = fgetl(fid);
end
fclose(fid);
model_stim(1).cond='W';
model_stim(1).stim=w_cell;
% sentences 
d_s=dir([data_path,'/models_output/*sentences.txt']);
i=1;
fid=fopen(strcat(d_s(i).folder,'/',d_s(i).name));
tline = fgetl(fid);
s_cell = cell(0,1);
while ischar(tline)
    s_cell{end+1,1} = strtrim(tline);
    tline = fgetl(fid);
end
fclose(fid);
model_stim(2).cond='S';
model_stim(2).stim=s_cell;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compare against network 
% load the pickle data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% bert base 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d_model= dir([data_path,'/models_output/*_model.mat']);
fprintf(' %d .pkl files were found \n', length(d_model));
models={'bert','gpt2','roberta','transfo'};
C={'S','W'};
conditions={'sentence','wordlist'};
all_model_dat=struct;
for i=1:length(models)
    Y=contains({d_model(:).name},strcat('_',models{i}));
    for k=1:length(C)
        X=contains({d_model(:).name},C{k});
        Z=contains({model_stim(:).cond},C{k});
        idx=X & Y;
        if i==4 % transfo is different 
            M=load([data_path,'/models_output/',d_model(idx).name]);
            M=M.data;
            M1=mat2cell(M,ones(1,size(M,1)),size(M,2),size(M,3));
            M=cellfun(@squeeze,M1,'UniformOutput',false);
        else 
            M=load([data_path,'/models_output/',d_model(idx).name]);
            M=transpose(M.data);
        end 
        all_model_dat(i).model=models{i};
        all_model_dat(i).(strcat(conditions{k},'_output'))=M;
        all_model_dat(i).(strcat(conditions{k},'_output_len'))=cellfun(@(x) size(x,1),M);
        all_model_dat(i).(strcat(conditions{k},'_stim'))=model_stim(Z).stim;
    end
end
%% find shared stim between model and data
stim_type={'sentence','wordlist'};
for k=1:size(stim_type,2)
    fprintf(['comparing ',stim_type{k},' to model stim\n'])
    eval(strcat('A=all_sub_dat(1).sess_',stim_type{k},'_stim_shared;')) ;
    for kk=1:size(all_model_dat,2)
        eval(strcat('B=all_model_dat(',num2str(kk),').',stim_type{k},'_stim;')) ;
        eval(strcat('M=all_model_dat(',num2str(kk),').model;')) ;
        A1=intersect(lower(A),lower(B));
        share_in_model=cell2mat(cellfun(@(x) find(contains(lower(B),x)),A1,'UniformOutput',false));
        eval(strcat('C=all_model_dat(',num2str(kk),').',stim_type{k},'_output_len;')) ;
        mdl_token_len=C(share_in_model);
        A2=B(share_in_model(mdl_token_len==8));
    for p=1:size(all_sub_dat,1)
        eval(strcat('all_sub_dat(p).sess_',stim_type{k},'_',M,'_stim_shared=A2;')) ;
    end
    eval(strcat('all_model_dat(',num2str(kk),').sess_',stim_type{k},'_stim_shared=A2;')) ;
    end 
end
%% contruct RDM between models and data 
%%
mdl_sub_rdm_dat=struct;
find_first_index= @(input_cell,val) find(cell2mat(cellfun(@(x) ~isempty(x), input_cell,'UniformOutput',false)),val,'first');
stim_type={'sentence','wordlist'};
model_types={all_model_dat(:).model};
for k=1:size(stim_type,2)
    fprintf(['adding ',stim_type{k},' electrode-stim mat\n'])
    
    for q=1:length(model_types)
        elec_stim_mat=[];
        elec_stim_str=[];
        elec_stim_wordpos=[];
    for kk=1:size(all_sub_dat,1)
        eval(strcat('A=all_sub_dat(',num2str(kk),').sess_',stim_type{k},'_stim;'));% get_stim
        eval(strcat('B=all_sub_dat(',num2str(kk),').sess_',stim_type{k},'_',model_types{q},'_stim_shared;'));% get_stim_shared with model
        idx=cellfun(@(x) (regexpi(A,x)),B,'UniformOutput',false);
        idx_1st=cellfun(@(x) find_first_index(x,1), idx, 'UniformOutput',false);
        if length(idx_1st)~=length(B)
            error('mismatch between index of data and model ') 
        end 
        %idx=cellfun(@(x) find_index(x), idx, 'UniformOutput',false);
        % get hilbert
        eval(strcat('Hilb=all_sub_dat(',num2str(kk),').sess_',stim_type{k},'_hilb_tensor;'));% get_stim_shared
        eval(strcat('Stim_L=all_sub_dat(',num2str(kk),').sess_stim_length;'));% get_word_length
        El=cellfun(@(x) double(squeeze(Hilb(:,x,:))),idx_1st','UniformOutput',false);
        El_demeaned=cellfun(@(x) transpose(bsxfun(@minus, x', mean(x'))),El,'UniformOutput',false);
        El_cut_demeaned=[];
        El_cut=[];
        a=size(El_demeaned{1},1);
        for p=1:size(El_demeaned,2)
            El_cut=[El_cut,mat2cell(El{p},a,Stim_L*ones(1,8))];
            El_cut_demeaned=[El_cut_demeaned,mat2cell(El_demeaned{p},a,Stim_L*ones(1,8))];
        end
        El_cut_ave=cell2mat(cellfun(@(x) nanmean(x,2),El_cut,'UniformOutput',false));
        El_cut_demeaned_ave=cell2mat(cellfun(@(x) nanmean(x,2),El_cut_demeaned,'UniformOutput',false));
        El_string=cellfun(@(x) strsplit(A{x},' '),idx_1st','UniformOutput',false);
        El_wordpos=cellfun(@(x) 1:size(x,2),El_string,'UniformOutput',false);
        El_string=[El_string{:}];
        El_wordpos=[El_wordpos{:}];
        elec_stim_mat=[elec_stim_mat;El_cut_demeaned_ave];      
        
    end
    % constrct mdl_stim_mat
    eval(strcat('mdl_sub_rdm_dat.sess_',stim_type{k},'_',model_types{q},'_elec_stim_mat=elec_stim_mat;')) ;
    eval(strcat('mdl_sub_rdm_dat.sess_',stim_type{k},'_',model_types{q},'_elec_stim_str=El_string;')) ;
    % do above for the model now 
    eval(strcat('mdl_type=all_model_dat(',num2str(q),').model;'));% get_stim
    if ~strcmp(mdl_type,model_types{q})
        error ('mismatch')
    end 
    eval(strcat('A_mdl=all_model_dat(',num2str(q),').',stim_type{k},'_stim;'));% get_stim
    eval(strcat('B_mdl=all_model_dat(',num2str(q),').sess_',stim_type{k},'_stim_shared;'));% get_stim_shared with model
    eval(strcat('Output_mdl=all_model_dat(',num2str(q),').',stim_type{k},'_output;'));% get_stim_shared with model
    idx=cellfun(@(x) (regexpi(A_mdl,x)),B_mdl,'UniformOutput',false);
        idx_1st=cellfun(@(x) find_first_index(x,1), idx, 'UniformOutput',false);
        if length(idx_1st)~=length(B_mdl)
            error('mismatch between index of data and model ') 
        end
        
        Mdl_El=cellfun(@double,Output_mdl([idx_1st{:}]),'UniformOutput',false);
        Mdl_El_demeaned=cellfun(@(x) bsxfun(@minus, x, mean(x)),Mdl_El,'UniformOutput',false);    
        mdl_El_cut_demeaned=[];
        for p=1:size(Mdl_El_demeaned,1)
            mdl_El_cut_demeaned=[mdl_El_cut_demeaned,Mdl_El_demeaned{p}'];
        end
        eval(strcat('mdl_sub_rdm_dat.mdl_',stim_type{k},'_',model_types{q},'_elec_stim_mat=mdl_El_cut_demeaned;')) ;
        
    end 
     % make word 2 vec

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% construct RDMS
conditions={'sentence','wordlist'};
model_types={all_model_dat(:).model};
all_mdl_sub_RDMs=struct;
A=fieldnames(mdl_sub_rdm_dat);
for i=1:length(conditions)
    for k=1:length(model_types)
    
    idx_sub=strcmp(A,strcat('sess_',conditions{i},'_',model_types{k},'_elec_stim_mat'));
    idx_mdl=strcmp(A,strcat('mdl_',conditions{i},'_',model_types{k},'_elec_stim_mat'));
    idx_mdl_str=strcmp(A,strcat('sess_',conditions{i},'_',model_types{k},'_elec_stim_str'));
    % 
    B_sub=transpose(mdl_sub_rdm_dat.(A{idx_sub}));
    RDM_sub=squareform(pdist(B_sub,'correlation'));
    %
    B_mdl=transpose(mdl_sub_rdm_dat.(A{idx_mdl}));
    RDM_mdl=squareform(pdist(B_mdl,'correlation'));
    %
    
    %
    flatten_RDM_sub=flatten_RDM(RDM_sub,1);
    flatten_RDM_mdl=flatten_RDM(RDM_mdl,1);
    [R_sub_mdl,p_sub_mdl] = corrcoef(flatten_RDM_sub,flatten_RDM_mdl);

    % 
    all_mdl_sub_RDMs(i).cond=conditions{i};
    all_mdl_sub_RDMs(i).(strcat(model_types{k},'_mdl_rdm'))=RDM_mdl;
    all_mdl_sub_RDMs(i).(strcat(model_types{k},'_sess_rdm'))=RDM_sub;
    all_mdl_sub_RDMs(i).(strcat(model_types{k},'_sess_corr'))=R_sub_mdl(1,2);

    fprintf(['added ',conditions{i},'-',model_types{k},' rdm \n'])
    end 
    w2v_RDM=all_RDMs(i).cond_w2v_rdm;
    flatten_RDM_w2v=flatten_RDM(w2v_RDM,1);
    [R_sub_w2v,p_sub_mdl] = corrcoef(flatten_RDM_sub,flatten_RDM_w2v);
    all_mdl_sub_RDMs(i).(strcat('w2v','_sess_corr'))=R_sub_w2v(1,2);
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

%%
%
i=1;

NN_sentences=S_NN_cell(H);
sentence_length=cell2mat(cellfun(@(x) size(x,1),NN_sentences,'UniformOutput',false));
H_intersect=H(sentence_length==8);
NN_sentences_new=S_NN_cell(H_intersect);
NN_sentences_new=cellfun(@double,NN_sentences_new,'UniformOutput',false);
bert_sentnce_mat=transpose(cell2mat(NN_sentences_new));
bert_sentence_RDM=squareform(pdist(bert_sentnce_mat','correlation'));

look_out_window=135; %200 ms
offset=0; % 100 ms
word_range=[1:8];
electrode_words_in_sentence_mat=[];
all_shared_sentences={};
for k=1:length(subject_ids)
    num_electrodes=size(all_subjects_data{k,2},1);
    electrode_sentence_time_tensor=all_subjects_data{k,2};
    word_length=all_subjects_data{k,7};
    subject_sentences=all_subjects_data{k,4};
    sentence_locations=cellfun(@(x) find(all_subjects_data{k,5}==x),mat2cell(H_intersect,ones(1,size(H_intersect,1)),1),'UniformOutput',false);
    shared_sentences=cellfun(@(x) subject_sentences{x(1)},sentence_locations,'UniformOutput',false);
    all_shared_sentences=[all_shared_sentences,shared_sentences];
    for i=1:num_electrodes
        a=[];b=[];
        electrode_response=squeeze(electrode_sentence_time_tensor(i,:,:));  
        electrode_sentence_time_cell=mat2cell(electrode_response,ones(1,size(electrode_response,1)),word_length*ones(1,8));
        electrode_shared_sentence_time_cell=cellfun(@(x) electrode_sentence_time_cell(x,:),sentence_locations,'UniformOutput',false);
        a=cellfun(@(x) cellfun(@(y) double(nanmean(y(offset+[1:look_out_window]))),x),electrode_shared_sentence_time_cell,'UniformOutput',false);
        b=cell2mat(cellfun(@(x) nanmean(x,1),a,'UniformOutput',false));
        b=transpose(b);
        electrode_words_in_sentence_mat=[electrode_words_in_sentence_mat;reshape(b,1,[])];
    end 
end 
sentence_vs_bert_mat=electrode_words_in_sentence_mat;
sentence_vs_bert_RDM=squareform(pdist(sentence_vs_bert_mat','correlation'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GPT2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=2 ;
S_gpt2_cell=load([data_path,'/models_output/',d_model(i).name]);
S_NN_cell=S_gpt2_cell.data;
S_NN_cell=transpose(S_NN_cell);
NN_sentences=S_NN_cell(H);
sentence_length=cell2mat(cellfun(@(x) size(x,1),NN_sentences,'UniformOutput',false));
H_intersect=H(sentence_length==8);
NN_sentences_new=S_NN_cell(H_intersect);
NN_sentences_new=cellfun(@double,NN_sentences_new,'UniformOutput',false);
gpt2_sentnce_mat=transpose(cell2mat(NN_sentences_new));
gpt2_sentence_RDM=squareform(pdist(gpt2_sentnce_mat','correlation'));

look_out_window=135; %200 ms
offset=0; % 100 ms
word_range=[1:8];
electrode_words_in_sentence_mat=[];
all_shared_sentences={};
for k=1:length(subject_ids)
    num_electrodes=size(all_subjects_data{k,2},1);
    electrode_sentence_time_tensor=all_subjects_data{k,2};
    word_length=all_subjects_data{k,7};
    subject_sentences=all_subjects_data{k,4};
    sentence_locations=cellfun(@(x) find(all_subjects_data{k,5}==x),mat2cell(H_intersect,ones(1,size(H_intersect,1)),1),'UniformOutput',false);
    shared_sentences=cellfun(@(x) subject_sentences{x(1)},sentence_locations,'UniformOutput',false);
    all_shared_sentences=[all_shared_sentences,shared_sentences];
    for i=1:num_electrodes
        a=[];b=[];
        electrode_response=squeeze(electrode_sentence_time_tensor(i,:,:));  
        electrode_sentence_time_cell=mat2cell(electrode_response,ones(1,size(electrode_response,1)),word_length*ones(1,8));
        electrode_shared_sentence_time_cell=cellfun(@(x) electrode_sentence_time_cell(x,:),sentence_locations,'UniformOutput',false);
        a=cellfun(@(x) cellfun(@(y) double(nanmean(y(offset+[1:look_out_window]))),x),electrode_shared_sentence_time_cell,'UniformOutput',false);
        b=cell2mat(cellfun(@(x) nanmean(x,1),a,'UniformOutput',false));
        b=transpose(b);
        electrode_words_in_sentence_mat=[electrode_words_in_sentence_mat;reshape(b,1,[])];
    end 
end 
sentence_vs_gpt2_mat=electrode_words_in_sentence_mat;
sentence_vs_gpt2_RDM=squareform(pdist(sentence_vs_gpt2_mat','correlation'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Roberta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=3;
NN_sentences_new=[];
S_roberta_cell=load([data_path,'/models_output/',d_model(i).name]);
S_NN_cell=transpose(S_roberta_cell.data);
NN_sentences=S_NN_cell(H);
sentence_length=cell2mat(cellfun(@(x) size(x,1),NN_sentences,'UniformOutput',false));
H_intersect=H(sentence_length==8);
NN_sentences_new=S_NN_cell(H_intersect);
NN_sentences_new=cellfun(@double,NN_sentences_new,'UniformOutput',false);
roberta_sentnce_mat=transpose(cell2mat(NN_sentences_new));
roberta_sentence_RDM=squareform(pdist(roberta_sentnce_mat','correlation'));

look_out_window=135; %200 ms
offset=0; % 100 ms
word_range=[1:8];
electrode_words_in_sentence_mat=[];
all_shared_sentences={};
for k=1:length(subject_ids)
    num_electrodes=size(all_subjects_data{k,2},1);
    electrode_sentence_time_tensor=all_subjects_data{k,2};
    word_length=all_subjects_data{k,7};
    subject_sentences=all_subjects_data{k,4};
    sentence_locations=cellfun(@(x) find(all_subjects_data{k,5}==x),mat2cell(H_intersect,ones(1,size(H_intersect,1)),1),'UniformOutput',false);
    shared_sentences=cellfun(@(x) subject_sentences{x(1)},sentence_locations,'UniformOutput',false);
    all_shared_sentences=[all_shared_sentences,shared_sentences];
    for i=1:num_electrodes
        a=[];b=[];
        electrode_response=squeeze(electrode_sentence_time_tensor(i,:,:));  
        electrode_sentence_time_cell=mat2cell(electrode_response,ones(1,size(electrode_response,1)),word_length*ones(1,8));
        electrode_shared_sentence_time_cell=cellfun(@(x) electrode_sentence_time_cell(x,:),sentence_locations,'UniformOutput',false);
        a=cellfun(@(x) cellfun(@(y) double(nanmean(y(offset+[1:look_out_window]))),x),electrode_shared_sentence_time_cell,'UniformOutput',false);
        b=cell2mat(cellfun(@(x) nanmean(x,1),a,'UniformOutput',false));
        b=transpose(b);
        electrode_words_in_sentence_mat=[electrode_words_in_sentence_mat;reshape(b,1,[])];
    end 
end 
sentence_vs_roberta_mat=electrode_words_in_sentence_mat;
sentence_vs_roberta_RDM=squareform(pdist(sentence_vs_roberta_mat','correlation'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% transfo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=4 ;
NN_sentences_new=[];
S_transfo_cell=load([data_path,'/models_output/',d_model(i).name]);
S_NN_tensor=S_transfo_cell.data;
NN_sentences=S_NN_tensor(H,:,:);
NN_sentences=double(NN_sentences);
NN_sentences_new=permute(NN_sentences,[2,1,3]);
NN_sentences_mat=reshape(NN_sentences_new,size(NN_sentences_new,2)*size(NN_sentences_new,1),[]);
transfo_sentence_RDM=squareform(pdist(NN_sentences_mat,'correlation'));
sentence_vs_transfo_RDM=sentence_RDM;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% wordlist

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BERT
i=5 ;
W_bert_cell=load([data_path,'/models_output/',d_model(i).name]);
W_NN_cell=transpose(W_bert_cell.data);
NN_words=W_NN_cell(H2);
word_length=cell2mat(cellfun(@(x) size(x,1),NN_words,'UniformOutput',false));
H2_intersect=H2(word_length==8);
NN_words_new=W_NN_cell(H2_intersect);
NN_words_new=cellfun(@double,NN_words_new,'UniformOutput',false);
bert_words_mat=transpose(cell2mat(NN_words_new));
bert_words_RDM=squareform(pdist(bert_words_mat','correlation'));

electrode_words_in_wordlist_mat=[];
for k=1:length(subject_ids)
    num_electrodes=size(all_subjects_data{k,2},1);
    electrode_wordlist_time_tensor=all_subjects_data{k,3};
    word_length=all_subjects_data{k,7};
    subject_words=all_subjects_data{k,6};
    wordlist_locations=cellfun(@(x) find(all_subjects_data{k,9}==x),mat2cell(H2_intersect,ones(1,size(H2_intersect,1)),1),'UniformOutput',false);
    shared_wordlist=cellfun(@(x) subject_words{x(1)},wordlist_locations,'UniformOutput',false);
    a=[];b=[];
    for i=1:num_electrodes
        electrode_response=squeeze(electrode_wordlist_time_tensor(i,:,:));
        electrode_wordlist_time_cell=mat2cell(electrode_response,ones(1,size(electrode_response,1)),word_length*ones(1,8));
        electrode_shared_wordlist_time_cell=cellfun(@(x) electrode_wordlist_time_cell(x,:),wordlist_locations,'UniformOutput',false);
        a=cellfun(@(x) cellfun(@(y) double(nanmean(y(offset+[1:look_out_window]))),x),electrode_shared_wordlist_time_cell,'UniformOutput',false);
        b=cell2mat(cellfun(@(x) nanmean(x,1),a,'UniformOutput',false));
        b=transpose(b);
        electrode_words_in_wordlist_mat=[electrode_words_in_wordlist_mat;reshape(b,1,[])];
    end 
end 
wordlist_vs_bert_mat=electrode_words_in_wordlist_mat;
wordlist_vs_bert_RDM=squareform(pdist(wordlist_vs_bert_mat','correlation'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GPT2
i=6 ;
W_gpt2_cell=load([data_path,'/models_output/',d_model(i).name]);
W_NN_cell=transpose(W_gpt2_cell.data);
NN_words=W_NN_cell(H2);
word_length=cell2mat(cellfun(@(x) size(x,1),NN_words,'UniformOutput',false));
H2_intersect=H2(word_length==8);
NN_words_new=W_NN_cell(H2_intersect);
gpt2_words_mat=transpose(cell2mat(NN_words_new));
gpt2_words_RDM=squareform(pdist(gpt2_words_mat','correlation'));

electrode_words_in_wordlist_mat=[];
for k=1:length(subject_ids)
    num_electrodes=size(all_subjects_data{k,2},1);
    electrode_wordlist_time_tensor=all_subjects_data{k,3};
    word_length=all_subjects_data{k,7};
    subject_words=all_subjects_data{k,6};
    wordlist_locations=cellfun(@(x) find(all_subjects_data{k,9}==x),mat2cell(H2_intersect,ones(1,size(H2_intersect,1)),1),'UniformOutput',false);
    shared_wordlist=cellfun(@(x) subject_words{x(1)},wordlist_locations,'UniformOutput',false);
    a=[];b=[];
    for i=1:num_electrodes
        electrode_response=squeeze(electrode_wordlist_time_tensor(i,:,:));
        electrode_wordlist_time_cell=mat2cell(electrode_response,ones(1,size(electrode_response,1)),word_length*ones(1,8));
        electrode_shared_wordlist_time_cell=cellfun(@(x) electrode_wordlist_time_cell(x,:),wordlist_locations,'UniformOutput',false);
        a=cellfun(@(x) cellfun(@(y) double(nanmean(y(offset+[1:look_out_window]))),x),electrode_shared_wordlist_time_cell,'UniformOutput',false);
        b=cell2mat(cellfun(@(x) nanmean(x,1),a,'UniformOutput',false));
        b=transpose(b);
        electrode_words_in_wordlist_mat=[electrode_words_in_wordlist_mat;reshape(b,1,[])];
    end 
end 
wordlist_vs_gpt2_mat=electrode_words_in_wordlist_mat;
wordlist_vs_gpt2_RDM=squareform(pdist(wordlist_vs_gpt2_mat','correlation'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Roberta

i=7 ;
W_roberta_cell=load([data_path,'/models_output/',d_model(i).name]);
W_NN_cell=transpose(W_roberta_cell.data);
NN_words=W_NN_cell(H2);
word_length=cell2mat(cellfun(@(x) size(x,1),NN_words,'UniformOutput',false));
H2_intersect=H2(word_length==8);
NN_words_new=W_NN_cell(H2_intersect);
roberta_words_mat=transpose(cell2mat(NN_words_new));
roberta_words_RDM=squareform(pdist(roberta_words_mat','correlation'));

electrode_words_in_wordlist_mat=[];
for k=1:length(subject_ids)
    num_electrodes=size(all_subjects_data{k,2},1);
    electrode_wordlist_time_tensor=all_subjects_data{k,3};
    word_length=all_subjects_data{k,7};
    subject_words=all_subjects_data{k,6};
    wordlist_locations=cellfun(@(x) find(all_subjects_data{k,9}==x),mat2cell(H2_intersect,ones(1,size(H2_intersect,1)),1),'UniformOutput',false);
    shared_wordlist=cellfun(@(x) subject_words{x(1)},wordlist_locations,'UniformOutput',false);
    a=[];b=[];
    for i=1:num_electrodes
        electrode_response=squeeze(electrode_wordlist_time_tensor(i,:,:));
        electrode_wordlist_time_cell=mat2cell(electrode_response,ones(1,size(electrode_response,1)),word_length*ones(1,8));
        electrode_shared_wordlist_time_cell=cellfun(@(x) electrode_wordlist_time_cell(x,:),wordlist_locations,'UniformOutput',false);
        a=cellfun(@(x) cellfun(@(y) double(nanmean(y(offset+[1:look_out_window]))),x),electrode_shared_wordlist_time_cell,'UniformOutput',false);
        b=cell2mat(cellfun(@(x) nanmean(x,1),a,'UniformOutput',false));
        b=transpose(b);
        electrode_words_in_wordlist_mat=[electrode_words_in_wordlist_mat;reshape(b,1,[])];
    end 
end 
wordlist_vs_roberta_mat=electrode_words_in_wordlist_mat;
wordlist_vs_roberta_RDM=squareform(pdist(electrode_words_in_wordlist_mat','correlation'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% transfo 
i=8 ;
NN_words_mat=[];
W_transfo_cell=load([data_path,'/models_output/',d_model(i).name]);
W_NN_tensor=W_transfo_cell.data;
NN_words=W_NN_tensor(H2,:,:);
NN_words=double(NN_words);
NN_words_new=permute(NN_words,[2,1,3]);
NN_words_mat=reshape(NN_words_new,size(NN_words_new,2)*size(NN_words_new,1),[]);
transfo_words_RDM=squareform(pdist(NN_words_mat,'correlation'));

electrode_words_in_wordlist_mat=[];
for k=1:length(subject_ids)
    num_electrodes=size(all_subjects_data{k,2},1);
    electrode_wordlist_time_tensor=all_subjects_data{k,3};
    word_length=all_subjects_data{k,7};
    subject_words=all_subjects_data{k,6};
    wordlist_locations=cellfun(@(x) find(all_subjects_data{k,9}==x),mat2cell(H2,ones(1,size(H2,1)),1),'UniformOutput',false);
    shared_wordlist=cellfun(@(x) subject_words{x(1)},wordlist_locations,'UniformOutput',false);
    a=[];b=[];
    for i=1:num_electrodes
        electrode_response=squeeze(electrode_wordlist_time_tensor(i,:,:));
        electrode_wordlist_time_cell=mat2cell(electrode_response,ones(1,size(electrode_response,1)),word_length*ones(1,8));
        electrode_shared_wordlist_time_cell=cellfun(@(x) electrode_wordlist_time_cell(x,:),wordlist_locations,'UniformOutput',false);
        a=cellfun(@(x) cellfun(@(y) double(nanmean(y(offset+[1:look_out_window]))),x),electrode_shared_wordlist_time_cell,'UniformOutput',false);
        b=cell2mat(cellfun(@(x) nanmean(x,1),a,'UniformOutput',false));
        b=transpose(b);
        electrode_words_in_wordlist_mat=[electrode_words_in_wordlist_mat;reshape(b,1,[])];
    end 
end 
wordlist_vs_transfo_mat=electrode_words_in_wordlist_mat;
wordlist_vs_transfo_RDM=squareform(pdist(wordlist_vs_transfo_mat','correlation'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
figure;
set(gcf,'position',[808 19 752 1319])
subplot(4,2,1)
colors = cbrewer('div', 'RdBu', 128);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);
imagesc(bert_sentence_RDM);
daspect([1,1,1])
set(gca, 'ydir', 'reverse','box','off');
title(['bert-Sentences-',num2str(size(bert_sentnce_mat,2))]);
handles = colorbar;
handles.TickDirection = 'out';
handles.Box = 'off';
drawnow;
subplot(4,2,2)
colors = cbrewer('div', 'RdBu', 128);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);
imagesc(sentence_vs_bert_RDM);
daspect([1,1,1])
set(gca, 'ydir', 'reverse','box','off');
title(['ecog-Sentences-',num2str(size(sentence_vs_bert_mat,2))]);
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
title(['gpt2-Sentences-',num2str(size(gpt2_sentnce_mat,2))]);
handles = colorbar;
handles.TickDirection = 'out';
handles.Box = 'off';
drawnow;
subplot(4,2,4)
colors = cbrewer('div', 'RdBu', 128);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);
imagesc(sentence_vs_gpt2_RDM);
daspect([1,1,1])
set(gca, 'ydir', 'reverse','box','off');
title(['ecog-Sentences-',num2str(size(sentence_vs_gpt2_mat,2))]);
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
title(['roberta-Sentences-',num2str(size(roberta_sentnce_mat,2))]);
handles = colorbar;
handles.TickDirection = 'out';
handles.Box = 'off';
drawnow;
subplot(4,2,6)
colors = cbrewer('div', 'RdBu', 128);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);
imagesc(sentence_vs_roberta_RDM);
daspect([1,1,1])
set(gca, 'ydir', 'reverse','box','off');
title(['ecog-Sentences-',num2str(size(sentence_vs_roberta_mat,2))]);
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
imagesc(sentence_vs_transfo_RDM);
daspect([1,1,1])
set(gca, 'ydir', 'reverse','box','off');
title('ecog-sentence');
handles = colorbar;
handles.TickDirection = 'out';
handles.Box = 'off';
drawnow;
print(gcf,'-fillpage', '-dpdf', strcat(analysis_path,'/','NN_vs_ECOG_RDM','.pdf'));
%% 
figure;
set(gcf,'position',[808 19 752 1319])
subplot(4,2,1)
colors = cbrewer('div', 'RdBu', 128);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);
imagesc(bert_sentence_RDM);
daspect([1,1,1])
set(gca, 'ydir', 'reverse','box','off');
title(['bert-Sentences-',num2str(size(bert_sentence_RDM,2))]);
handles = colorbar;
handles.TickDirection = 'out';
handles.Box = 'off';
drawnow;
subplot(4,2,2)
colors = cbrewer('div', 'RdBu', 128);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);
imagesc(bert_words_RDM);
daspect([1,1,1])
set(gca, 'ydir', 'reverse','box','off');
title(['bert-words-',num2str(size(bert_words_RDM,2))]);
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
title(['gpt2-Sentences-',num2str(size(gpt2_sentnce_mat,2))]);
handles = colorbar;
handles.TickDirection = 'out';
handles.Box = 'off';
drawnow;
subplot(4,2,4)
colors = cbrewer('div', 'RdBu', 128);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);
imagesc(gpt2_words_RDM);
daspect([1,1,1])
set(gca, 'ydir', 'reverse','box','off');
title(['gpt2-words-',num2str(size(gpt2_words_RDM,2))]);
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
title(['roberta-Sentences-',num2str(size(roberta_sentence_RDM,2))]);
handles = colorbar;
handles.TickDirection = 'out';
handles.Box = 'off';
drawnow;
subplot(4,2,6)
colors = cbrewer('div', 'RdBu', 128);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);
imagesc(roberta_words_RDM);
daspect([1,1,1])
set(gca, 'ydir', 'reverse','box','off');
title(['roberta-words-',num2str(size(roberta_words_RDM,2))]);
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
title(['transfo-sentence-',num2str(size(transfo_sentence_RDM,2))]);
handles = colorbar;
handles.TickDirection = 'out';
handles.Box = 'off';
drawnow;
subplot(4,2,8)
colors = cbrewer('div', 'RdBu', 128);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);
imagesc(transfo_words_RDM);
daspect([1,1,1])
set(gca, 'ydir', 'reverse','box','off');
title(['transfo-words-',num2str(size(transfo_words_RDM,2))]);

handles = colorbar;
handles.TickDirection = 'out';
handles.Box = 'off';
drawnow;
print(gcf,'-fillpage', '-dpdf', strcat(analysis_path,'/','NN_sentence_vs_words_RDM','.pdf'));
%% do correlation between RDMS 
flatten_bert_sentence_RDM=flatten_RDM(bert_sentence_RDM,1);
flatten_sentence_vs_bert_RDM=flatten_RDM(sentence_vs_bert_RDM,1);
flatten_bert_words_RDM=flatten_RDM(bert_words_RDM,1);
flatten_words_vs_bert_RDM=flatten_RDM(wordlist_vs_bert_RDM,1);
% 
flatten_gpt2_sentence_RDM=flatten_RDM(gpt2_sentence_RDM,1);
flatten_sentence_vs_gpt2_RDM=flatten_RDM(sentence_vs_gpt2_RDM,1);
flatten_gpt2_word_RDM=flatten_RDM(gpt2_words_RDM,1);
flatten_words_vs_gpt2_RDM=flatten_RDM(wordlist_vs_gpt2_RDM,1);
%
flatten_roberta_sentence_RDM=flatten_RDM(roberta_sentence_RDM,1);
flatten_sentence_vs_roberta_RDM=flatten_RDM(sentence_vs_roberta_RDM,1);
flatten_roberta_word_RDM=flatten_RDM(roberta_words_RDM,1);
flatten_words_vs_roberta_RDM=flatten_RDM(wordlist_vs_roberta_RDM,1);
% 
flatten_transfo_sentence_RDM=flatten_RDM(transfo_sentence_RDM,1);
flatten_sentence_vs_transfo_RDM=flatten_RDM(sentence_vs_transfo_RDM,1);
flatten_transfo_words_RDM=flatten_RDM(transfo_words_RDM,1);
flatten_words_vs_transfo_RDM=flatten_RDM(wordlist_vs_transfo_RDM,1);
% 
flatten_word2vec_sentence_RDM=flatten_RDM(sentence_vs_word2vec_RDM,1);
flatten_sentence_vs_w2vec_RDM=flatten_RDM(sentence_RDM,1);

flatten_word2vec_words_RDM=flatten_RDM(words_vs_word2vec_RDM,1);
flatten_words_vs_word2vec_RDM=flatten_RDM(wordlist_RDM,1);
%% compute a correlation between ecog data and networks 
[R_ecog_bert_sentence,p_s_bert] = corrcoef(flatten_sentence_vs_bert_RDM,flatten_bert_sentence_RDM);
[R_ecog_gpt2_sentence,p_s_gpt2] = corrcoef(flatten_sentence_vs_gpt2_RDM,flatten_gpt2_sentence_RDM);
[R_ecog_roberta_sentence,p_s_roberta] = corrcoef(flatten_sentence_vs_roberta_RDM,flatten_roberta_sentence_RDM);
[R_ecog_transfo_sentence,p_s_transfo] = corrcoef(flatten_sentence_vs_transfo_RDM,flatten_transfo_sentence_RDM);
[R_ecog_w2vec_sentence,p_s_w2vec] = corrcoef(flatten_sentence_vs_w2vec_RDM,flatten_word2vec_sentence_RDM);

% 
[R_ecog_bert_wordlist,p_w_bert] = corrcoef(flatten_words_vs_bert_RDM,flatten_bert_words_RDM);
[R_ecog_gpt2_wordlist,p_w_gpt2] = corrcoef(flatten_words_vs_gpt2_RDM,flatten_gpt2_word_RDM);
[R_ecog_roberta_wordlist,p_w_robert] = corrcoef(flatten_words_vs_roberta_RDM,flatten_roberta_word_RDM);
[R_ecog_transfo_wordlist,p_w_transfo] = corrcoef(flatten_words_vs_transfo_RDM,flatten_transfo_words_RDM);
[R_ecog_w2vec_wordlist,p_w_w2vec] = corrcoef(flatten_words_vs_word2vec_RDM,flatten_word2vec_words_RDM);

y=[R_ecog_bert_sentence(1,2),R_ecog_bert_wordlist(1,2);...
   R_ecog_gpt2_sentence(1,2),R_ecog_gpt2_wordlist(1,2);...
   R_ecog_roberta_sentence(1,2),R_ecog_roberta_wordlist(1,2);...
   R_ecog_transfo_sentence(1,2),R_ecog_transfo_wordlist(1,2);...
   R_ecog_w2vec_sentence(1,2),R_ecog_w2vec_wordlist(1,2)];

y_sig=[p_s_bert(1,2),p_w_bert(1,2);...
   p_s_gpt2(1,2),p_w_gpt2(1,2);...
   p_s_roberta(1,2),p_w_robert(1,2);...
   p_s_transfo(1,2),p_w_transfo(1,2);...
   p_s_w2vec(1,2),p_w_w2vec(1,2)]<0.05;


% 
figure;
bl=bar(y,'Facecolor','flat');
set(bl(1),'Linestyle','none','BarWidth',.8);
set(bl(2),'Linestyle','none','BarWidth',.8);
set(bl(1),'Displayname','Sentences');
set(bl(2),'Displayname','Wordlist');
set(gca,'xticklabel',{'bert','gpt2','roberta','transfo','word2vec'});
set(gca,'FontSize',14)
set(gca,'box','off')
legend('location','northeastoutside')
title('Correlation between ECoG RDM and Network RDMs');
ylabel('Correlation');

print(gcf,'-bestfit', '-dpdf', strcat(analysis_path,'/','ECoG_NN_RDM_correlation','.pdf'));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
shared_word_2_vec_sentences=cell2mat(cellfun(@(y) double(cell2mat(cellfun(@(x) word2vec(emb,lower(x)),strsplit(strrep(y,"'S",""),' '),'UniformOutput',false)')),shared_sentences,'UniformOutput',false));

shared_sentence_words=(cellfun(@(y) strsplit(strrep(y,"'S",""),' '),shared_sentences,'UniformOutput',false));
shared_sentence_words=[shared_sentence_words{:}]';
%
Y_sentences = pdist(shared_word_2_vec_sentences,'cosine');


shared_word_2_vec_wordlist=cell2mat(cellfun(@(y) double(cell2mat(cellfun(@(x) word2vec(emb,lower(x)),strsplit(strrep(strrep(y,"TASHA","SASHA"),"'S",""),' '),'UniformOutput',false)')),shared_wordlist,'UniformOutput',false));

shared_wordlist_words=(cellfun(@(y) strsplit(strrep(strrep(y,"TASHA","SASHA"),"'S",""),' '),shared_wordlist,'UniformOutput',false));
shared_wordlist_words=[shared_wordlist_words{:}]';
%
Y_wordlist = pdist(shared_word_2_vec_wordlist,'cosine');
save(strcat(analysis_path,'/','electrode_words_in_wordlist_matrix.mat'),'electrode_words_in_wordlist_mat');
wordlist_RDM=squareform(pdist(electrode_words_in_wordlist_mat','correlation'));
