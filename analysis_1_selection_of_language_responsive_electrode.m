% selecting langauge responsive electrodes and add them to the subject
% info. 
% tested for subject 1 and validated with Terry's data . 
%% step 0: prepare the data 
clear all 
save_path='/Users/eghbalhosseiniasl1/MyData/ECoG-sentence/crunched/';
data_path='/Users/eghbalhosseiniasl1/MyData/ecog-sentence';
d= dir([data_path,'/**/AMC026*_crunched.mat']);
fprintf(' %d .mat files were found \n', length(d));
gamma_band_index=4;
num_of_permutation=1000; 
p_threshold=0.01;
%% 
sentence_gamma_ave_all=[];
gamma_ave=[];
gamma_condition=[];
nonwords_gamma_ave_all=[];
for k=1:2:length(d)
    fprintf('adding %s from %s \n',d(k).name, strcat(d(k).folder,'/',d(k).name));
    subj=load(strcat(d(k).folder,'/',d(k).name));
    subj_id=fieldnames(subj);
    %subj=subj.(subj_id{1});
    data=subj.data;
    info=subj.info;
    
% step 1: compute mean across the word positions in each trial 

    sentence_trial_index=~cellfun(@isempty,cellfun(@ (x) strfind(x,'S'),info.word_type,'UniformOutput',false));
    % sentence 
    sentences=[data{sentence_trial_index}];
    sentence_gamma_band_ave_envelope=cellfun(@(x) x(1:8,gamma_band_index),{sentences.signal_ave_envelope_downsample_parsed},'UniformOutput',false);
    sentence_gamma_band_ave_envelope=[sentence_gamma_band_ave_envelope{:,:}];
    % convert the gamma into a tensor of shape : channel*trial*words
    sentence_gamma_band_ave_envelope_tensor=cell2mat(permute(sentence_gamma_band_ave_envelope,[3,2,1]));
    %words=~cell2mat(cellfun(@isempty,cellfun(@(x) strfind(x,'word'),[sentences.stimuli_type],'UniformOutput',false),'UniformOutput',false));
    sentence_gamma_ave=nanmean(sentence_gamma_band_ave_envelope_tensor,3);
    %  nonword trials 
    nonword_trial_index=~cellfun(@isempty,cellfun(@ (x) strfind(x,'N'),info.word_type,'UniformOutput',false));
    nonwords=[data{nonword_trial_index}];
    nonword_gamma_band_ave_envelope=cellfun(@(x) x(1:8,gamma_band_index),{nonwords.signal_ave_envelope_downsample_parsed},'UniformOutput',false);
    nonword_gamma_band_ave_envelope=[nonword_gamma_band_ave_envelope{:,:}];
    % convert the gamma into a tensor of shape : channel*trial*words
    nonwords_gamma_band_ave_envelope_tensor=cell2mat(permute(nonword_gamma_band_ave_envelope,[3,2,1]));
    %words=~cell2mat(cellfun(@isempty,cellfun(@(x) strfind(x,'word'),[sentences.stimuli_type],'UniformOutput',false),'UniformOutput',false));
    nonwords_gamma_ave=nanmean(nonwords_gamma_band_ave_envelope_tensor,3);
    % add the new ones to the list 
    sentence_gamma_ave_all=[sentence_gamma_ave_all,sentence_gamma_ave];
    nonwords_gamma_ave_all=[nonwords_gamma_ave_all,nonwords_gamma_ave];
    gamma_ave=[gamma_ave,[sentence_gamma_ave,nonwords_gamma_ave]];
    gamma_condition=[gamma_condition,[sentence_gamma_ave*0+1,nonwords_gamma_ave*0-1]];
    fprintf('added %s from %s \n',d(k).name, strcat(d(k).folder,'/',d(k).name));
end
    
%% step 2: compute a correlation between trial means and vector of condition labels )
% sentences = 1, nonword-lists =-1 
gamma_ave=transpose(gamma_ave);
gamma_condition=transpose(gamma_condition);
gamma_ave=double(gamma_ave);
gamma_condition=double(gamma_condition);
[RHO,PVAL] = corr(gamma_ave,gamma_condition,'Type','Spearman');
rho_original=diag(RHO);
rho_original_repeat=repmat(rho_original,[1,1000]);
rho_positive=rho_original>0;
%% step 3: random permutation of conditions labels, repeat 1000 times 
rho_permuted=[];
for k=1:1000
     if ~mod(k,100)
    fprintf(' %d \n',k );
    end
    random_index=randperm(size(gamma_condition,1));
    [RHO_rand,~] = corr(gamma_ave,gamma_condition(random_index,:),'Type','Spearman');
    rho_permuted=[rho_permuted,diag(RHO_rand)];
end 

%% step 4 : compute fraction of correlations in step 3 that produce higher correlation that step 2. 
p_fraction=sum(rho_permuted>rho_original_repeat,2)./size(rho_permuted,2);
p_significant=p_fraction<p_threshold;
channel_significant=p_significant.*rho_positive;

%% step 5: add the language electrode to back to the data 
for k=1:length(d)
    fprintf('adding language electrodes to %s \n', strcat(d(k).folder,'/',d(k).name));
    subj=load(strcat(d(k).folder,'/',d(k).name));
    subj_id=fieldnames(subj);
    %subj=subj.(subj_id{1});
    dat=subj.data;
    info=subj.info;
    info.language_responsive_electrodes=channel_significant;
    subject_name=info.subject;
    session_name=info.session_name;
    eval(strcat(subject_name,'_',session_name,'.data=dat;')) ;
    eval(strcat(subject_name,'_',session_name,'.info=info;'));
    
    save(strcat(save_path,subject_name,'_',session_name,'_crunched.mat'),strcat(subject_name,'_',session_name),'-v7.3');
end 
