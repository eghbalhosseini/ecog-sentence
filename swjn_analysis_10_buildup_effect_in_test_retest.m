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
analysis_path='/Users/eghbalhosseiniasl1/MyData/ECoG-sentence/analysis/swjn_analysis_10_buidlup_effect_in_test_retest/';
%
if 1
    fprintf('adding basic ecog tools to path \n');
    addpath('~/MyCodes/basic-ecog-tools/');
    addpath('~/MyCodes/ecog-sentence/');
    addpath(genpath('~/MyCodes/basic-ecog-tools/activeBrain'));
    addpath(genpath('~/MyCodes/basic-ecog-tools/ecog-filters'));
    addpath(genpath('~/MyCodes/basic-ecog-tools/mex'));
end

find_index= @(input_cell) find(cell2mat(cellfun(@(x) ~isempty(x), input_cell,'UniformOutput',false)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  STEP 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
session_sentence_examples={};
session_wordlist_examples={};
session_sentence_hilbert_band_envelope_tensor=[];
session_words_hilbert_band_envelope_tensor=[];
%
for i=1:length(d)
    subj=load(strcat(d(i).folder,'/',d(i).name));
    subj_id=fieldnames(subj);
    subj=subj.(subj_id{1});
    data=subj.data;
    info=subj.info;
    %language_electrode=info.language_responsive_electrodes_hilbert_odd;
    language_electrode=info.ramp_electrodes_hilbert_odd;
    language_electrode_num=find(language_electrode);
    % step 1: extract electrodes with siginificant language response
    sentence_trial_index=~cellfun(@isempty,cellfun(@ (x) strfind(x,'S'),info.word_type,'UniformOutput',false));
    %%%%%%%%%%%%%%% sentence
    hilbert_band_envelope=[];
    sentences=[data{sentence_trial_index}];%
    word_length=sentences(1).signal_range_downsample(1,2)-sentences(1).signal_range_downsample(1,1)+1;
    % creat a cell with wordposition(row)*time in trial(column) structure
    % hilbert mean (changlab)
    hilbert_band_envelope=cellfun(@(x) x(1:8),{sentences.signal_hilbert_zs_downsample_parsed},'UniformOutput',false);
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
    hilbert_band_envelope=cellfun(@(x) x(1:8),{words.signal_hilbert_zs_downsample_parsed},'UniformOutput',false);
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
    %
    fprintf('added %s from %s \n',d(i).name, strcat(d(i).folder,'/',d(i).name));
    clear data subj
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% extract location of sentence and wordlist repetitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sentence
session_sentence_examples_1=cellfun(@(x) convertCharsToStrings(x),session_sentence_examples,'UniformOutput',false);
session_sentence_examples_split=transpose(cellfun(@(x) strsplit(x,' '),session_sentence_examples,'UniformOutput',false));
unique_sentences=table2cell(unique(cell2table(session_sentence_examples_1)));
unique_sentences=cellfun(@(x) x{:},unique_sentences,'UniformOutput',false);
unique_sentence_location=cellfun(@(x) (regexpi(session_sentence_examples,x)),unique_sentences,'UniformOutput',false);
unique_sentence_location=cellfun(@(x) find_index(x), unique_sentence_location, 'UniformOutput',false);
% wordlist
session_wordlist_examples_1=cellfun(@(x) convertCharsToStrings(x),session_wordlist_examples,'UniformOutput',false);
session_wordlist_examples_split=transpose(cellfun(@(x) strsplit(x,' '),session_wordlist_examples,'UniformOutput',false));
unique_wordlist=table2cell(unique(cell2table(session_wordlist_examples_1)));
unique_wordlist=cellfun(@(x) x{:},unique_wordlist,'UniformOutput',false);
unique_wordlist_location=cellfun(@(x) (regexpi(session_wordlist_examples,x)),unique_wordlist,'UniformOutput',false);
unique_wordlist_location=cellfun(@(x) find_index(x), unique_wordlist_location, 'UniformOutput',false);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  get common sentence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
sentence_repetition_location=cell2mat(unique_sentence_location(cellfun(@length,unique_sentence_location)>1)');
first_repeat_sentences=sentence_repetition_location(1,:);
second_repeat_sentences=sentence_repetition_location(2,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot average for pre-post in word closing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
condition='mean';
%
% first sentence
idx=first_repeat_sentences;
A=double((session_sentence_hilbert_band_envelope_tensor(:,idx,:)));
B=reshape(A,size(A,1)*size(A,2),size(A,3));
electrode_sentence_mean=mat2cell(B,size(A,2)*ones(1,size(A,1)),size(A,3))
A=transpose(electrode_sentence_response);
electrode_sentence_response=transpose(bsxfun(@minus, A, mean(A)));
[a,b]=size(electrode_sentence_response);
electrode_sentence_response_cell=mat2cell(electrode_sentence_response,a,word_length*ones(1,8));
electrode_sentence_response_cell_1=[electrode_sentence_response_cell;session_sentence_examples_split{idx}];
all_first_sentence_instances_signal=[all_first_sentence_instances_signal,electrode_sentence_response_cell_1];
%second sentence
idx=second_repeat_index;
electrode_sentence_response=double(squeeze(session_sentence_hilbert_band_envelope_tensor(:,idx,:)));
A=transpose(electrode_sentence_response);
electrode_sentence_response=transpose(bsxfun(@minus, A, mean(A)));
[a,b]=size(electrode_sentence_response);
electrode_sentence_response_cell=mat2cell(electrode_sentence_response,a,word_length*ones(1,8));
electrode_sentence_response_cell_1=[electrode_sentence_response_cell;session_sentence_examples_split{idx}];
all_second_sentence_instances_signal=[all_second_sentence_instances_signal,electrode_sentence_response_cell_1];


idx=first_repeat_index;
electrode_sentence_response=double(squeeze(session_sentence_hilbert_band_envelope_tensor(:,idx,:)));
A=transpose(electrode_sentence_response);
electrode_sentence_response=transpose(bsxfun(@minus, A, mean(A)));
[a,b]=size(electrode_sentence_response);
electrode_sentence_response_cell=mat2cell(electrode_sentence_response,a,word_length*ones(1,8));
electrode_sentence_response_cell_mean=cellfun(@(x) mean(x,2),electrode_sentence_response_cell,'UniformOutput',false);
electrode_sentence_response_cell_1=[electrode_sentence_response_cell_mean;session_sentence_examples_split{idx}];
all_first_sentence_instances_signal=[all_first_sentence_instances_signal,electrode_sentence_response_cell_1];
%second sentence
idx=second_repeat_index;
electrode_sentence_response=double(squeeze(session_sentence_hilbert_band_envelope_tensor(:,idx,:)));
A=transpose(electrode_sentence_response);
electrode_sentence_response=transpose(bsxfun(@minus, A, mean(A)));
[a,b]=size(electrode_sentence_response);
electrode_sentence_response_cell=mat2cell(electrode_sentence_response,a,word_length*ones(1,8));
electrode_sentence_response_cell_mean=cellfun(@(x) mean(x,2),electrode_sentence_response_cell,'UniformOutput',false);
electrode_sentence_response_cell_1=[electrode_sentence_response_cell_mean;session_sentence_examples_split{idx}];
all_second_sentence_instances_signal=[all_second_sentence_instances_signal,electrode_sentence_response_cell_1];

%%
first_second_rep_corr=[];
first_first_rep_corr=[];
second_second_repetition_correlation=[];


for i=1:size(all_first_sentence_instances_signal,2)
    between_rep_xcorr={};
    within_rep_xcorr={};
    A=all_first_sentence_instances_signal{1,i};
    B=all_second_sentence_instances_signal(1,:);
    C=all_first_sentence_instances_signal(1,:);
    % between
    between_rep_xcorr=cellfun(@(x) corrcoef(x(:),A(:)),B,'UniformOutput',false);
    between_rep_xcorr=cellfun(@(x) x(1,2),between_rep_xcorr);
    % within
    within_rep_xcorr=cellfun(@(x) corrcoef(x(:),A(:)),C,'UniformOutput',false);
    within_rep_xcorr=cellfun(@(x) x(1,2),within_rep_xcorr);
    % add them to the matrix 
    first_second_rep_corr=[first_second_rep_corr;between_rep_xcorr];
    first_first_rep_corr=[first_first_rep_corr;within_rep_xcorr];
end
% 
A=first_first_rep_corr;
diag_mean=nanmean(diag(A));
A_off_diag=A+diag(nan*ones(1,size(A,1)));
off_diag_mean=nanmean(A_off_diag(:));


%
colors = cbrewer('div', 'RdBu', 128);
colors=flipud(colors);
colormap(colors);
subplot(2,1,1)
imagesc(A);
daspect([1,1,1]);
handles = colorbar;
handles.TickDirection = 'out';
handles.Box = 'off';
ylabel('first repetition');
xlabel('first repetition');
title({'Sentences : ',condition ,'; S>N+buildup',...
    ['diag mean: ',num2str(diag_mean)],['off diag mean: ',num2str(off_diag_mean)]})

%
A=first_second_rep_corr;
diag_mean=nanmean(diag(A));
A_off_diag=A+diag(nan*ones(1,size(A,1)));
off_diag_mean=nanmean(A_off_diag(:));
subplot(2,1,2);
imagesc(A);
daspect([1,1,1]);
handles = colorbar;
handles.TickDirection = 'out';
handles.Box = 'off';
ylabel('first repetition')
xlabel('second repetition')
title({'Sentences : ',condition ,'; S>N+buildup',...
    ['diag mean: ',num2str(diag_mean)],['off diag mean: ',num2str(off_diag_mean)]})

%
if ~exist(strcat(analysis_path,info.subject))
    mkdir(strcat(analysis_path,info.subject))
end

print(gcf, '-fillpage','-dpdf', strcat(analysis_path,info.subject,'/',info.subject,'_reliability_to_repetition_buildup_',condition,'.pdf'));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot average for pre-post in word closing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







