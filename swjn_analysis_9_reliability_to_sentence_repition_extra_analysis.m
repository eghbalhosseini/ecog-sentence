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
analysis_path='/Users/eghbalhosseiniasl1/MyData/ECoG-sentence/analysis/swjn_analysis_9_reliability_to_repetition/';
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
%%  STEP 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
find_index= @(input_cell) find(cell2mat(cellfun(@(x) ~isempty(x), input_cell,'UniformOutput',false)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
electrode_group='language';

session_sentence_examples={};
session_wordlist_examples={};
session_nonwords_examples={};
session_jabberwocky_examples={};
session_sentence_hilbert_band_envelope_tensor=[];
session_words_hilbert_band_envelope_tensor=[];
session_nonwords_hilbert_band_envelope_tensor=[];
session_jabberwocky_hilbert_band_envelope_tensor=[];
%
for i=1:length(d)
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
% nonwords 
session_nonwords_examples_1=cellfun(@(x) convertCharsToStrings(x),session_nonwords_examples,'UniformOutput',false);
session_nonwords_examples_split=transpose(cellfun(@(x) strsplit(x,' '),session_nonwords_examples,'UniformOutput',false));
unique_nonwords=table2cell(unique(cell2table(session_nonwords_examples_1)));
unique_nonwords=cellfun(@(x) x{:},unique_nonwords,'UniformOutput',false);
unique_nonwords_location=cellfun(@(x) (regexpi(session_nonwords_examples,x)),unique_nonwords,'UniformOutput',false);
unique_nonwords_location=cellfun(@(x) find_index(x), unique_nonwords_location, 'UniformOutput',false);
% jabberwocky
session_jabberwocky_examples_1=cellfun(@(x) convertCharsToStrings(x),session_jabberwocky_examples,'UniformOutput',false);
session_jabberwocky_examples_split=transpose(cellfun(@(x) strsplit(x,' '),session_jabberwocky_examples,'UniformOutput',false));
unique_jabberwocky=table2cell(unique(cell2table(session_jabberwocky_examples_1)));
unique_jabberwocky=cellfun(@(x) x{:},unique_jabberwocky,'UniformOutput',false);
unique_jabberwocky_location=cellfun(@(x) (regexpi(session_jabberwocky_examples,x)),unique_jabberwocky,'UniformOutput',false);
unique_jabberwocky_location=cellfun(@(x) find_index(x), unique_jabberwocky_location, 'UniformOutput',false);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  define condition to run analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
close all;
condition='sentence';
window_condition='mean';
sampling_coeff=450/135;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch condition
    case 'sentence'
        A=unique_sentence_location;
        B=cell2mat(A(cellfun(@length,A)>1)');
        first_repeat=B(1,:);
        second_repeat=B(2,:);
        session_hilbert_envelope=session_sentence_hilbert_band_envelope_tensor;
        session_examples_split=session_sentence_examples_split;
    case 'nonwords'
        A=unique_nonwords_location;
        B=cell2mat(A(cellfun(@length,A)>1)');
        first_repeat=B(1,:);
        second_repeat=B(2,:);
        session_hilbert_envelope=session_nonwords_hilbert_band_envelope_tensor;
        session_examples_split=session_nonwords_examples_split;
    case 'wordlist'
        A=unique_wordlist_location;
        B=cell2mat(A(cellfun(@length,A)>1)');
        first_repeat=B(1,:);
        second_repeat=B(2,:);
        session_hilbert_envelope=session_words_hilbert_band_envelope_tensor;
        session_examples_split=session_wordlist_examples_split;
    case 'jabberwocky'
        A=unique_jabberwocky_location;
        B=cell2mat(A(cellfun(@length,A)>1)');
        first_repeat=B(1,:);
        second_repeat=B(2,:);
        session_hilbert_envelope=session_jabberwocky_hilbert_band_envelope_tensor;
        session_examples_split=session_jabberwocky_examples_split;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get signals
% 
all_second_instances_signal=[];
all_first_instances_signal=[];
for i=1:length(first_repeat)
    first_repeat_index=first_repeat(i);
    second_repeat_index=second_repeat(i);    
    % first instance
    idx=first_repeat_index;
    electrode_response=double(squeeze(session_hilbert_envelope(:,idx,:)));
    A=transpose(electrode_response);
    electrode_response=transpose(bsxfun(@minus, A, mean(A)));
    [a,b]=size(electrode_response);
    electrode_response_cell=mat2cell(electrode_response,a,word_length*ones(1,8));
    if strcmp(window_condition,'full')
        electrode_response_cell=electrode_response_cell;
    elseif strcmp(window_condition,'mean')
        electrode_response_cell=cellfun(@(x) mean(x,2),electrode_response_cell,'UniformOutput',false);
    end
    electrode_response_cell_1=[electrode_response_cell;session_examples_split{idx}];
    all_first_instances_signal=[all_first_instances_signal,electrode_response_cell_1];
    %second instance
    idx=second_repeat_index;
    electrode_response=double(squeeze(session_hilbert_envelope(:,idx,:)));
    A=transpose(electrode_response);
    electrode_response=transpose(bsxfun(@minus, A, mean(A)));
    [a,b]=size(electrode_response);
    electrode_response_cell=mat2cell(electrode_response,a,word_length*ones(1,8));
    if strcmp(window_condition,'full')
        electrode_response_cell=electrode_response_cell;
    elseif strcmp(window_condition,'mean')
        electrode_response_cell=cellfun(@(x) mean(x,2),electrode_response_cell,'UniformOutput',false);
    end
    electrode_response_cell_1=[electrode_response_cell;session_examples_split{idx}];
    all_second_instances_signal=[all_second_instances_signal,electrode_response_cell_1];
 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
first_instance=all_first_instances_signal;
second_instance=all_second_instances_signal;
first_second_rep_corr=[];
first_first_rep_corr=[];
second_second_rep_corr=[];
% 
for i=1:size(first_instance,2)
    between_rep_xcorr={};
    within_rep_xcorr={};
    A=first_instance{1,i};
    B=second_instance(1,:);
    C=first_instance(1,:);
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
figure;
set(gcf,'position',[1436,374,981,960])
colors = cbrewer('div', 'RdBu', 128);
colors=flipud(colors);colormap(colors);
subplot(2,1,1)
imagesc(A);
daspect([1,1,1]);
handles = colorbar;
handles.TickDirection = 'out';
handles.Box = 'off';
ylabel('first repetition');
xlabel('first repetition');
title({condition ,'; electrodes: ',electrode_group,...
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
title({condition ,'; electrodes: ',electrode_group,...
    ['diag mean: ',num2str(diag_mean)],['off diag mean: ',num2str(off_diag_mean)]})

%
if ~exist(strcat(analysis_path,info.subject))
    mkdir(strcat(analysis_path,info.subject))
end
print(gcf, '-fillpage','-dpdf', strcat(analysis_path,info.subject,'/',info.subject,'_reliability_to_repetition_',electrode_group,'_elec_',condition,'.pdf'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RDM comparison 

first_instance=all_first_instances_signal;
second_instance=all_second_instances_signal;

A=first_instance(1,:);
B=transpose(cell2mat(A));
first_RDM=squareform(pdist(B,'correlation'));
A=second_instance(1,:);
B=transpose(cell2mat(A));
second_RDM=squareform(pdist(B,'correlation'));
%
flatten_first_RDM=flatten_RDM(first_RDM,1);
flatten_second_RDM=flatten_RDM(second_RDM,1);
[R_corr,p_corr] = corrcoef(flatten_first_RDM,flatten_second_RDM);
first_second_rep_corr=R_corr(1,2);
% do permutation 
r = arrayfun(@(x) transpose(randperm(length(flatten_first_RDM))),[1:1000],'UniformOutput',false);
permuted_flatten_second=cellfun(@(x) flatten_second_RDM(x),r,'UniformOutput',false);
permuted_corr=cellfun(@(x) corrcoef(flatten_first_RDM,x),permuted_flatten_second,'UniformOutput',false);
permuted_corr=cellfun(@(x) x(1,2),permuted_corr);
figure;
set(gcf,'position',[1436,374,981,960])
colors = cbrewer('div', 'RdBu', 256);
colors=flipud(colors);
colormap(colors);
%
ax=axes('position',[.1,.5,.3,.3]);
img=imagesc(first_RDM,[0,2]);
daspect([1,1,1]);
handles = colorbar;
handles.TickDirection = 'out';
handles.Box = 'off';
ylabel('first repetition');
xlabel('first repetition');
title({strcat(condition,' : First RDM'), 'electrodes: ',electrode_group});
%
ax=axes('position',[.1,.1,.3,.3]);
imagesc(second_RDM,[0,2]);
daspect([1,1,1]);
handles = colorbar;
handles.TickDirection = 'out';
handles.Box = 'off';
ylabel('second repetition');
xlabel('second repetition');
title({strcat(condition,' : second RDM'), 'electrodes: ',electrode_group});
%
ax1=subplot('position',[.5,.1,.5,.6]);
b_jitter_permuted=beeswarm((permuted_corr*0+1)',permuted_corr','hex','none',.3);
bl=get(ax1,'Children');
ax1.YLimMode='auto';
bl.MarkerFaceColor=[.5,.5,.5]
bl.DisplayName=['Permutated RDM']
hold on
e=errorbar(1,mean(permuted_corr),std(permuted_corr));
set(e,'Capsize',4);set(e,'Linewidth',2);set(e,'marker','o');set(e,'Color','k')
hAnnotation=arrayfun(@(x) get(x,'Annotation'),e);
hLegendEntry = arrayfun(@(x) get(x,'LegendInformation'),hAnnotation);
arrayfun(@(x) set(x,'IconDisplayStyle','off'),hLegendEntry);
bl1=scatter(1,first_second_rep_corr,50,'filled','DisplayName','RDM Correlation');
bl1.MarkerFaceColor=[0,0,0]; bl1.MarkerFaceAlpha=1; bl1.MarkerEdgeColor=[0,0,0]; bl1.LineWidth=2;
ax1.XAxis.Visible = 'off';
ax1.YLabel.String='Correlation';
ax1.Title.String={'Correlation between 1st and 2nd RDM','vs 1st and permuted 2nd RDM'}
y_lim=ax1.YLim;
x_lim=ax1.XLim;
ax1.YLim=[-max(y_lim),max(y_lim)];
ax1.XLim=[.7,1.3];
legend('location','northeastoutside')
print(gcf, '-bestfit','-dpdf', strcat(analysis_path,info.subject,'/',info.subject,'_RDM_correlation_',electrode_group,'_elec_',condition,'.pdf'));
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute one shot within and between correlations 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=all_first_instances_signal(1,:);
B=all_second_instances_signal(1,:);
% first 
C=cellfun(@(x) cellfun(@(y) corrcoef(x,y),A','UniformOutput',false),A,'UniformOutput',false);
first_first_corr=cell2mat(cellfun(@(x) x(1,2), [C{:,:}],'UniformOutput',false));
% second 
C=cellfun(@(x) cellfun(@(y) corrcoef(x,y),B','UniformOutput',false),B,'UniformOutput',false);
second_second_corr=cell2mat(cellfun(@(x) x(1,2), [C{:,:}],'UniformOutput',false));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% do the analysis for all condition. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
conditions={'sentence','wordlist','jabberwocky','nonwords'};
window_condition='mean';
all_RDMs=struct;
for k=1:length(conditions)
    condition=conditions{k};
    switch condition
    case 'sentence'
        A=unique_sentence_location;
        B=cell2mat(A(cellfun(@length,A)>1)');
        first_repeat=B(1,:);
        second_repeat=B(2,:);
        session_hilbert_envelope=session_sentence_hilbert_band_envelope_tensor;
        session_examples_split=session_sentence_examples_split;
    case 'nonwords'
        A=unique_nonwords_location;
        B=cell2mat(A(cellfun(@length,A)>1)');
        first_repeat=B(1,:);
        second_repeat=B(2,:);
        session_hilbert_envelope=session_nonwords_hilbert_band_envelope_tensor;
        session_examples_split=session_nonwords_examples_split;
    case 'wordlist'
        A=unique_wordlist_location;
        B=cell2mat(A(cellfun(@length,A)>1)');
        first_repeat=B(1,:);
        second_repeat=B(2,:);
        session_hilbert_envelope=session_words_hilbert_band_envelope_tensor;
        session_examples_split=session_wordlist_examples_split;
    case 'jabberwocky'
        A=unique_jabberwocky_location;
        B=cell2mat(A(cellfun(@length,A)>1)');
        first_repeat=B(1,:);
        second_repeat=B(2,:);
        session_hilbert_envelope=session_jabberwocky_hilbert_band_envelope_tensor;
        session_examples_split=session_jabberwocky_examples_split;
    end 
    % 
    % get signals
% 
    all_second_instances_signal=[];
    all_first_instances_signal=[];
    for i=1:length(first_repeat)
    first_repeat_index=first_repeat(i);
    second_repeat_index=second_repeat(i);    
    % first instance
    idx=first_repeat_index;
    electrode_response=double(squeeze(session_hilbert_envelope(:,idx,:)));
    A=transpose(electrode_response);
    electrode_response=transpose(bsxfun(@minus, A, mean(A)));
    [a,b]=size(electrode_response);
    electrode_response_cell=mat2cell(electrode_response,a,word_length*ones(1,8));
    if strcmp(window_condition,'full')
        electrode_response_cell=electrode_response_cell;
    elseif strcmp(window_condition,'mean')
        electrode_response_cell=cellfun(@(x) mean(x,2),electrode_response_cell,'UniformOutput',false);
    end
    electrode_response_cell_1=[electrode_response_cell;session_examples_split{idx}];
    all_first_instances_signal=[all_first_instances_signal,electrode_response_cell_1];
    %second instance
    idx=second_repeat_index;
    electrode_response=double(squeeze(session_hilbert_envelope(:,idx,:)));
    A=transpose(electrode_response);
    electrode_response=transpose(bsxfun(@minus, A, mean(A)));
    [a,b]=size(electrode_response);
    electrode_response_cell=mat2cell(electrode_response,a,word_length*ones(1,8));
    if strcmp(window_condition,'full')
        electrode_response_cell=electrode_response_cell;
    elseif strcmp(window_condition,'mean')
        electrode_response_cell=cellfun(@(x) mean(x,2),electrode_response_cell,'UniformOutput',false);
    end
    electrode_response_cell_1=[electrode_response_cell;session_examples_split{idx}];
    all_second_instances_signal=[all_second_instances_signal,electrode_response_cell_1];
 
    end
    % 
    first_instance=all_first_instances_signal;
    second_instance=all_second_instances_signal;
    
    A=first_instance(1,:);
    B=transpose(cell2mat(A));
    first_RDM=squareform(pdist(B,'correlation'));
    A=second_instance(1,:);
    B=transpose(cell2mat(A));
    second_RDM=squareform(pdist(B,'correlation'));
    %
    flatten_first_RDM=flatten_RDM(first_RDM,1);
    flatten_second_RDM=flatten_RDM(second_RDM,1);
    [R_corr,p_corr] = corrcoef(flatten_first_RDM,flatten_second_RDM);
    first_second_rep_corr=R_corr(1,2);
    % do permutation
    r = arrayfun(@(x) transpose(randperm(length(flatten_first_RDM))),[1:1000],'UniformOutput',false);
    permuted_flatten_second=cellfun(@(x) flatten_second_RDM(x),r,'UniformOutput',false);
    permuted_corr=cellfun(@(x) corrcoef(flatten_first_RDM,x),permuted_flatten_second,'UniformOutput',false);
    permuted_corr=cellfun(@(x) x(1,2),permuted_corr);
    all_RDMs(k).condtion=(condition);
    all_RDMs(k).first_sec_cor=first_second_rep_corr;
    all_RDMs(k).permutted=permuted_corr; 
    
end
% 
figure;set(gcf,'position',[1436,374,981,960]);
%
    colors = cbrewer('div', 'RdBu', 256);
    colors=flipud(colors);
    colormap(colors);
    %
    ax1=axes('position',[.1,.1,.7,.7]);
    permuted_corr_cell=[{all_RDMs(:).permutted}];
    y=mat2cell(1:length(permuted_corr_cell),1,ones(1,4));
    y=mat2cell([1:.5:2.5],1,ones(1,4));
    b_jitter_cell=cellfun(@(x,z) beeswarm((x*0+z)',x','hex','none',.3),permuted_corr_cell,y,'UniformOutput',false);
    
    bl=get(ax1,'Children');
    cla
    ax1.XLimMode='auto';
    ax1.XLimMode='auto';
    hl=plot([y{1}-.5,y{end}+.5],[0,0],'k--')
    hAnnotation=arrayfun(@(x) get(x,'Annotation'),hl);
    hLegendEntry = arrayfun(@(x) get(x,'LegendInformation'),hAnnotation);
    arrayfun(@(x) set(x,'IconDisplayStyle','off'),hLegendEntry);
    hold on
    bl=cellfun(@(x,y) scatter(x,y,5,'filled','markerfacecolor',[.5,.5,.5]),b_jitter_cell,permuted_corr_cell)
    hAnnotation=arrayfun(@(x) get(x,'Annotation'),bl);
    hLegendEntry = arrayfun(@(x) get(x,'LegendInformation'),hAnnotation);
    arrayfun(@(x) set(x,'IconDisplayStyle','off'),hLegendEntry);
    
    ax1.YLimMode='auto';
    hAnnotation=arrayfun(@(x) get(x,'Annotation'),bl);
    hLegendEntry = arrayfun(@(x) get(x,'LegendInformation'),hAnnotation);
    arrayfun(@(x) set(x,'IconDisplayStyle','off'),hLegendEntry);
    
    e=cellfun(@(x,y) errorbar(x,mean(y),std(y),'capsize',4,'LineWidth',2,'marker','o','color','k'),y,permuted_corr_cell);
    hAnnotation=arrayfun(@(x) get(x,'Annotation'),e);
    hLegendEntry = arrayfun(@(x) get(x,'LegendInformation'),hAnnotation);
    arrayfun(@(x) set(x,'IconDisplayStyle','off'),hLegendEntry);
    first_second_corr_cell=[{all_RDMs(:).first_sec_cor}];
    bl1=cellfun(@(x,y) scatter(x,y,50,'filled','markerfacecolor',[0,0,0],'markeredgecolor',[0,0,0],'linewidth',2),y,first_second_corr_cell);
    hAnnotation=arrayfun(@(x) get(x,'Annotation'),bl1);
    hLegendEntry = arrayfun(@(x) get(x,'LegendInformation'),hAnnotation);
    arrayfun(@(x) set(x,'IconDisplayStyle','off'),hLegendEntry);
    set(ax1,'xtick',cell2mat(y))
    ax1.XAxis.Color=[1,1,1];
    ax1.XAxis.TickLabels=conditions;
    ax1.XAxis.Visible = 'off';
    ax1.YLabel.String='Correlation';
    ax1.Title.String={'Correlation between 1st and 2nd RDM','vs 1st and permuted 2nd RDM',electrode_group};
   
    y_lim=ax1.YLim;
    x_lim=ax1.XLim;
    ax1.YLim=[-1.05*max(y_lim),1.05*max(y_lim)];
    cellfun(@(x,y) text(x,min(ax1.YLim),y,'horizontalalignment','center','fontsize',12,'fontweight','bold'),y,conditions);
    ax1.Box='off'
    print(gcf, '-bestfit','-dpdf', strcat(analysis_path,info.subject,'/',info.subject,'_SWJN_RDM_correlation_',electrode_group,'_elec.pdf'));
    



