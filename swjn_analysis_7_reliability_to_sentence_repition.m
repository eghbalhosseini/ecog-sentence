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
analysis_path='/Users/eghbalhosseiniasl1/MyData/ECoG-sentence/analysis/swjn_analysis_7_reliability_to_sentence_repetition/';
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

session_sentence_examples={};
session_wordlist_examples={};
sentence_sentence_locations=[];
session_sentence_hilbert_band_envelope_buildup_elec_tensor=[];
session_sentence_hilbert_zs_band_envelope_buildup_elec_tensor=[];
session_nonwords_hilbert_band_envelope_buildup_elec_tensor=[];
session_words_hilbert_band_envelope_buildup_elec_tensor=[];
session_jabberwocky_hilbert_band_envelope_buildup_elec_tensor=[];

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
    hilbert_band_envelope_buildup_elec_tensor=hilbert_band_envelope_tensor(find(language_electrode),:,:);

    session_sentence_hilbert_band_envelope_buildup_elec_tensor=cat(2,session_sentence_hilbert_band_envelope_buildup_elec_tensor,hilbert_band_envelope_buildup_elec_tensor);
    % find sentence locations 
    example_sentence=cellfun(@(x) x(2:end),{sentences(:).trial_string},'UniformOutput',false);
    session_sentence_examples=[session_sentence_examples;example_sentence'];
    example_sentence_locations=cellfun(@(x) (regexpi([node_cell{:,1}],x)),example_sentence,'UniformOutput',false);
    example_sentence_locations=cell2mat(cellfun(@(x) find_index(x), example_sentence_locations, 'UniformOutput',false));
    sentence_sentence_locations=[sentence_sentence_locations;example_sentence_locations'];
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
    hilbert_band_envelope_buildup_elec_tensor=hilbert_band_envelope_tensor(find(language_electrode),:,:);
    session_words_hilbert_band_envelope_buildup_elec_tensor=cat(2,session_words_hilbert_band_envelope_buildup_elec_tensor,hilbert_band_envelope_buildup_elec_tensor);
    % 
    example_words=cellfun(@(x) x(2:end),{words(:).trial_string},'UniformOutput',false);
    session_wordlist_examples=[session_wordlist_examples;example_words'];
    % 
    fprintf('added %s from %s \n',d(i).name, strcat(d(i).folder,'/',d(i).name));
    clear data subj
end
%% 
session_sentence_examples_1=cellfun(@(x) convertCharsToStrings(x),session_sentence_examples,'UniformOutput',false);
session_sentence_examples_split=transpose(cellfun(@(x) strsplit(x,' '),session_sentence_examples,'UniformOutput',false));

unique_sentences=table2cell(unique(cell2table(session_sentence_examples_1)));
unique_sentences=cellfun(@(x) x{:},unique_sentences,'UniformOutput',false);
unique_sentence_location=cellfun(@(x) (regexpi(session_sentence_examples,x)),unique_sentences,'UniformOutput',false);
unique_sentence_location=cellfun(@(x) find_index(x), unique_sentence_location, 'UniformOutput',false);

%% 
session_wordlist_examples_1=cellfun(@(x) convertCharsToStrings(x),session_wordlist_examples,'UniformOutput',false);
session_wordlist_examples_split=transpose(cellfun(@(x) strsplit(x,' '),session_wordlist_examples,'UniformOutput',false));

unique_wordlist=table2cell(unique(cell2table(session_wordlist_examples_1)));
unique_wordlist=cellfun(@(x) x{:},unique_wordlist,'UniformOutput',false);
unique_wordlist_location=cellfun(@(x) (regexpi(session_wordlist_examples,x)),unique_wordlist,'UniformOutput',false);
unique_wordlist_location=cellfun(@(x) find_index(x), unique_wordlist_location, 'UniformOutput',false);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  get common sentence and wordlists 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
first_repeat_sentences=[9:16];
second_repeat_sentences=[0:7]+17;
sentence_sentence_locations(second_repeat_sentences);
% 
session_sentence_examples_split=transpose(cellfun(@(x) strsplit(x,' '),session_sentence_examples,'UniformOutput',false));
temp=session_sentence_hilbert_band_envelope_buildup_elec_tensor;
session_sentence_hilbert_band_envelope_buildup_elec_tensor=temp;

%session_sentence_hilbert_band_envelope_buildup_elec_tensor=randn(size(temp,1),size(temp,2),size(temp,3));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot average for pre-post in word closing 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sampling_coeff=450/135;
win_width=400:50:450;%ms
win_width=450;%ms
win_step=0;%ms
win_width_sample=floor(win_width./sampling_coeff);%ms
win_step_sample=floor(win_step./sampling_coeff);%ms

num_rows=2;
num_columns=length(win_width);

figure;
set(gcf,'position',[1442 1 1275 1021])

for k=1:length(win_width)
    all_second_sentence_instances_ave=[];
    all_first_sentence_instances_ave=[];
    for i=1:length(first_repeat_sentences)
        first_repeat_index=first_repeat_sentences(i);
        second_repeat_index=second_repeat_sentences(i);
        if sentence_sentence_locations(first_repeat_index)==sentence_sentence_locations(second_repeat_index)
            electrode_first_sentence_response=double(squeeze(session_sentence_hilbert_band_envelope_buildup_elec_tensor(:,first_repeat_index,:)));
            electrode_first_sentence_response_cell=mat2cell(electrode_first_sentence_response,ones(1,size(electrode_first_sentence_response,1)),size(electrode_first_sentence_response,2)/8*ones(1,8));
            %electrode_first_sentence_response_cell_window=cellfun(@(x) double(x(win_width_sample(k)+[0:win_step_sample])),electrode_first_sentence_response_cell,'UniformOutput',false);
            electrode_first_sentence_response_cell_window=cellfun(@(x) double(x(:)),electrode_first_sentence_response_cell,'UniformOutput',false);
            electrode_first_sentence_response_cell_ave=cellfun(@(x) nanmean(x),electrode_first_sentence_response_cell_window,'UniformOutput',false);
            electrode_first_sentence_response_cell_ave_cell=mat2cell(cell2mat(electrode_first_sentence_response_cell_ave),size(electrode_first_sentence_response_cell_ave,1),ones(1,size(electrode_first_sentence_response_cell_ave,2)));
            electrode_first_sentence_response_cell_ave_cell=[electrode_first_sentence_response_cell_ave_cell;session_sentence_examples_split{first_repeat_index}];
            %
            electrode_second_sentence_response=double(squeeze(session_sentence_hilbert_band_envelope_buildup_elec_tensor(:,second_repeat_index,:)));
            electrode_second_sentence_response_cell=mat2cell(electrode_second_sentence_response,ones(1,size(electrode_second_sentence_response,1)),size(electrode_second_sentence_response,2)/8*ones(1,8));
            %electrode_second_sentence_response_cell_window=cellfun(@(x) double(x(win_width_sample(k)+[0:win_step_sample])),electrode_second_sentence_response_cell,'UniformOutput',false);
            electrode_second_sentence_response_cell_window=cellfun(@(x) double(x(:)),electrode_second_sentence_response_cell,'UniformOutput',false);
            electrode_second_sentence_response_cell_ave=cellfun(@(x) nanmean(x),electrode_second_sentence_response_cell_window,'UniformOutput',false);
            electrode_second_sentence_response_cell_ave_cell=mat2cell(cell2mat(electrode_second_sentence_response_cell_ave),size(electrode_second_sentence_response_cell_ave,1),ones(1,size(electrode_second_sentence_response_cell_ave,2)));
            electrode_second_sentence_response_cell_ave_cell=[electrode_second_sentence_response_cell_ave_cell;session_sentence_examples_split{second_repeat_index}];
            all_first_sentence_instances_ave=[all_first_sentence_instances_ave,electrode_first_sentence_response_cell_ave_cell];
            all_second_sentence_instances_ave=[all_second_sentence_instances_ave,electrode_second_sentence_response_cell_ave_cell];
            
        end      
    end
% 
first_second_repetition_correlation=[];
first_first_repetition_correlation=[];
second_second_repetition_correlation=[];


for i=1:size(all_first_sentence_instances_ave,2)
    cross_correlation={};
    within_cross_correlation={};
    cross_correlation=cellfun(@(x) corrcoef(x,all_first_sentence_instances_ave{1,i}),all_second_sentence_instances_ave(1,:),'UniformOutput',false);
    within_cross_correlation=cellfun(@(x) corrcoef(x,all_first_sentence_instances_ave{1,i}),all_first_sentence_instances_ave(1,:),'UniformOutput',false);
    second_within_cross_correlation=cellfun(@(x) corrcoef(x,all_second_sentence_instances_ave{1,i}),all_second_sentence_instances_ave(1,:),'UniformOutput',false);
    correlation_row=cellfun(@(x) x(1,2),cross_correlation);
    within_correlation_row=cellfun(@(x) x(1,2),within_cross_correlation);
    second_within_correlation_row=cellfun(@(x) x(1,2),second_within_cross_correlation);
    first_second_repetition_correlation=[first_second_repetition_correlation;correlation_row];
    first_first_repetition_correlation=[first_first_repetition_correlation;within_correlation_row];
    second_second_repetition_correlation=[second_second_repetition_correlation;second_within_correlation_row];
end


axis_min=min([first_first_repetition_correlation(:);first_second_repetition_correlation(:)]);
axis_max=max([first_first_repetition_correlation(:);first_second_repetition_correlation(:)]);

diag_mean=nanmean(diag(first_first_repetition_correlation));
tria_first_first_repetition_correlation=triu(first_first_repetition_correlation,1);
tria_ones=triu(ones(size(tria_first_first_repetition_correlation)),1);
tria_ones(tria_ones==0)=nan;
first_first_repetition_correlation_non_diag=first_first_repetition_correlation.*tria_ones;
off_diag_mean=nanmean(first_first_repetition_correlation_non_diag(:));

colors = cbrewer('div', 'RdBu', 128);
colors=flipud(colors);
colormap(colors);
subplot(2,length(win_width),(k-1)+1)
imagesc(first_first_repetition_correlation,[axis_min,axis_max]);
daspect([1,1,1]);
handles = colorbar;
handles.TickDirection = 'out';
handles.Box = 'off';

ylabel('first repetition');
xlabel('first repetition');
title({'Sentences : full window mean','S>N+buildup',['diag mean: ',num2str(diag_mean)],['off diag mean: ',num2str(off_diag_mean)]})

%diag_mean=nanmean(diag(first_second_repetition_correlation));
%tria_first_second_repetition_correlation=triu(first_second_repetition_correlation,1);
%tria_ones=triu(ones(size(tria_first_second_repetition_correlation)),1);
%tria_ones(tria_ones==0)=nan;

%first_second_repetition_correlation_non_diag=first_second_repetition_correlation.*tria_ones;
%off_diag_mean=nanmean(first_second_repetition_correlation_non_diag(:));
% 
diag_mean=nanmean(diag(first_second_repetition_correlation));
tria_ones=ones(size(first_second_repetition_correlation));
tria_ones(find(eye(size(tria_ones,1))))=nan;
first_second_repetition_correlation_non_diag=first_second_repetition_correlation.*tria_ones;
off_diag_mean=nanmean(first_second_repetition_correlation_non_diag(:));


subplot(2,length(win_width),length(win_width)+(k-1)+1)

imagesc(first_second_repetition_correlation,[axis_min,axis_max]);
daspect([1,1,1]);
handles = colorbar;
handles.TickDirection = 'out';
handles.Box = 'off';
ylabel('first repetition')
xlabel('second repetition')
title({'Sentences : full window mean','S>N+buildup',['diag mean: ',num2str(diag_mean)],['off diag mean: ',num2str(off_diag_mean)]})


     if ~exist(strcat(analysis_path,info.subject))
        mkdir(strcat(analysis_path,info.subject))
    end
    
    print(gcf, '-fillpage','-dpdf', strcat(analysis_path,info.subject,'/',info.subject,'_reliability_sentence_full_window_mean_buildup_elec.pdf')); 

    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot average for pre-post in word closing 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sampling_coeff=450/135;
win_width=400:50:450;%ms
win_width=450;%ms
win_step=0;%ms
win_width_sample=floor(win_width./sampling_coeff);%ms
win_step_sample=floor(win_step./sampling_coeff);%ms

num_rows=2;
num_columns=length(win_width);

figure;
set(gcf,'position',[1442 1 1275 1021])

for k=1:length(win_width)
    all_second_sentence_instances_signal=[];
    all_first_sentence_instances_signal=[];
    for i=1:length(first_repeat_sentences)
        first_repeat_index=first_repeat_sentences(i);
        second_repeat_index=second_repeat_sentences(i);
        if sentence_sentence_locations(first_repeat_index)==sentence_sentence_locations(second_repeat_index)
            electrode_first_sentence_response=double(squeeze(session_sentence_hilbert_band_envelope_buildup_elec_tensor(:,first_repeat_index,:)));
            electrode_first_sentence_response_cell=mat2cell(electrode_first_sentence_response,size(electrode_first_sentence_response,1),size(electrode_first_sentence_response,2)/8*ones(1,8));
            electrode_first_sentence_response_cell_1=[electrode_first_sentence_response_cell;session_sentence_examples_split{first_repeat_index}];
            all_first_sentence_instances_signal=[all_first_sentence_instances_signal,electrode_first_sentence_response_cell_1];
            %
            electrode_second_sentence_response=double(squeeze(session_sentence_hilbert_band_envelope_buildup_elec_tensor(:,second_repeat_index,:)));
            electrode_second_sentence_response_cell=mat2cell(electrode_second_sentence_response,size(electrode_second_sentence_response,1),size(electrode_second_sentence_response,2)/8*ones(1,8));
            electrode_second_sentence_response_cell_1=[electrode_second_sentence_response_cell;session_sentence_examples_split{second_repeat_index}];
            all_second_sentence_instances_signal=[all_second_sentence_instances_signal,electrode_second_sentence_response_cell_1];
            
        end      
    end
% 
first_second_repetition_correlation=[];
first_first_repetition_correlation=[];


for i=1:size(all_first_sentence_instances_signal,2)
    cross_correlation={};
    within_cross_correlation={};
    cross_correlation=cellfun(@(x) corrcoef(x,all_first_sentence_instances_signal{1,i}),all_second_sentence_instances_signal(1,:),'UniformOutput',false);
    within_cross_correlation=cellfun(@(x) corrcoef(x,all_first_sentence_instances_signal{1,i}),all_first_sentence_instances_signal(1,:),'UniformOutput',false);
    correlation_row=cellfun(@(x) x(1,2),cross_correlation);
    within_correlation_row=cellfun(@(x) x(1,2),within_cross_correlation);
    first_second_repetition_correlation=[first_second_repetition_correlation;correlation_row];
    first_first_repetition_correlation=[first_first_repetition_correlation;within_correlation_row];
end

axis_min=min([first_first_repetition_correlation(:);first_second_repetition_correlation(:)]);
axis_max=max([first_first_repetition_correlation(:);first_second_repetition_correlation(:)]);

diag_mean=nanmean(diag(first_first_repetition_correlation));
tria_first_first_repetition_correlation=triu(first_first_repetition_correlation,1);
tria_ones=triu(ones(size(tria_first_first_repetition_correlation)),1);
tria_ones(tria_ones==0)=nan;

first_first_repetition_correlation_non_diag=first_first_repetition_correlation.*tria_ones;
off_diag_mean=nanmean(first_first_repetition_correlation_non_diag(:));

colors = cbrewer('div', 'RdBu', 128);
colors=flipud(colors);
colormap(colors);
subplot(2,length(win_width),(k-1)+1)
imagesc(first_first_repetition_correlation,[axis_min,axis_max]);
daspect([1,1,1]);
handles = colorbar;
handles.TickDirection = 'out';
handles.Box = 'off';
ylabel('first repetition');
xlabel('first repetition');
title({'Sentences : full window','S>N+buildup',['diag mean: ',num2str(diag_mean)],['off diag mean: ',num2str(off_diag_mean)]})

% diag_mean=nanmean(diag(first_second_repetition_correlation));
% tria_first_second_repetition_correlation=triu(first_second_repetition_correlation,1);
% tria_ones=triu(ones(size(tria_first_second_repetition_correlation)),1);
% tria_ones(tria_ones==0)=nan;
% first_second_repetition_correlation_non_diag=first_second_repetition_correlation.*tria_ones;
% off_diag_mean=nanmean(first_second_repetition_correlation_non_diag(:));

diag_mean=nanmean(diag(first_second_repetition_correlation));
tria_ones=ones(size(first_second_repetition_correlation));
tria_ones(find(eye(size(tria_ones,1))))=nan;

first_second_repetition_correlation_non_diag=first_second_repetition_correlation.*tria_ones;
off_diag_mean=nanmean(first_second_repetition_correlation_non_diag(:));

subplot(2,length(win_width),length(win_width)+(k-1)+1)

imagesc(first_second_repetition_correlation,[axis_min,axis_max]);
daspect([1,1,1]);
handles = colorbar;
handles.TickDirection = 'out';
handles.Box = 'off';
ylabel('first repetition')
xlabel('second repetition')
title({'Sentences : full window','S>N+buildup',['diag mean: ',num2str(diag_mean)],['off diag mean: ',num2str(off_diag_mean)]})


 if ~exist(strcat(analysis_path,info.subject))
        mkdir(strcat(analysis_path,info.subject))
    end
    
    print(gcf, '-fillpage','-dpdf', strcat(analysis_path,info.subject,'/',info.subject,'_reliability_sentence_full_window_buildup_elec.pdf')); 
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% wordlist repetition 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot average for pre-post in word closing 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
unique_wordlist_location
sampling_coeff=450/135;
win_width=400:50:450;%ms
win_width=450;%ms
win_step=0;%ms
win_width_sample=floor(win_width./sampling_coeff);%ms
win_step_sample=floor(win_step./sampling_coeff);%ms

num_rows=2;
num_columns=length(win_width);

figure;
set(gcf,'position',[1442 1 1275 1021])
word_repetition_location=cell2mat(unique_wordlist_location(cellfun(@length,unique_wordlist_location)>1)');
first_repeat_wordss=word_repetition_location(1,:);
second_repeat_wordss=word_repetition_location(2,:);
for k=1:length(win_width)
    all_second_words_instances_ave=[];
    all_first_words_instances_ave=[];
    for i=1:length(first_repeat_wordss)
        first_repeat_index=first_repeat_wordss(i);
        second_repeat_index=second_repeat_wordss(i);
        
        electrode_first_words_response=double(squeeze(session_words_hilbert_band_envelope_buildup_elec_tensor(:,first_repeat_index,:)));
        electrode_first_words_response_cell=mat2cell(electrode_first_words_response,ones(1,size(electrode_first_words_response,1)),size(electrode_first_words_response,2)/8*ones(1,8));
        %electrode_first_words_response_cell_window=cellfun(@(x) double(x(win_width_sample(k)+[0:win_step_sample])),electrode_first_words_response_cell,'UniformOutput',false);
        electrode_first_words_response_cell_window=cellfun(@(x) double(x(:)),electrode_first_words_response_cell,'UniformOutput',false);
        electrode_first_words_response_cell_ave=cellfun(@(x) nanmean(x),electrode_first_words_response_cell_window,'UniformOutput',false);
        electrode_first_words_response_cell_ave_cell=mat2cell(cell2mat(electrode_first_words_response_cell_ave),size(electrode_first_words_response_cell_ave,1),ones(1,size(electrode_first_words_response_cell_ave,2)));
        electrode_first_words_response_cell_ave_cell=[electrode_first_words_response_cell_ave_cell;session_wordlist_examples_split{first_repeat_index}];
        %
        electrode_second_words_response=double(squeeze(session_words_hilbert_band_envelope_buildup_elec_tensor(:,second_repeat_index,:)));
        electrode_second_words_response_cell=mat2cell(electrode_second_words_response,ones(1,size(electrode_second_words_response,1)),size(electrode_second_words_response,2)/8*ones(1,8));
        %electrode_second_words_response_cell_window=cellfun(@(x) double(x(win_width_sample(k)+[0:win_step_sample])),electrode_second_words_response_cell,'UniformOutput',false);
        electrode_second_words_response_cell_window=cellfun(@(x) double(x(:)),electrode_second_words_response_cell,'UniformOutput',false);
        electrode_second_words_response_cell_ave=cellfun(@(x) nanmean(x),electrode_second_words_response_cell_window,'UniformOutput',false);
        electrode_second_words_response_cell_ave_cell=mat2cell(cell2mat(electrode_second_words_response_cell_ave),size(electrode_second_words_response_cell_ave,1),ones(1,size(electrode_second_words_response_cell_ave,2)));
        electrode_second_words_response_cell_ave_cell=[electrode_second_words_response_cell_ave_cell;session_wordlist_examples_split{second_repeat_index}];
        all_first_words_instances_ave=[all_first_words_instances_ave,electrode_first_words_response_cell_ave_cell];
        all_second_words_instances_ave=[all_second_words_instances_ave,electrode_second_words_response_cell_ave_cell];
        
    
    end
% 
first_second_repetition_correlation=[];
first_first_repetition_correlation=[];


for i=1:size(all_first_words_instances_ave,2)
    cross_correlation={};
    within_cross_correlation={};
    cross_correlation=cellfun(@(x) corrcoef(x,all_first_words_instances_ave{1,i}),all_second_words_instances_ave(1,:),'UniformOutput',false);
    within_cross_correlation=cellfun(@(x) corrcoef(x,all_first_words_instances_ave{1,i}),all_first_words_instances_ave(1,:),'UniformOutput',false);
    correlation_row=cellfun(@(x) x(1,2),cross_correlation);
    within_correlation_row=cellfun(@(x) x(1,2),within_cross_correlation);
    first_second_repetition_correlation=[first_second_repetition_correlation;correlation_row];
    first_first_repetition_correlation=[first_first_repetition_correlation;within_correlation_row];
end

axis_min=min([first_first_repetition_correlation(:);first_second_repetition_correlation(:)]);
axis_max=max([first_first_repetition_correlation(:);first_second_repetition_correlation(:)]);

diag_mean=nanmean(diag(first_first_repetition_correlation));
tria_first_first_repetition_correlation=triu(first_first_repetition_correlation,1);
tria_ones=triu(ones(size(tria_first_first_repetition_correlation)),1);
tria_ones(tria_ones==0)=nan;

first_first_repetition_correlation_non_diag=first_first_repetition_correlation.*tria_ones;
off_diag_mean=nanmean(first_first_repetition_correlation_non_diag(:));

colors = cbrewer('div', 'RdBu', 128);
colors=flipud(colors);
colormap(colors);
subplot(2,length(win_width),(k-1)+1)
imagesc(first_first_repetition_correlation,[axis_min,axis_max]);
daspect([1,1,1]);
handles = colorbar;
handles.TickDirection = 'out';
handles.Box = 'off';

ylabel('first repetition');
xlabel('first repetition');
title({'Word list : full window mean','S>N+buildup',['diag mean: ',num2str(diag_mean)],['off diag mean: ',num2str(off_diag_mean)]})

diag_mean=nanmean(diag(first_second_repetition_correlation));
tria_ones=ones(size(first_second_repetition_correlation));
tria_ones(find(eye(size(tria_ones,1))))=nan;

first_second_repetition_correlation_non_diag=first_second_repetition_correlation.*tria_ones;
off_diag_mean=nanmean(first_second_repetition_correlation_non_diag(:));


subplot(2,length(win_width),length(win_width)+(k-1)+1)

imagesc(first_second_repetition_correlation,[axis_min,axis_max]);
daspect([1,1,1]);
handles = colorbar;
handles.TickDirection = 'out';
handles.Box = 'off';
ylabel('first repetition')
xlabel('second repetition')
title({'Word list : full window mean','S>N+buildup',['diag mean: ',num2str(diag_mean)],['off diag mean: ',num2str(off_diag_mean)]})

 if ~exist(strcat(analysis_path,info.subject))
        mkdir(strcat(analysis_path,info.subject))
    end
    
    print(gcf, '-fillpage','-dpdf', strcat(analysis_path,info.subject,'/',info.subject,'_reliability_wordlist_full_window_mean_buildup_elec.pdf')); 

    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot across all 
sampling_coeff=450/135;
win_width=400:50:450;%ms
win_width=450;%ms
win_step=0;%ms
win_width_sample=floor(win_width./sampling_coeff);%ms
win_step_sample=floor(win_step./sampling_coeff);%ms

num_rows=2;
num_columns=length(win_width);

figure;
set(gcf,'position',[1442 1 1275 1021])

for k=1:length(win_width)
    all_second_words_instances_signal=[];
    all_first_words_instances_signal=[];
    for i=1:length(first_repeat_wordss)
        first_repeat_index=first_repeat_wordss(i);
        second_repeat_index=second_repeat_wordss(i);
        
        electrode_first_words_response=double(squeeze(session_words_hilbert_band_envelope_buildup_elec_tensor(:,first_repeat_index,:)));
        electrode_first_words_response_cell=mat2cell(electrode_first_words_response,size(electrode_first_words_response,1),size(electrode_first_words_response,2)/8*ones(1,8));
        electrode_first_words_response_cell_1=[electrode_first_words_response_cell;session_wordlist_examples_split{first_repeat_index}];
        all_first_words_instances_signal=[all_first_words_instances_signal,electrode_first_words_response_cell_1];
        %
        electrode_second_words_response=double(squeeze(session_words_hilbert_band_envelope_buildup_elec_tensor(:,second_repeat_index,:)));
        electrode_second_words_response_cell=mat2cell(electrode_second_words_response,size(electrode_second_words_response,1),size(electrode_second_words_response,2)/8*ones(1,8));
        electrode_second_words_response_cell_1=[electrode_second_words_response_cell;session_wordlist_examples_split{second_repeat_index}];
        all_second_words_instances_signal=[all_second_words_instances_signal,electrode_second_words_response_cell_1];
        
    end
% 
first_second_repetition_correlation=[];
first_first_repetition_correlation=[];


for i=1:size(all_first_words_instances_signal,2)
    cross_correlation={};
    within_cross_correlation={};
    cross_correlation=cellfun(@(x) corrcoef(x,all_first_words_instances_signal{1,i}),all_second_words_instances_signal(1,:),'UniformOutput',false);
    within_cross_correlation=cellfun(@(x) corrcoef(x,all_first_words_instances_signal{1,i}),all_first_words_instances_signal(1,:),'UniformOutput',false);
    correlation_row=cellfun(@(x) x(1,2),cross_correlation);
    within_correlation_row=cellfun(@(x) x(1,2),within_cross_correlation);
    first_second_repetition_correlation=[first_second_repetition_correlation;correlation_row];
    first_first_repetition_correlation=[first_first_repetition_correlation;within_correlation_row];
end

axis_min=min([first_first_repetition_correlation(:);first_second_repetition_correlation(:)]);
axis_max=max([first_first_repetition_correlation(:);first_second_repetition_correlation(:)]);

diag_mean=nanmean(diag(first_first_repetition_correlation));
tria_first_first_repetition_correlation=triu(first_first_repetition_correlation,1);
tria_ones=triu(ones(size(tria_first_first_repetition_correlation)),1);
tria_ones(tria_ones==0)=nan;

first_first_repetition_correlation_non_diag=first_first_repetition_correlation.*tria_ones;
off_diag_mean=nanmean(first_first_repetition_correlation_non_diag(:));

colors = cbrewer('div', 'RdBu', 128);
colors=flipud(colors);
colormap(colors);
subplot(2,length(win_width),(k-1)+1)
imagesc(first_first_repetition_correlation,[axis_min,axis_max]);
daspect([1,1,1]);
handles = colorbar;
handles.TickDirection = 'out';
handles.Box = 'off';

ylabel('first repetition');
xlabel('first repetition');
title({'Word list: full window','S>N+buildup',['diag mean: ',num2str(diag_mean)],['off diag mean: ',num2str(off_diag_mean)]})

diag_mean=nanmean(diag(first_second_repetition_correlation));
tria_ones=ones(size(first_second_repetition_correlation));
tria_ones(find(eye(size(tria_ones,1))))=nan;

first_second_repetition_correlation_non_diag=first_second_repetition_correlation.*tria_ones;
off_diag_mean=nanmean(first_second_repetition_correlation_non_diag(:));

subplot(2,length(win_width),length(win_width)+(k-1)+1)

imagesc(first_second_repetition_correlation,[axis_min,axis_max]);
daspect([1,1,1]);
handles = colorbar;
handles.TickDirection = 'out';
handles.Box = 'off';

ylabel('first repetition')
xlabel('second repetition')
title({'Word list : full window','S>N+buildup',['diag mean: ',num2str(diag_mean)],['off diag mean: ',num2str(off_diag_mean)]})


     if ~exist(strcat(analysis_path,info.subject))
        mkdir(strcat(analysis_path,info.subject))
    end
    
    print(gcf, '-fillpage','-dpdf', strcat(analysis_path,info.subject,'/',info.subject,'_reliability_wordlist_full_window_buildup_elec.pdf')); 

    
end







