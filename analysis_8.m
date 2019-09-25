clear all 
close all 
home
%% 
if 1
    fprintf('adding basic ecog tools to path \n');
    addpath('~/MyCodes/basic-ecog-tools/');
    addpath('~/MyCodes/ecog-sentence/');
    addpath(genpath('~/MyCodes/basic-ecog-tools/activeBrain'));
    addpath(genpath('~/MyCodes/basic-ecog-tools/ecog-filters'));
    addpath(genpath('~/MyCodes/basic-ecog-tools/mex'));
end 

%% 
data_path='/Users/eghbalhosseiniasl1/MyData/ecog-sentence';
d= dir([data_path,'/**/AMC026*_crunched.mat']);
fprintf(' %d .mat files were found \n', length(d));

% find the node information
d_nmerge= dir([data_path,'/**/*nmerge.txt']);
fprintf(' %d nmerge files were found \n', length(d_nmerge));
i=1;
% load node info
[node_cell,node_table]=generate_node_table_eh(strcat(d_nmerge(i).folder,'/',d_nmerge(i).name));
% 
plot_width=.9/length(d);
plot_length=.87;
%% there are 31 unique patterns, make a table of each sentence and for each word location what pattern it is member of 
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

%% for the location in the trials for all unique patterns 
open_nodes_cell=cellfun(@(x) x(1:8), node_table.open_nodes,'UniformOutput',false);
unique_pattern_locations={};
for p=1:size(unique_open_node_patterns,1)
    unique_open_node_pattern=unique_open_node_patterns(p,:)
    pattern_locations={};
    for i=1:8
        a=unique_open_node_pattern(1:i);
        pattern_locations=[pattern_locations,
        find(~cell2mat(cellfun(@isempty,cellfun(@(x) strfind(num2str(x),num2str(a)),open_nodes_cell,'UniformOutput',false),'UniformOutput',false)))];
    end
    unique_pattern_locations=[unique_pattern_locations;[unique_open_node_pattern,transpose(pattern_locations)]];
end
%% 
find_index= @(input_cell) find(cell2mat(cellfun(@(x) ~isempty(x), input_cell,'UniformOutput',false)));
all_example_sentence_locations=[];
all_sentence_pattern=[];
all_closing_pattern=[];
sentence_electrode_with_langauge_accross_sessions=[];
words_electrode_with_langauge_accross_sessions=[];
all_sentence_pattern_id=[];
for i=1:length(d)
    example_sentence_locations=[];
    fprintf('adding %s from %s \n',d(i).name, strcat(d(i).folder,'/',d(i).name));
    subj=load(strcat(d(i).folder,'/',d(i).name));
    subj_id=fieldnames(subj);
    subj=subj.(subj_id{1});
    data=subj.data;
    info=subj.info;
    language_electrode=info.ramp_electrodes;
     
    language_electrode_num=find(language_electrode);
    % step 1: extract electrodes with siginificant language response
    sentence_trial_index=~cellfun(@isempty,cellfun(@ (x) strfind(x,'S'),info.word_type,'UniformOutput',false));
    % sentence
    sentences=[data{sentence_trial_index}];% 
    sentence_gamma_band_envelope=cellfun(@(x) x(1:10),{sentences.signal_hilbert_zs_downsample_parsed},'UniformOutput',false);
    % creat a cell with wordposition(row)*time in trial(column) structure
    sentence_gamma_band_envelope=[sentence_gamma_band_envelope{:,:}];
    sentence_gamma_band_envelope=cellfun(@transpose,sentence_gamma_band_envelope,'UniformOutput',false);
    % make a words positions*channel* trial tensor 
    sentence_gamma_band_ave_envelope_tensor=cell2mat(permute(sentence_gamma_band_envelope,[1,3,2]));
    %append to langauge_channel*trial*words positions
    sentence_gamma_band_ave_envelope_tensor=permute(sentence_gamma_band_ave_envelope_tensor,[2,3,1]);
    electrodes_with_language_response=sentence_gamma_band_ave_envelope_tensor(find(language_electrode),:,:);
    sentence_electrode_with_langauge_accross_sessions=cat(2,sentence_electrode_with_langauge_accross_sessions,electrodes_with_language_response);
    % 
    
    
    % find sentence locations 
    example_sentence={sentences(:).trial_string};
    example_sentence=cellfun(@(x) x(2:end),example_sentence,'UniformOutput',false);
    % find trial location in node cell 
    example_sentence_locations=cellfun(@(x) (regexpi([node_cell{:,1}],x)),example_sentence,'UniformOutput',false);
    example_sentence_locations=cell2mat(cellfun(@(x) find_index(x), example_sentence_locations, 'UniformOutput',false));
 
    all_sentence_pattern=[all_sentence_pattern;
        cell2mat(cellfun(@(x) x(1:8),node_table.open_nodes(example_sentence_locations),'UniformOutput',false))];
    all_sentence_pattern_id=[all_sentence_pattern_id;cell2mat(cellfun(@(x) find(ismember(unique_open_node_patterns,x,'rows')),...
        cellfun(@(x) x(1:8),node_table.open_nodes(example_sentence_locations),'UniformOutput',false),'UniformOutput',false))];
    all_closing_pattern=[all_closing_pattern;
        cell2mat(cellfun(@(x) x(1:8),node_table.closed_nodes(example_sentence_locations),'UniformOutput',false))];
    %
    words_trial_index=~cellfun(@isempty,cellfun(@ (x) strfind(x,'W'),info.word_type,'UniformOutput',false));
    % words
    words=[data{words_trial_index}];
    words_gamma_band_envelope=cellfun(@(x) x(1:10),{words.signal_hilbert_zs_downsample_parsed},'UniformOutput',false);
    % creat a cell with wordposition(row)*time in trial(column) structure
    words_gamma_band_envelope=[words_gamma_band_envelope{:,:}];
    words_gamma_band_envelope=cellfun(@transpose,words_gamma_band_envelope,'UniformOutput',false);
    % make a words positions*channel* trial tensor
    words_gamma_band_ave_envelope_tensor=cell2mat(permute(words_gamma_band_envelope,[1,3,2]));
    %append to langauge_channel*trial*words positions
    words_gamma_band_ave_envelope_tensor=permute(words_gamma_band_ave_envelope_tensor,[2,3,1]);
    electrodes_with_language_response=words_gamma_band_ave_envelope_tensor(find(language_electrode),:,:);
    words_electrode_with_langauge_accross_sessions=cat(2,words_electrode_with_langauge_accross_sessions,electrodes_with_language_response);
    
    
    fprintf('added %s from %s \n',d(i).name, strcat(d(i).folder,'/',d(i).name));
    
end

%% 
% do it on the 
unique_subject_node_patterns=unique(all_sentence_pattern,'row');
unique_subject_node_patterns_cell=mat2cell(all_sentence_pattern,...
    ones(1,size(all_sentence_pattern,1)),8);
pattern_transitions=cell(size(unique_open_node_pattern_cell,2)-1);
for i=2:8
    step_unique_pattern=unique(unique_subject_node_patterns(:,1:i-1),'row');
    for k=1:size(step_unique_pattern,1)
        transitions=find(ismember(unique_subject_node_patterns(:,1:i-1),step_unique_pattern(k,:),'rows'));
        pattern_transitions{k,i-1}=[unique_subject_node_patterns(transitions,1:i),transitions];
    end
end 

subject_sentence_pattern_id=cell2mat(cellfun(@(x) find(ismember(unique_subject_node_patterns,x,'rows')),...
    unique_subject_node_patterns_cell,'UniformOutput',false))

%% 
num_rows=7;
num_columns=max(sum(~cellfun(@isempty,pattern_transitions(:,1:end-1))));
total_plots=num_rows*num_columns;
plot_height=(.8-.05)./num_rows;
plot_width=(0.8-.05)./num_columns;
look_out_window=100; %200 ms
anchor_points=1:135:size(sentence_electrode_with_langauge_accross_sessions,3);

colors = cbrewer('seq','YlGnBu',  5);
colors=flipud(colors);

for i=1:length(language_electrode_num)
    electrode_response=squeeze(sentence_electrode_with_langauge_accross_sessions(i,:,:));
    electrode_response_ave=nanmean(electrode_response,1);
    close all 
    figure('position',[-2261 452 2218 1291]);

    for k=1:size(pattern_transitions,2)-1
        pattern_transition=pattern_transitions(:,k);
        plot_num=find(~cellfun(@isempty,pattern_transition));
        %all_patterns=cell2mat(cellfun(@(x) x(:,end),pattern_transitions(plot_num,k),'UniformOutput',false));
        %pattern_to_take_mean=electrode_response(ismember(subject_sentence_pattern_id,all_patterns),...
       %             anchor_points(k-1)+[0:look_out_window]);
       ax={};
       y_lims=[];
       for l=plot_num'
            
            ax{l}=axes('position',[.02+(l-1)*plot_width+0.01*(l-1),.03*(k-1)+(k-1)*plot_height+.03,plot_width,plot_height-.01]);
            pattern=pattern_transition{l};
            tansitions_types=unique(pattern(:,end-1));
            drift=[0:3:20];
            for kk=1:size(tansitions_types,1)
                
                pattern_idx_to_take_mean=pattern(pattern(:,end-1)==tansitions_types(kk),end);
                pattern_to_take_mean=electrode_response(ismember(subject_sentence_pattern_id,pattern_idx_to_take_mean),...
                    anchor_points(size(pattern,2)-1)+[0:look_out_window]);
                bl = plot(anchor_points(size(pattern,2)-1)+[0:look_out_window],...
                    nanmean(pattern_to_take_mean,1), 'color', colors(kk,:),'linewidth',3);
                hold on
                a=nanmean(pattern_to_take_mean,1);
                text(anchor_points(size(pattern,2)-1)-drift(kk),...
                    double(a(1)),num2str(tansitions_types(kk)),...
                    'FontSize',12,'color',colors(kk,:),'VerticalAlignment','bottom','HorizontalAlignment','right','FontWeight','bold');   
            end
            y_lims=[y_lims;get(gca,'ylim')];
            bl = plot(anchor_points(size(pattern,2)-1)+[0:look_out_window],...
                    electrode_response_ave(anchor_points(size(pattern,2)-1)+[0:look_out_window]), 'k--');
            set(gca,'box','off');
            if l==1
             ax{l}.YAxis.Visible = 'on';
             
             %ax{l}.XTick=[anchor_points(size(pattern,2))+[0,look_out_window]];
             %ax{l}.XTickLabel= [anchor_points(size(pattern,2))+[0,look_out_window]];  
            else 
            ax{l}.YAxis.Visible = 'off';
            ax{l}.XTick=[];
            end 
            title(num2str((unique(pattern(:,1:end-2),'row'))));
       end
       new_bounds=quantile(y_lims(:),15);
       for l=plot_num'
           ax{l}.YLim=[new_bounds(1),new_bounds(end)];
       end 
    end
    dim = [.2 .1 .5, .1 ];
    str = strrep(strcat(info.subject,'_channel_',num2str(language_electrode_num(i))),'_', ' ');
    h=annotation('textbox',dim,'String',str,'Fontsize',14,'FontWeight','bold','linestyle','none');
    if ~exist(strcat('/Users/eghbalhosseiniasl1/MyData/ecog-sentence/transition-response/',info.subject))
        mkdir(strcat('/Users/eghbalhosseiniasl1/MyData/ecog-sentence/transition_response/',info.subject))
    end 
    print(gcf, '-depsc', strcat('/Users/eghbalhosseiniasl1/MyData/ecog-sentence/transition_response/',info.subject,'/',info.subject,'_channel_',num2str(language_electrode_num(i)),'_response_trnasition','.eps')); 
 
end 
