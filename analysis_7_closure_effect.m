
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
i=1 
% load node info
[node_cell,node_table]=generate_node_table_eh(strcat(d_nmerge(i).folder,'/',d_nmerge(i).name));
% 
plot_width=.9/length(d);
plot_length=.87;
%%  find unique pattern, make a table of each sentence and for each word location what pattern it is member of 


nodes_closings_cell=cellfun(@(x) x(1:8), node_table.closed_nodes,'UniformOutput',false);
unique_nodes_closings_patterns=unique(cell2mat(nodes_closings_cell),'rows');
all_pattern_locations={};
for p=1:size(nodes_closings_cell,1)
    node_closing_pattern=nodes_closings_cell{p,:};
    
    pattern_locations={};
    for i=1:8
       
        [~,pattern_memebership_id]=ismember(unique_nodes_closings_patterns(:,1:i),node_closing_pattern(1:i),'rows');
        pattern_locations=[pattern_locations,find(pattern_memebership_id)];
    end
    all_pattern_locations=[all_pattern_locations;[{node_closing_pattern},pattern_locations]];
end
%% 
%% 
find_index= @(input_cell) find(cell2mat(cellfun(@(x) ~isempty(x), input_cell,'UniformOutput',false)));
all_example_sentence_locations=[];
all_sentence_pattern=[];
all_closing_pattern=[];
sentence_electrode_with_langauge_accross_sessions=[];
words_electrode_with_langauge_accross_sessions=[];

for i=1:length(d)
    example_sentence_locations=[];
    fprintf('adding %s from %s \n',d(i).name, strcat(d(i).folder,'/',d(i).name));
    subj=load(strcat(d(i).folder,'/',d(i).name));
    subj_id=fieldnames(subj);
    subj=subj.(subj_id{1});
    data=subj.data;
    info=subj.info;
    language_electrode=info.language_responsive_electrodes;
     
    language_electrode_num=find(language_electrode);
    % step 1: extract electrodes with siginificant language response
    sentence_trial_index=~cellfun(@isempty,cellfun(@ (x) strfind(x,'S'),info.word_type,'UniformOutput',false));
    % sentence
    sentences=[data{sentence_trial_index}];
    sentence_gamma_band_envelope=cellfun(@(x) x(1:8),{sentences.signal_hilbert_zs_downsample_parsed},'UniformOutput',false);
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

%% fill the empy locations with nan 
electrode_pattern_section_cell_modified=electrode_pattern_section_cell;
empty_cell=find(cell2mat(cellfun(@isempty,electrode_pattern_section_cell,'UniformOutput',false)));
template_mat=electrode_pattern_section_cell{1,1,1}*nan;
for nn=empty_cell'
    electrode_pattern_section_cell_modified{nn}=template_mat;
end

max_size=max(max(cell2mat(cellfun(@(x) size(x,3),electrode_pattern_section_cell_modified,'UniformOutput',false))));
    % add nan to fill the matrix 
    % find non_max cells 
    non_max_cell=find(cell2mat(cellfun(@(x) size(x,3)~=max_size, electrode_pattern_section_cell_modified,'UniformOutput',false))) ;
for nn=non_max_cell'
        cell_response=[];
        cell_response=electrode_pattern_section_cell_modified{nn};
        cell_dummy_response=nan*ones(size(cell_response,1),size(cell_response,2),max_size-size(cell_response,3));
        electrode_pattern_section_cell_modified{nn}=cat(3,cell_response,cell_dummy_response);
end

%% create a distribution plot of channel/time/opening-closing

one_closing=(all_closing_pattern')'==1;
no_node_closing=((all_closing_pattern')'==0);
% first column doesnt count 
no_node_closing(:,1)=0;
% 
look_out_window=60; %200 ms
anchor_points=1:135:size(sentence_electrode_with_langauge_accross_sessions,3);
all_electrode_closings={};
% find
for k=1:size(sentence_electrode_with_langauge_accross_sessions,1)
electrode_closing=[];
electrode_pre_post=[];
for kk=1:size(sentence_electrode_with_langauge_accross_sessions,2)    
    for kkk=find(one_closing(kk,:))
        electrode_closing=[electrode_closing,...
            squeeze(sentence_electrode_with_langauge_accross_sessions(k,kk,anchor_points(kkk)+[-look_out_window+1:+look_out_window]))];
        electrode_pre_post=[electrode_pre_post;...
                            mean(squeeze(sentence_electrode_with_langauge_accross_sessions(k,kk,anchor_points(kkk)+[-look_out_window+1:0]))),...
                            mean(squeeze(sentence_electrode_with_langauge_accross_sessions(k,kk,anchor_points(kkk)+[1:look_out_window])))];
    end 
end 
all_electrode_closings=[all_electrode_closings;[{electrode_closing},{electrode_pre_post}]];
end


    
% 
all_electrode_no_closings={};
% find
for k=1:size(sentence_electrode_with_langauge_accross_sessions,1)
electrode_closing=[];
electrode_pre_post=[];
for kk=1:size(sentence_electrode_with_langauge_accross_sessions,2)    
    for kkk=find(no_node_closing(kk,:))
        electrode_closing=[electrode_closing,...
            squeeze(sentence_electrode_with_langauge_accross_sessions(k,kk,anchor_points(kkk)+[-look_out_window+1:+look_out_window]))];
        electrode_pre_post=[electrode_pre_post;...
                            mean(squeeze(sentence_electrode_with_langauge_accross_sessions(k,kk,anchor_points(kkk)+[-look_out_window+1:0]))),...
                            mean(squeeze(sentence_electrode_with_langauge_accross_sessions(k,kk,anchor_points(kkk)+[1:look_out_window])))];
    end 
end 
all_electrode_no_closings=[all_electrode_no_closings;[{electrode_closing},{electrode_pre_post}]];
end
%% for each trial find cases of closing and non closing and plot them on a 2d 
all_electrode_closing_no_closings={};
for k=1:size(sentence_electrode_with_langauge_accross_sessions,1)
electrode_closing=[];
electrode_no_closing=[];
for kk=1:size(sentence_electrode_with_langauge_accross_sessions,2)
    % find index for closing an non closing 
    closings=find(all_closing_pattern(kk,:)==1);
    closings=closings(closings<8);
    closings=closings(randperm(length(closings)));
    non_closings=find(all_closing_pattern(kk,:)==0);
    non_closings=non_closings(non_closings>1);
    non_closings=non_closings(randperm(length(non_closings)));
    if ~isempty(closings) & ~isempty(non_closings)
    closing_section=squeeze(sentence_electrode_with_langauge_accross_sessions(k,kk,anchor_points(closings(1))+[-look_out_window+1:+look_out_window]));
    non_closing_section=squeeze(sentence_electrode_with_langauge_accross_sessions(k,kk,anchor_points(non_closings(1))+[-look_out_window+1:+look_out_window]));
    electrode_closing=[electrode_closing,closing_section];
    electrode_no_closing=[electrode_no_closing,non_closing_section];
    end 
end 
all_electrode_closing_no_closings=[all_electrode_closing_no_closings;[{electrode_closing},{electrode_no_closing}]];
end

%% for each electrode plot 

close all;
colors=cbrewer('qual','Set1',10);
num_rows=4;
num_columns=4;
total_plots=num_rows*num_columns;
for i=1:size(all_electrode_closing_no_closings,1)
    channel_response=all_electrode_closing_no_closings(i,:);
    channel_pre=cellfun(@(x) mean(x(1:60,:),1), channel_response,'UniformOutput',false);
    channel_post=cellfun(@(x) mean(x(61:end,:),1), channel_response,'UniformOutput',false);
    figure(fix((i-1)/total_plots)+1);
    set(gcf,'Position',[-882 449 877 1325]);
    sub_title=sprintf('subj: %s,\n electrodes: %d ', info.subject ,(i));
    a=subplot(num_rows,num_columns,i-total_plots*fix((i-1)/total_plots));
   
    h=scatter(channel_pre{1},channel_post{1},10);
    plot(channel_pre{1},channel_post{1},'LineStyle','none','Marker','o',...
    'MarkerSize',5,'MarkerFaceColor','r','MarkerEdgeColor',[1,1,1]);
    hold on;; 
    set(gca,'XAxisLocation','origin','YAxisLocation','origin','box','off');
%
    hold on 
    plot(channel_pre{2},channel_post{2},'LineStyle','none','Marker','o',...
    'MarkerSize',5,'MarkerFaceColor',[.5,.5,1],'MarkerEdgeColor',[1,1,1]);
    
   daspect([1,1,1]);
    plot(get(gca,'xlim'),get(gca,'xlim'),'k--')
    
    set(a,'color','none','layer','top')

end  
