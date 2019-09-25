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
d= dir([data_path,'/**/AMC*_crunched.mat']);
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
sentence_ave_electrode_with_langauge_accross_sessions=[];
sentence_ave_hilbert_electrode_with_langauge_accross_sessions=[];
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
    sentence_gamma_band_envelope=[sentence_gamma_band_envelope{:,:}];
    sentence_gamma_band_envelope=cellfun(@transpose,sentence_gamma_band_envelope,'UniformOutput',false);
    sentence_gamma_band_envelope_tensor=cell2mat(permute(sentence_gamma_band_envelope,[1,3,2]));
    sentence_gamma_band_envelope_tensor=permute(sentence_gamma_band_envelope_tensor,[2,3,1]);
    electrodes_with_language_response=sentence_gamma_band_envelope_tensor(find(language_electrode),:,:);
    sentence_electrode_with_langauge_accross_sessions=cat(2,sentence_electrode_with_langauge_accross_sessions,electrodes_with_language_response);
    % 
    sentence_gamma_band_envelope_ave=cellfun(@(x) x(1:8),{sentences.signal_ave_hilbert_zs_downsample_parsed},'UniformOutput',false);
    sentence_gamma_band_ave_envelope=[sentence_gamma_band_envelope_ave{:,:}];
    sentence_gamma_band_ave_envelope_tensor=cell2mat(permute(sentence_gamma_band_ave_envelope,[3,2,1]));
    ave_electrodes_with_language_response=sentence_gamma_band_ave_envelope_tensor(find(language_electrode),:,:);
    sentence_ave_electrode_with_langauge_accross_sessions=cat(2,sentence_ave_electrode_with_langauge_accross_sessions,ave_electrodes_with_language_response);
    % 
    sentence_hilbert_band_envelope_ave=cellfun(@(x) x(1:8),{sentences.signal_ave_hilbert_downsample_parsed},'UniformOutput',false);
    sentence_hilbert_band_ave_envelope=[sentence_hilbert_band_envelope_ave{:,:}];
    sentence_hilbert_band_ave_envelope_tensor=cell2mat(permute(sentence_hilbert_band_ave_envelope,[3,2,1]));
    ave_hilbert_electrodes_with_language_response=sentence_hilbert_band_ave_envelope_tensor(find(language_electrode),:,:);
    sentence_ave_hilbert_electrode_with_langauge_accross_sessions=cat(2,sentence_ave_hilbert_electrode_with_langauge_accross_sessions,ave_hilbert_electrodes_with_language_response);
    
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
    
    fprintf('added %s from %s \n',d(i).name, strcat(d(i).folder,'/',d(i).name));
    
end

%% get the nnmf of the pattern 
%[coeff,score,latent,tsquared,explained,mu] = pca(all_sentence_pattern);

[W_pattern,H_pattern,D_pattern] = nnmf(all_sentence_pattern,3);



%% 
colors=cbrewer('qual','Set1',10);

for i=1:length(language_electrode_num)
    electrode_response=squeeze(sentence_ave_hilbert_electrode_with_langauge_accross_sessions(i,:,:));
    sub_title=sprintf('subj: %s,\n electrodes: %d ', info.subject ,language_electrode_num(i));
    word_position=repmat(1:size(electrode_response,2),[size(electrode_response,1),1]);
    perturbed_word_position=word_position+.2*rand(size(electrode_response,1),size(electrode_response,2))-.1;
    figure('position',[-1685 519 492 1131])
    a=subplot(4,3,[1:3]);
    h=scatter(perturbed_word_position(:),electrode_response(:),4);
    set(h,'CData',colors(6,:),'MarkerFacecolor',colors(1,:),'MarkerFaceAlpha',.3)
    hold on
    e=plot(mean(word_position,1),mean(electrode_response,1),...
        '-s','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',2,'color','k');
    xlim([min(min(word_position))-.5,max(max(word_position))+.5])
    set(a,'XTick',[1:8]);
    e.Color=[0,0,0];
    y_quantile=quantile(electrode_response(:),50);
    ylim([y_quantile(1),y_quantile(end)]);
    title(sub_title)
    [coeff,score,latent,tsquared,explained,mu] = pca(double(electrode_response-repmat(mean(electrode_response,2),1,size(electrode_response,2))))
    for kk=1:3
    a=subplot(4,3,[3+kk]);
    electrode_reconstructed=score(:,kk)*coeff(:,kk)'+repmat(mu,size(score,1),1);
    h=scatter(perturbed_word_position(:),electrode_reconstructed(:),4);
    set(h,'CData',colors(6,:),'MarkerFacecolor',colors(1,:),'MarkerFaceAlpha',.3)
    hold on
    e=plot(mean(word_position,1),mean(electrode_reconstructed,1),...
        '-s','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',2,'color','k');
    xlim([min(min(word_position))-.5,max(max(word_position))+.5])
    set(a,'XTick',[1:8])
    e.Color=[0,0,0];
    end 
    
    
    
    a=subplot(4,3,[7:9]);
    electrode_response=squeeze(sentence_ave_hilbert_electrode_with_langauge_accross_sessions(i,:,:));
    [W_data,H_data,D_data] = nnmf(electrode_response,3); 
    electrode_reconstructed=W_data*H_data;
    h=scatter(perturbed_word_position(:),electrode_reconstructed(:),4);
    set(h,'CData',colors(6,:),'MarkerFacecolor',colors(1,:),'MarkerFaceAlpha',.3)
    hold on
    e=plot(mean(word_position,1),mean(electrode_reconstructed,1),...
        '-s','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',2,'color','k');
    xlim([min(min(word_position))-.5,max(max(word_position))+.5])
    set(a,'XTick',[1:8])
    e.Color=[0,0,0];
    
    %
    for kk=1:size(H_data,1)
    a=subplot(4,3,[9+kk]);
    electrode_reconstructed=W_data(:,kk)*H_data(kk,:);
    h=scatter(perturbed_word_position(:),electrode_reconstructed(:),4);
    set(h,'CData',colors(6,:),'MarkerFacecolor',colors(1,:),'MarkerFaceAlpha',.3)
    hold on
    e=plot(mean(word_position,1),mean(electrode_reconstructed,1),...
        '-s','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',2,'color','k');
    xlim([min(min(word_position))-.5,max(max(word_position))+.5])
    set(a,'XTick',[1:8])
    e.Color=[0,0,0];
    end 
    
    
    
    
end 


