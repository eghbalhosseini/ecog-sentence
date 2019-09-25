clear all ;
close all ;
home;
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
analysis_path='/Users/eghbalhosseiniasl1/MyData/ecog-sentence/analysis/analysis_3_visualize_ch_response/';
subject_id='AMC026';
d= dir([data_path,strcat('/**/',subject_id,'*_crunched.mat')]);
fprintf(' %d .mat files were found \n', length(d));
i
% find the node information
d_nmerge= dir([data_path,'/**/*nmerge.txt']);
fprintf(' %d nmerge files were found \n', length(d_nmerge));
% load node info
i=1 
[node_cell,node_table]=generate_node_table_eh(strcat(d_nmerge(i).folder,'/',d_nmerge(i).name));
% 
plot_width=.9/length(d);
plot_length=.87;
%% 
%figure;
%set(gcf,'OuterPosition',[-1615 441 694 1337]);
%plot_width=.9/length(d);
%plot_length=.87;
find_index= @(input_cell) find(cell2mat(cellfun(@(x) ~isempty(x), input_cell,'UniformOutput',false)));
colors = cbrewer('seq','YlGnBu',  15);
colors=flipud(colors);

for i=1:length(d)
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
    sentences=[data{sentence_trial_index}];
    sentence_gamma_band_envelope=cellfun(@(x) x(1:8),{sentences.signal_hilbert_zs_downsample_parsed},'UniformOutput',false);
    % 
    example_sentence={sentences(:).trial_string};
    example_sentence=cellfun(@(x) x(2:end),example_sentence,'UniformOutput',false);
    % find trial location in node cell 
    example_sentence_locations=cellfun(@(x) (regexpi([node_cell{:,1}],x)),example_sentence,'UniformOutput',false);
    example_sentence_locations=cell2mat(cellfun(@(x) find_index(x), example_sentence_locations, 'UniformOutput',false));
    sentences_node_closing=node_cell(example_sentence_locations,3);
    sentences_open_nodes=node_cell(example_sentence_locations,2);
    % add the 
    % creat a cell with wordposition(row)*time in trial(column) structure
    sentence_gamma_band_envelope=[sentence_gamma_band_envelope{:,:}];
    % convert the gamma into a tensor of shape : channel*time in trial*
    % trial number
    sentence_gamma_band_envelope_tensor=cell2mat(permute(sentence_gamma_band_envelope,[3,1,2]));
    %append to langauge_channel*trial*words positions
    sentence_gamma_band_envelope_tensor=sentence_gamma_band_envelope_tensor(find(language_electrode),:,:);
    %  words trials 
    words_trial_index=~cellfun(@isempty,cellfun(@ (x) strfind(x,'W'),info.word_type,'UniformOutput',false));
    words=[data{words_trial_index}];
    words_gamma_band_ave_envelope=cellfun(@(x) x(1:8),{words.signal_hilbert_zs_downsample_parsed},'UniformOutput',false);
    words_gamma_band_ave_envelope=[words_gamma_band_ave_envelope{:,:}];
    % convert the gamma into a tensor of shape : channel*trial*words
    words_gamma_band_ave_envelope_tensor=cell2mat(permute(words_gamma_band_ave_envelope,[3,2,1]));

    
    
    fprintf('added %s from %s \n',d(i).name, strcat(d(i).folder,'/',d(i).name));
    
    
    for k=1:2%size(sentence_gamma_band_envelope_tensor,1)
        figure('position',[-1856 442 717 1309]);
        axes('position',[.2,.2,.6,.6])
        anchor_points=0.0:135:size(sentence_gamma_band_envelope_tensor,2)-1;
        jump_val=floor(max(max(sentence_gamma_band_envelope_tensor(k,:,:)))-min(min(sentence_gamma_band_envelope_tensor(k,:,:))));
        for kk=1:size(sentence_gamma_band_envelope_tensor,3)
            plot(0*sentence_gamma_band_envelope_tensor(k,:,kk)+(kk-1)*jump_val,'k-','linewidth',1);
            hold on
            plot(sentence_gamma_band_envelope_tensor(k,:,kk)+(kk-1)*jump_val,'-','linewidth',2,'color',colors(kk,:));
            %plot(words_gamma_band_ave_envelope_tensor(k,:,kk)+(kk-1)*jump_val,'-','linewidth',.5,'color',colors(kk,:));
            plot([0,0,0]-5,[-2,0,2]'+(kk-1)*jump_val,'-','color',[0,0,0])
            text([0,0,0]-10,double([-2,0,2]+(kk-1)*jump_val),strsplit(num2str([-2,0,2]),' '),...
                    'FontSize',10,'color',[0,0,0],'VerticalAlignment','middle','HorizontalAlignment','right','FontWeight','bold');   
        
            
            ylim=get(gca,'ylim');
            a=anchor_points*0+min(sentence_gamma_band_envelope_tensor(k,:,kk))+(kk-1)*jump_val;
                text(anchor_points+5,double(a),strsplit(example_sentence{1,kk},' '),...
                    'FontSize',8,'color',colors(kk,:),'VerticalAlignment','top','HorizontalAlignment','left','FontWeight','bold');   
        end 
        set(gca, 'ydir', 'normal','box','off','xtick',...
            0.0:135:size(sentence_gamma_band_envelope_tensor,2)-1,'XTickLabel',...
            [1:size(sentence_gamma_band_envelope_tensor,3)],'xlim',[-10,size(sentence_gamma_band_envelope_tensor,2)+10]);
        plot(repmat(transpose(0.0:135:size(sentence_gamma_band_envelope_tensor,2)-1),1,2)',...
            repmat(get(gca,'ylim'),size(sentence_gamma_band_envelope_tensor,3),1)','--','color',[.5,.5,.5])
        set(gca,'ycolor',[1 1 1 ])
        xlabel('word position','FontWeight','bold');
       text(-25,0,'z-score','FontSize',10,'color',[0,0,0],'VerticalAlignment','middle','HorizontalAlignment','right','FontWeight','bold');   
        title(strcat(info.subject,'-session ',num2str(i),'- electrode ',num2str(language_electrode_num(k))))
      if ~exist(strcat(analysis_path,info.subject))
        mkdir(strcat(analysis_path,info.subject))
      end 
      print(gcf, '-depsc', strcat(analysis_path,info.subject,'/',info.subject,'_channel_',num2str(language_electrode_num(k)),'_visualization','.eps')); 
      close(gcf)
    end 
    
  
        
end 



%% 
open_nodes=cell2mat(cellfun(@(x) x(1:8), node_table.open_nodes,'UniformOutput',false));
unique_open_node_patterns=unique(open_nodes,'rows');
pattern_cell={};
for nn=1:8
    pattern=unique(unique_open_node_patterns(:,nn));
    pattern_cell=[pattern_cell,pattern];
end 
%% there are 31 unique patterns, make a table of each sentence and for each word location what pattern it is member of 
open_nodes_cell=cellfun(@(x) x(1:8), node_table.open_nodes,'UniformOutput',false);
all_pattern_locations={};
for p=1:size(unique_open_node_patterns,1)
    unique_open_node_pattern=unique_open_node_patterns(p,:)
    pattern_locations={};
    for i=1:8
        a=unique_open_node_pattern(1:i);
        pattern_locations=[pattern_locations,
        find(~cell2mat(cellfun(@isempty,cellfun(@(x) strfind(num2str(x),num2str(a)),open_nodes_cell,'UniformOutput',false),'UniformOutput',false)))];
    end
    all_pattern_locations=[all_pattern_locations;[unique_open_node_pattern,transpose(pattern_locations)]];
end

%% 
%% there are 31 unique patterns, make a table of each sentence and for each word location what pattern it is member of 
open_nodes_cell=cellfun(@(x) x(1:8), node_table.open_nodes,'UniformOutput',false);
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
%% for each electrode and all pattern condition , find relevant sentences and add them to the so you have num_electrode*pattern_condition*time tensor in the end

find_index= @(input_cell) find(cell2mat(cellfun(@(x) ~isempty(x), input_cell,'UniformOutput',false)));

for i=1:length(d)
    fprintf('adding %s from %s \n',d(i).name, strcat(d(i).folder,'/',d(i).name));
    subj=load(strcat(d(i).folder,'/',d(i).name));
    subj_id=fieldnames(subj);
    subj=subj.(subj_id{1});
    data=subj.data;
    info=subj.info;
    language_electrode=info.language_responsive_electrodes;
    if ~exist('electrode_pattern_section_tensor')
        electrode_pattern_section_tensor=nan*ones(sum(language_electrode),size(unique_open_node_patterns,1),8*135);
    end 
    language_electrode_num=find(language_electrode);
    % step 1: extract electrodes with siginificant language response
    sentence_trial_index=~cellfun(@isempty,cellfun(@ (x) strfind(x,'S'),info.word_type,'UniformOutput',false));
    % sentence
    sentences=[data{sentence_trial_index}];
    sentence_gamma_band_envelope=cellfun(@(x) x(1:5),{sentences.signal_hilbert_zs_downsample_parsed},'UniformOutput',false);
    example_sentence={sentences(:).trial_string};
    example_sentence=cellfun(@(x) x(2:end),example_sentence,'UniformOutput',false);
    % find trial location in node cell 
    example_sentence_locations=cellfun(@(x) (regexpi([node_cell{:,1}],x)),example_sentence,'UniformOutput',false);
    example_sentence_locations=cell2mat(cellfun(@(x) find_index(x), example_sentence_locations, 'UniformOutput',false));
    
    % creat a cell with wordposition(row)*time in trial(column) structure
    sentence_gamma_band_envelope=[sentence_gamma_band_envelope{:,:}];
    sentence_lang_elec_gamma_band_envelope=cellfun(@(x) x(language_electrode_num,:),sentence_gamma_band_envelope,'UniformOutput',false);
    time_vector=reshape(1:8*135,[],8);
    for k=1:length(example_sentence_locations)
        
        sentence_elec_gamma_bands=sentence_lang_elec_gamma_band_envelope(:,k);
        sentence_pattern_membership=all_pattern_locations(example_sentence_locations(k),2:end);
        total_iterations=sum(cell2mat(cellfun(@length,sentence_pattern_membership,'UniformOutput',false)));
        fprintf('doing sentence %d total iterations %d \n',example_sentence_locations(k),total_iterations);
        for kk=1:length(sentence_lang_elec_gamma_band_envelope,1)
            section_membership=sentence_pattern_membership{kk};
            fprintf('doing section %d \n',kk);
            for kkk=1:length(section_membership)
                electrode_pattern_section_tensor(:,section_membership(kkk),time_vector(:,kk))=sentence_lang_elec_gamma_band_envelope(kk,k);
                % assign the channel response to the secion
            end 
        end
        
    end 
   
    
    % convert the gamma into a tensor of shape : channel*time in trial*
    % trial number
    %sentence_lang_elec_gamma_band_envelope_tensor=cell2mat(permute(sentence_lang_elec_gamma_band_envelope,[3,1,2]));

    fprintf('added %s from %s \n',d(i).name, strcat(d(i).folder,'/',d(i).name));
    
end
%% num_electrode*pattern_condition*time
a=nanmean(electrode_pattern_section_tensor,2);


%% 
i=1;
strings = textscan(extractFileText(strcat(d(i).folder,'/',d(i).name)),'%s','Delimiter','\n')';
strings=[strings{:}];
strings=cellfun(@(x) erase(x,'(. .)'),strings,'UniformOutput',false);
a=cellfun(@(x) regexpi(x,'\)'), strings,'UniformOutput',false);
b=cellfun(@(x) regexpi(x,'\('), strings,'UniformOutput',false);
closings=a{1};
openings=b{1};
closings=[closings;-sign(closings)];
openings=[openings;sign(openings)];
[~,sort_id]=sort([openings(1,:),closings(1,:)]);
open_close=[openings(2,:),closings(2,:)];
open_close=open_close(sort_id);
a=make_hierarchical_tree(open_close);


b1=mean(a{1},1);
domain_val=[];
domain=[];
ranges=transpose(a{2});
for i=1:3
    domain(i,:)=(b1<ranges(i,2) & b1>ranges(i,1));
    domain_val(i)=mean(b1(b1<ranges(i,2) & b1>ranges(i,1)),2);
end
domain_remove=~sum(domain,1);
new_b=b1(find(domain_remove));
c1=sort([domain_val,new_b]);

ranges=transpose(a{3});

domain=[];
for i=1:size(ranges,1)
    domain(i,:)=(c1<ranges(i,2) & c1>ranges(i,1));
    
end




