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
      print(gcf, '-dtiffn', strcat(analysis_path,info.subject,'/',info.subject,'_channel_',num2str(language_electrode_num(k)),'_visualization')); 
      print(gcf, '-dpng', strcat(analysis_path,info.subject,'/',info.subject,'_channel_',num2str(language_electrode_num(k)),'_visualization')); 
      close(gcf)
    end        
end 



