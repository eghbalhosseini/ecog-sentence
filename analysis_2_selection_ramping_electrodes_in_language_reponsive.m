% selecting langauge responsive electrodes and add them to the subject
% info. 
% tested for subject 1 and validated with Terry's data . 
%% step 0: prepare the data 
clear all 
close all 
home 
%% 
data_path='/Users/eghbalhosseiniasl1/MyData/ecog-sentence';
d= dir([data_path,'/**/AMC026*_crunched.mat']);
fprintf(' %d .mat files were found \n', length(d));
gamma_band_index=4;
p_threshold=0.01;
%% 
electrode_with_langauge_accross_sessions=[];
for i=2:2:length(d)
    fprintf('adding %s from %s \n',d(i).name, strcat(d(i).folder,'/',d(i).name));
    subj=load(strcat(d(i).folder,'/',d(i).name));
    subj_id=fieldnames(subj);
    subj=subj.(subj_id{1});
    data=subj.data;
    info=subj.info;
    language_electrode=info.language_responsive_electrodes;

% step 1: extract electrodes with siginificant language response
    
    sentence_trial_index=~cellfun(@isempty,cellfun(@ (x) strfind(x,'S'),info.word_type,'UniformOutput',false));
    % sentence 
    sentences=[data{sentence_trial_index}];
    sentence_gamma_band_ave_envelope=cellfun(@(x) x(1:8,gamma_band_index),{sentences.signal_ave_envelope_downsample_parsed},'UniformOutput',false);
    % creat a cell with wordposition(row)*trial(column) structure
    sentence_gamma_band_ave_envelope=[sentence_gamma_band_ave_envelope{:,:}];
    % convert the gamma into a tensor of shape : channel*trial*words
    sentence_gamma_band_ave_envelope_tensor=cell2mat(permute(sentence_gamma_band_ave_envelope,[3,2,1]));
    %append to langauge_channel*trial*words positions
    electrodes_with_language_response=sentence_gamma_band_ave_envelope_tensor(find(language_electrode),:,:);
    electrode_with_langauge_accross_sessions=cat(2,electrode_with_langauge_accross_sessions,electrodes_with_language_response);
    % 
    fprintf('added %s from %s \n',d(i).name, strcat(d(i).folder,'/',d(i).name));
end   
%% find ramping electrodes 
electrode_num=find(language_electrode);
sentence_all=[];
for i=1:length(electrode_num)
    channel_response=(double(squeeze(electrode_with_langauge_accross_sessions(i,:,:))));
    
    word_position=(repmat(1:size(channel_response,2),[size(channel_response,1),1]));
    [r_sentence,p_sent]=corr(channel_response(:),word_position(:),'type','Spearman','rows','complete');
    sentence_all=[sentence_all;...
        [r_sentence,p_sent]];
end 

ramp_electrodes=electrode_num(find(sentence_all(:,2) < 0.001 & sentence_all(:,1) > 0));
ramp_electrodes_location=(sentence_all(:,1)>0 & sentence_all(:,2)<0.01);
% 
electrode_with_ramp_across_sessions=electrode_with_langauge_accross_sessions(ramp_electrodes_location,:,:);
channel_ramp=zeros(size(language_electrode,1),1);
channel_ramp(ramp_electrodes)=1;
%% plot electrodes
close all;
electrode_num=find(language_electrode);
colors=cbrewer('qual','Set1',10);
num_rows=4;
num_columns=2;
total_plots=num_rows*num_columns;
for i=1:length(electrode_num)
    channel_response=double(squeeze(electrode_with_langauge_accross_sessions(i,:,:)));
    word_position=repmat(1:size(channel_response,2),[size(channel_response,1),1]);
    perturbed_word_position=word_position+.2*rand(size(word_position,1),size(word_position,2))-.25;
  
    figure(fix((i-1)/total_plots)+1);
    set(gcf,'Position',[-882 449 877 1325]);
    sub_title=sprintf('subj: %s,\n electrodes: %d ', info.subject ,electrode_num(i));
    a=subplot(num_rows,num_columns,i-total_plots*fix((i-1)/total_plots));
     ah=violinPlot(channel_response, 'histOri', 'right', 'widthDiv', [2 2], 'showMM', 0, ...
    'color',  mat2cell(colors(1, : ), 1),'histOpt',0);
    hold on 
    h=scatter(perturbed_word_position(:),channel_response(:),2);
    set(h,'CData',colors(6,:),'MarkerFacecolor',colors(1,:),'MarkerFaceAlpha',.3)
   
    e=plot(mean(word_position,1),mean(channel_response,1),...
        '-s','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',2,'color','k');
    xlim([min(min(word_position))-.5,max(max(word_position))+.5])
    set(a,'XTick',[1:8])
    e.Color=[0,0,0];
    title(sub_title);
    set(a,'color','none','layer','top')

end  

%% 
% sometimes, axes are plotted in a dark grey thats not exactly black (which
% % I find annoying). Make sure this doesnt happen.
% figHandles = get(groot, 'Children');
% for b= 1:length(figHandles)
%     axes = findobj(figHandles(b), 'type', 'axes');
%     for a = 1:length(axes),
%         
%         if axes(a).YColor &amp;amp;amp;amp;lt; [1 1 1],
%             axes(a).YColor = [0 0 0];
%         end
%         if axes(a).XColor &amp;amp;amp;amp;lt; [1 1 1],
%             axes(a).XColor = [0 0 0];
%         end
%     end
%     
%     print(figHandles(b), '-dpdf', strcat(subj_id,'_ramp_effect','.pdf'));
% end
 
% when you plotted several subplots but want them to have shared axes, use
% suplabel

% save to pdf
% see also export_fig from the file exchange

