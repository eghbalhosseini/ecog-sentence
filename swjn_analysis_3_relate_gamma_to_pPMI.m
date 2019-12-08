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
analysis_path='/Users/eghbalhosseiniasl1/MyData/ECoG-sentence/analysis/swjn_analysis_4_relate_gamma_to_pPMI/';
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
%% find the pMI information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d_pmi= dir([data_path,'/**/*pmi_sentences.csv']);
fprintf(' %d pmi files were found \n', length(d_pmi));
i=1;
[pmi_cell]=generate_pmi_table_eh(strcat(d_pmi(i).folder,'/',d_pmi(i).name));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  STEP 1 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
find_index= @(input_cell) find(cell2mat(cellfun(@(x) ~isempty(x), input_cell,'UniformOutput',false)));
all_sentence_pattern=[];
all_pmi_pattern=[];
all_sentence_pattern_id=[];

sentence_electrode_with_langauge_accross_sessions=[];
words_electrode_with_langauge_accross_sessions=[];


session_sentence_hilbert_band_envelope_lang_elec_tensor=[];
session_sentence_hilbert_zs_band_envelope_lang_elec_tensor=[];
session_nonwords_hilbert_band_envelope_lang_elec_tensor=[];
session_words_hilbert_band_envelope_lang_elec_tensor=[];
session_jabberwocky_hilbert_band_envelope_lang_elec_tensor=[];

for i=1:length(d)
    subj=load(strcat(d(i).folder,'/',d(i).name));
    subj_id=fieldnames(subj);
    subj=subj.(subj_id{1});
    data=subj.data;
    info=subj.info;
    language_electrode=info.ramp_electrodes_hilbert_odd;
    language_electrode_num=find(language_electrode);
    
    
    
    % step 1: extract electrodes with siginificant language response
    sentence_trial_index=~cellfun(@isempty,cellfun(@ (x) strfind(x,'S'),info.word_type,'UniformOutput',false));
    %%%%%%%%%%%%%%% sentence
    sentences=[data{sentence_trial_index}];% 
    word_length=sentences(1).signal_range_downsample(1,2)-sentences(1).signal_range_downsample(1,1)+1;
    % creat a cell with wordposition(row)*time in trial(column) structure
    %  
    hilbert_band_envelope=cellfun(@(x) x(1:8),{sentences.signal_hilbert_zs_downsample_parsed},'UniformOutput',false);
    hilbert_band_envelope=[hilbert_band_envelope{:,:}];
    hilbert_band_envelope=cellfun(@transpose,hilbert_band_envelope,'UniformOutput',false);
    % make a words positions*channel* trial tensor 
    hilbert_band_envelope_tensor=cell2mat(permute(hilbert_band_envelope,[1,3,2]));
    %append to langauge_channel*trial*words positions
    hilbert_band_envelope_tensor=permute(hilbert_band_envelope_tensor,[2,3,1]);
    hilbert_band_envelope_lang_elec_tensor=hilbert_band_envelope_tensor(find(language_electrode),:,:);

    session_sentence_hilbert_band_envelope_lang_elec_tensor=cat(2,session_sentence_hilbert_band_envelope_lang_elec_tensor,hilbert_band_envelope_lang_elec_tensor);
    % find sentence locations 
    example_sentence=cellfun(@(x) x(2:end),{sentences(:).trial_string},'UniformOutput',false);
    example_sentence_locations=cellfun(@(x) (regexpi([pmi_cell(:,1)],x)),example_sentence,'UniformOutput',false);
    example_sentence_locations=cell2mat(cellfun(@(x) find_index(x), example_sentence_locations, 'UniformOutput',false));
    all_sentence_pattern=[all_sentence_pattern;
        pmi_cell(example_sentence_locations,1)];
    all_sentence_pattern_id=[all_sentence_pattern_id;example_sentence_locations'];
    all_pmi_pattern=[all_pmi_pattern;
        pmi_cell(example_sentence_locations,2)];
    %%%%%%%%%%%%%%%%% words
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
    hilbert_band_envelope_lang_elec_tensor=hilbert_band_envelope_tensor(find(language_electrode),:,:);
    session_words_hilbert_band_envelope_lang_elec_tensor=cat(2,session_words_hilbert_band_envelope_lang_elec_tensor,hilbert_band_envelope_lang_elec_tensor);
    %%%%%%%%%%%%%%%%% Nonwords
    nonwords_trial_index=~cellfun(@isempty,cellfun(@ (x) strfind(x,'N'),info.word_type,'UniformOutput',false));
    % nonwords
    nonwords=[data{nonwords_trial_index}];

    
    hilbert_band_ave_envelope_pre_trial=cellfun(@(x) reshape(x,size(x,1),1,size(x,2)),{nonwords.signal_ave_pre_trial_hilbert_downsample},'UniformOutput',false);
    hilbert_band_ave_envelope_pre_trial=[hilbert_band_ave_envelope_pre_trial{:,:}];
    hilbert_band_envelope=cellfun(@(x) x(1:8),{nonwords.signal_hilbert_zs_downsample_parsed},'UniformOutput',false);
    % creat a cell with wordposition(row)*time in trial(column) structure
    hilbert_band_envelope=[hilbert_band_envelope{:,:}];
    hilbert_band_envelope=cellfun(@transpose,hilbert_band_envelope,'UniformOutput',false);
    % make a nonwords positions*channel* trial tensor
    hilbert_band_envelope_tensor=cell2mat(permute(hilbert_band_envelope,[1,3,2]));
    %append to langauge_channel*trial*nonwords positions
    hilbert_band_envelope_tensor=permute(hilbert_band_envelope_tensor,[2,3,1]);
    hilbert_band_envelope_lang_elec_tensor=hilbert_band_envelope_tensor(find(language_electrode),:,:);
    session_nonwords_hilbert_band_envelope_lang_elec_tensor=cat(2,session_nonwords_hilbert_band_envelope_lang_elec_tensor,hilbert_band_envelope_lang_elec_tensor);
    %%%%%%%%%%%%%%%%%Jabberwocky
    jabberwocky_trial_index=~cellfun(@isempty,cellfun(@ (x) strfind(x,'J'),info.word_type,'UniformOutput',false));
    jabberwocky=[data{jabberwocky_trial_index}];

    hilbert_band_ave_envelope_pre_trial=cellfun(@(x) reshape(x,size(x,1),1,size(x,2)),{jabberwocky.signal_ave_pre_trial_hilbert_downsample},'UniformOutput',false);
    hilbert_band_ave_envelope_pre_trial=[hilbert_band_ave_envelope_pre_trial{:,:}];
    hilbert_band_envelope=cellfun(@(x) x(1:8),{jabberwocky.signal_hilbert_zs_downsample_parsed},'UniformOutput',false);
    % creat a cell with wordposition(row)*time in trial(column) structure
    hilbert_band_envelope=[hilbert_band_envelope{:,:}];
    hilbert_band_envelope=cellfun(@transpose,hilbert_band_envelope,'UniformOutput',false);
    % make a nonwords positions*channel* trial tensor
    hilbert_band_envelope_tensor=cell2mat(permute(hilbert_band_envelope,[1,3,2]));
    %append to langauge_channel*trial*nonwords positions
    hilbert_band_envelope_tensor=permute(hilbert_band_envelope_tensor,[2,3,1]);
    hilbert_band_envelope_lang_elec_tensor=hilbert_band_envelope_tensor(find(language_electrode),:,:);
    session_jabberwocky_hilbert_band_envelope_lang_elec_tensor=cat(2,session_jabberwocky_hilbert_band_envelope_lang_elec_tensor,hilbert_band_envelope_lang_elec_tensor);
    
    fprintf('added %s from %s \n',d(i).name, strcat(d(i).folder,'/',d(i).name));
    clear data subj
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  find next word pmi  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
all_next_word_pmi=[];
for kk=1:size(all_pmi_pattern,1)
    pmi=all_pmi_pattern{kk};
    next_word_pmi=[];
    for kkk=1:size(pmi,1)-1
        next_word_pmi(kkk)=pmi(kkk,kkk+1);
    end 
    all_next_word_pmi=[all_next_word_pmi;next_word_pmi];
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  regression for word position and pmi 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all 
look_out_window=100; %200 ms
offset=30; % 100 ms
word_range=[1:8];
num_rows=5;
num_columns=2;
total_plots=num_rows*num_columns;
p=0;
for i=1:length(language_electrode_num)
    elec_response_to_next_word_pattern=[];
    electrode_response=squeeze(session_sentence_hilbert_band_envelope_lang_elec_tensor(i,:,:));
    electrode_response_sentence_cell=mat2cell(electrode_response,ones(1,size(electrode_response,1)),size(electrode_response,2)/8*ones(1,8));
    a=cell2mat(cellfun(@(x)  nanmean(x(offset+[1:look_out_window])),electrode_response_sentence_cell,'UniformOutput',false));
    a_diff=diff(a,1,2);
    y=double(a_diff(:));
    x=double(all_next_word_pmi(:));
    lm = fitlm(x,y);
    a_fit=lm.Fitted;
    
   
    
    f=figure(fix((i-1)/total_plots)+1);
    %set(f,'position',[1554 45 719 1281])
    sub_title=sprintf('subj: %s,\n electrodes: %d ', info.subject ,language_electrode_num(i));
    ax=subplot(num_rows,num_columns,i-total_plots*fix((i-1)/total_plots));
   
    bl=scatter(x,y,25,[1,.2,.2],'filled');
    bl.MarkerEdgeColor=[1,1,1];bl.MarkerEdgeAlpha=.5;bl.MarkerFaceAlpha=.7;
    hold on
    
    plot(x,a_fit,'k-','LineWidth',2)
    h1=plot(get(gca,'xlim'),0*get(gca,'ylim'),'k--');hAnnotation = get(h1,'Annotation');hLegendEntry = get(hAnnotation,'LegendInformation');set(hLegendEntry,'IconDisplayStyle','off');
    ax.XAxis.Visible = 'on'; 
    y_quantile=quantile(y,50);
    set(ax,'ydir', 'normal','box','off','ylim',[y_quantile(1),y_quantile(end)]);
    text(max(get(gca,'xlim')),min(get(gca,'ylim')),{strcat('slope: ',num2str(lm.Coefficients.Estimate(2))),strcat('R^2: ',num2str(lm.Rsquared.Ordinary(1))),...
        strcat('p_{value}: ',num2str(lm.Coefficients.pValue(2)))},...
        'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10,'fontweight','bold')
    title(sub_title);
    
    if ~mod(i,total_plots) | i==length(language_electrode_num)
        %legend('show','Location','northeastoutside')
        p=p+1;
        if ~exist(strcat(analysis_path,info.subject))
            mkdir(strcat(analysis_path,info.subject))
        end 
        xlabel('mutual info');
        window=sprintf('%d_%dms',floor(min(1000/info.downsample_sampling_rate*(offset+[1:look_out_window]))),...
            floor(max(1000/info.downsample_sampling_rate*(offset+[1:look_out_window]))));
        set(gcf,'PaperPosition',[.25 .25 8 6])
        print(gcf,'-fillpage', '-dpdf', strcat(analysis_path,info.subject,'/',info.subject,'_mutual_info',num2str(p),'_window_',window,'.pdf')); 
    end 
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  regression for word position and pmi 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


close all 
look_out_window=100; %200 ms
offset=30; % 100 ms
word_range=[1:8];
num_rows=5;
num_columns=2;
total_plots=num_rows*num_columns;
p=0;
for i=1:length(language_electrode_num)
    elec_response_to_next_word_pattern=[];
    electrode_response=squeeze(session_sentence_hilbert_band_envelope_lang_elec_tensor(i,:,:));
    electrode_response_sentence_cell=mat2cell(electrode_response,ones(1,size(electrode_response,1)),size(electrode_response,2)/8*ones(1,8));
    a=cell2mat(cellfun(@(x)  nanmean(x),electrode_response_sentence_cell,'UniformOutput',false));
    a_diff=diff(a,1,2);
    y=double(a_diff(:));
    x=double(all_next_word_pmi(:));
    lm = fitlm(x,y);
    a_fit=lm.Fitted;
    
   
    
    f=figure(fix((i-1)/total_plots)+1);
    %set(f,'position',[1554 45 719 1281])
    sub_title=sprintf('subj: %s,\n electrodes: %d ', info.subject ,language_electrode_num(i));
    ax=subplot(num_rows,num_columns,i-total_plots*fix((i-1)/total_plots));
   
    bl=scatter(x,y,25,[1,.2,.2],'filled');
    bl.MarkerEdgeColor=[1,1,1];bl.MarkerEdgeAlpha=.5;bl.MarkerFaceAlpha=.7;
    hold on
    
    plot(x,a_fit,'k-','LineWidth',2)
    h1=plot(get(gca,'xlim'),0*get(gca,'ylim'),'k--');hAnnotation = get(h1,'Annotation');hLegendEntry = get(hAnnotation,'LegendInformation');set(hLegendEntry,'IconDisplayStyle','off');
    ax.XAxis.Visible = 'on'; 
    y_quantile=quantile(y,50);
    set(ax,'ydir', 'normal','box','off','ylim',[y_quantile(1),y_quantile(end)]);
    text(max(get(gca,'xlim')),min(get(gca,'ylim')),{strcat('slope: ',num2str(lm.Coefficients.Estimate(2))),strcat('R^2: ',num2str(lm.Rsquared.Ordinary(1))),...
        strcat('p_{value}: ',num2str(lm.Coefficients.pValue(2)))},...
        'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10,'fontweight','bold')
    title(sub_title);
    
    if ~mod(i,total_plots) | i==length(language_electrode_num)
        %legend('show','Location','northeastoutside')
        p=p+1;
        if ~exist(strcat(analysis_path,info.subject))
            mkdir(strcat(analysis_path,info.subject))
        end 
        xlabel('mutual info');
        window=sprintf('%d_%dms',floor(min(1000/info.downsample_sampling_rate*(offset+[1:look_out_window]))),...
            floor(max(1000/info.downsample_sampling_rate*(offset+[1:look_out_window]))));
        set(gcf,'PaperPosition',[.25 .25 8 6])
        print(gcf,'-fillpage', '-dpdf', strcat(analysis_path,info.subject,'/',info.subject,'_mutual_info','_window_full_',num2str(p),'.pdf')); 
    end 
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  create a schematic of analysis pipline  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all 
look_out_window=100; %200 ms
offset=30; % 100 ms
word_range=[1:8];
num_rows=5;
num_columns=2;
total_plots=num_rows*num_columns;
p=0;
for i=1%:%length(language_electrode_num)
    elec_response_to_next_word_pattern=[];
    electrode_response=squeeze(session_sentence_hilbert_band_envelope_lang_elec_tensor(i,:,:));
    electrode_response_sentence_cell=mat2cell(electrode_response,ones(1,size(electrode_response,1)),size(electrode_response,2)/8*ones(1,8));
    electrode_1st_sent_cell=[electrode_response_sentence_cell(1,:)];
    a=cell2mat(cellfun(@(x)  nanmean(x),electrode_1st_sent_cell,'UniformOutput',false));
    b=ones(size(a,2));
    [row,col]=find(b);
    a_all_diffs=arrayfun(@(x,y) a(y)-a(x),row,col);
    a_all_diffs_upper=triu(reshape(a_all_diffs,size(a,2),[]),1);
    wordpos_all_diffs=arrayfun(@(x,y) y-x,row,col);
    wordpos_all_diffs_upper=triu(reshape(wordpos_all_diffs,size(a,2),[]),1);
    [row_p,col_p]=find(wordpos_all_diffs_upper);
    
    f=figure;
    set(f,'position',[-1341 334 960 960])
    ax=axes('position',[.05,.1,.25,.25])
    
    hold on 
    ax.XLim=[0,9];
    ax.YLim=[0,9];
    arrayfun(@(x,y,z) text(x,y,num2str(wordpos_all_diffs_upper(x,y)),'FontSize',12),row_p,col_p);
   
    ax.XAxis.Direction='reverse';
    ax.FontSize=16
    
    text(ax.XLim(2),ax.XLim(2)+1,'A','HorizontalAlignment','right','FontSize',20,'FontWeight','bold')
    daspect([1,1,1])
    ax.XTick=[1:8];ax.XTickLabelRotation=90;
    ax.XTickLabel=fliplr(strsplit(all_sentence_pattern{1}));
    ax.YTick=[1:8];ax.YTickLabelRotation=0;
    ax.YTickLabel=fliplr(strsplit(all_sentence_pattern{1}));
    ax.Title.String='Word Position';
    ax=axes('position',[.30,.1,.25,.25])
    hold on 
    ax.XLim=[0,9];
    ax.YLim=[0,9];
    pmi_sent=pmi_cell{all_sentence_pattern_id(1),2};
    arrayfun(@(x,y,z) text(y,x,num2str(pmi_sent(x,y),'%1.1f'),'FontSize',12),row_p,col_p);
    ax.YAxis.Direction='reverse';
    ax.XAxis.Direction='normal';
    daspect([1,1,1])
    ax.XTick=[1:8];ax.XTickLabelRotation=90;
    ax.XTickLabel=(strsplit(all_sentence_pattern{1}));
    ax.YTick=[1:8];ax.YTickLabelRotation=0;
    ax.YTickLabel=(strsplit(all_sentence_pattern{1}));
    ax.Title.String='Pointwise Mutual Info';
    text(ax.XLim(1),ax.XLim(1)-1,'B','HorizontalAlignment','right','FontSize',20,'FontWeight','bold')
    % 
    ax.YTick=[];
    ax.FontSize=16;
    ax=axes('position',[.55,.1,.25,.25]);
    hold on 
    ax.XLim=[0,9];
    ax.YLim=[0,9];
   
    arrayfun(@(x,y,z) text(y,x,num2str(a_all_diffs_upper(x,y),'%1.1f'),'FontSize',12,'horizontalalignment','center'),row_p,col_p);
    ax.YAxis.Direction='reverse';
    ax.XAxis.Direction='normal';
   daspect([1,1,1])
    fliplr(strsplit(all_sentence_pattern{1}))
    ax.XTick=[1:8];ax.XTickLabelRotation=90;
    ax.XTickLabel=(strsplit(all_sentence_pattern{1}));
    ax.YTick=[1:8];ax.YTickLabelRotation=0;
    ax.YTickLabel=(strsplit(all_sentence_pattern{1}));
     ax.Title.String={'Change in',' gamma power'};
     ax.YTick=[];
     ax.FontSize=16;
    text(ax.XLim(1),ax.XLim(1)-1,'C','HorizontalAlignment','right','FontSize',20,'FontWeight','bold')
    ax=axes('position',[.80,.1,.07,.25]);
    p=plot(a,[1:8],'Marker','.','color','k');
    p.MarkerSize=20;
    p.LineWidth=2
    ax.YLim=[0,9];
    ax.XLim=1.2*[min(a),max(a)];
    ax.XAxis.Direction='normal';
    ax.YAxis.Direction='reverse';
    ax.YAxis.Visible='off'
    ax.Box='off'
    ax.Title.String={'gamma','power'};
ax.FontSize=16;
        if ~exist(strcat(analysis_path,info.subject))
            mkdir(strcat(analysis_path,info.subject))
        end 
        %set(gcf,'PaperPosition',[.25 .25 8 6])
        print(gcf, '-djpeg', strcat(analysis_path,info.subject,'/',info.subject,'_representitive_fig')); 
    
end

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  create a schematic of analysis pipline for presentation   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all 
look_out_window=100; %200 ms
offset=30; % 100 ms
word_range=[1:8];
num_rows=5;
num_columns=2;
total_plots=num_rows*num_columns;
p=0;
for i=1%:%length(language_electrode_num)
    elec_response_to_next_word_pattern=[];
    electrode_response=squeeze(session_sentence_hilbert_band_envelope_lang_elec_tensor(i,:,:));
    electrode_response_sentence_cell=mat2cell(electrode_response,ones(1,size(electrode_response,1)),size(electrode_response,2)/8*ones(1,8));
    electrode_1st_sent_cell=[electrode_response_sentence_cell(1,:)];
    a=cell2mat(cellfun(@(x)  nanmean(x),electrode_1st_sent_cell,'UniformOutput',false));
    b=ones(size(a,2));
    [row,col]=find(b);
    a_all_diffs=arrayfun(@(x,y) a(y)-a(x),row,col);
    a_all_diffs_upper=triu(reshape(a_all_diffs,size(a,2),[]),1);
    wordpos_all_diffs=arrayfun(@(x,y) y-x,row,col);
    wordpos_all_diffs_upper=triu(reshape(wordpos_all_diffs,size(a,2),[]),1);
    [row_p,col_p]=find(wordpos_all_diffs_upper);
    
    close all 
f=figure;
aspect_ration=9.32./4.13;
y=500;
set(f,'position',[591 455 aspect_ration*y y]);
ax=axes('position',[.05,.15,.22,.22*aspect_ration])
    
    hold on 
    ax.XLim=[0,9];
    ax.YLim=[0,9];
    arrayfun(@(x,y,z) text(x,y,num2str(wordpos_all_diffs_upper(x,y)),'FontSize',12),row_p,col_p);
   
    ax.XAxis.Direction='reverse';
    ax.FontSize=14;
    ax.FontWeight='bold';
    
    %text(ax.XLim(2),ax.XLim(2)+1,'A','HorizontalAlignment','right','FontSize',20,'FontWeight','bold')
    daspect([1,1,1])
    ax.XTick=[1:8];ax.XTickLabelRotation=90;
    ax.XTickLabel=fliplr(strsplit(all_sentence_pattern{1}));
    ax.YTick=[1:8];ax.YTickLabelRotation=0;
    ax.YTickLabel=fliplr(strsplit(all_sentence_pattern{1}));
    ax.Title.String='Word Position (wp)';
    ax=axes('position',[.33,.15,.22,.22*aspect_ration])
    hold on 
    ax.XLim=[0,9];
    ax.YLim=[0,9];
    pmi_sent=pmi_cell{all_sentence_pattern_id(1),2};
    arrayfun(@(x,y,z) text(y,x,num2str(pmi_sent(x,y),'%1.1f'),'FontSize',12),row_p,col_p);
    ax.YAxis.Direction='reverse';
    ax.XAxis.Direction='normal';
    daspect([1,1,1])
    ax.XTick=[1:8];ax.XTickLabelRotation=90;
    ax.XTickLabel=(strsplit(all_sentence_pattern{1}));
    ax.YTick=[1:8];ax.YTickLabelRotation=0;
    ax.YTickLabel=(strsplit(all_sentence_pattern{1}));
    ax.Title.String='Pointwise Mutual Info (pmi)';
    %text(ax.XLim(1),ax.XLim(1)-1,'B','HorizontalAlignment','right','FontSize',20,'FontWeight','bold')
    % 
    %ax.YTick=[];
    ax.FontSize=14;
    ax.FontWeight='bold';
    
    ax=axes('position',[.61,.15,.22,.22*aspect_ration]);
    hold on 
    ax.XLim=[0,9];
    ax.YLim=[0,9];
   
    arrayfun(@(x,y,z) text(y,x,num2str(a_all_diffs_upper(x,y),'%1.1f'),'FontSize',12,'horizontalalignment','center'),row_p,col_p);
    ax.YAxis.Direction='reverse';
    ax.XAxis.Direction='normal';
   daspect([1,1,1])
    fliplr(strsplit(all_sentence_pattern{1}))
    ax.XTick=[1:8];ax.XTickLabelRotation=90;
    ax.XTickLabel=(strsplit(all_sentence_pattern{1}));
    ax.YTick=[1:8];ax.YTickLabelRotation=0;
    ax.YTickLabel=(strsplit(all_sentence_pattern{1}));
     ax.Title.String={'Change in',' gamma power'};
    % ax.YTick=[];
     ax.FontSize=14;
    ax.FontWeight='bold';
    
    %text(ax.XLim(1),ax.XLim(1)-1,'C','HorizontalAlignment','right','FontSize',20,'FontWeight','bold')
    ax=axes('position',[.85,.15,.1,.22*aspect_ration]);
    p=plot(a,[1:8],'Marker','.','color','k');
    p.MarkerSize=20;
    p.LineWidth=2
    ax.YLim=[0,9];
    ax.XLim=1.2*[min(a),max(a)];
    ax.XAxis.Direction='normal';
    ax.YAxis.Direction='reverse';
    ax.YAxis.Visible='off'
    ax.Box='off'
    ax.Title.String={'average','gamma','power'};
ax.FontSize=14;
    ax.FontWeight='bold';
   
    
    ax=axes('position',[.61,.82,.22,.05*aspect_ration]);
    b=barh(1,[.6],'DisplayName','wp','facecolor',[1,.3,.3])
    text(.65,1.2,'wp','FontSize',14,'FontWeight','bold')
    hold on
    b=barh(0,[.8],'DisplayName','pmi','facecolor',[.3,.3,1])
    text(.85,.2,'pmi','FontSize',14,'FontWeight','bold')
    %legend('position',[.81,.93,.05,.05])
     ax.XLim=[-.1,1.1]
    ax.YTick=[]
    ax.Box='off'
    ax.XAxis.Visible='off'
    ax.YAxis.Visible='off'
    ax.Title.String='R^2'
   ax.FontSize=14;

        if ~exist(strcat(analysis_path,info.subject))
            mkdir(strcat(analysis_path,info.subject))
        end 
        %set(gcf,'PaperPosition',[.25 .25 8 6])
        print(gcf, '-djpeg', strcat(analysis_path,info.subject,'/',info.subject,'_representitive_fig_presentation')); 
    
end
% 



%% 
%%  create a schematic of analysis pipline for presentation   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all 
look_out_window=100; %200 ms
offset=30; % 100 ms
word_range=[1:8];
num_rows=5;
num_columns=2;
total_plots=num_rows*num_columns;
p=0;
idx=3;
electrode_id=language_electrode_num(idx);
sentence_id=1;
for i=idx%:%length(language_electrode_num)
    elec_response_to_next_word_pattern=[];
    electrode_response=double(squeeze(session_sentence_hilbert_band_envelope_lang_elec_tensor(i,:,:)));
    electrode_response_sentence_cell=mat2cell(electrode_response,ones(1,size(electrode_response,1)),size(electrode_response,2)/8*ones(1,8));
    electrode_1st_sent_cell=[electrode_response_sentence_cell(sentence_id,:)];
    a=cell2mat(cellfun(@(x)  nanmean(x),electrode_1st_sent_cell,'UniformOutput',false));
   
    close all 
f=figure;
aspect_ration=9.32./4.13;
y=300;
set(f,'position',[591 455 aspect_ration*y y]);
ax=axes('position',[.1,.15,.38,.2*aspect_ration])
y1=electrode_response(sentence_id,:);
%y1=double(mean(elec_sentence_resp,2));
%e1=nanstd(elec_sentence_resp,0,2)./sqrt(sum(~isnan(elec_sentence_resp),2));
x=[1:length(y1)];
%[l,p] = boundedline(x, y1, e1, 'cmap',[1,0,0]);
l=plot(x,y1,'color',[.8,0,0])

l.DisplayName='S';
l.LineWidth=2;


ax.XLim=[-120,1080];
ax.XTick=[-120,1:135:1080];
ax.FontSize=14;
hold on 
ax.YLim=1.2*[min(y1),max(y1)];
ref=arrayfun(@(x)plot([x,x],ax.YLim,'k--'),ax.XTick(2:end));
Annotation=arrayfun(@(x) get(x,'Annotation'),ref);
hLegendEntry = arrayfun(@(x) get(x,'LegendInformation'),Annotation);
arrayfun(@(x) set(x,'IconDisplayStyle','off'),hLegendEntry);

ax.XTickLabel={'Word:','1','2','3','4','5','6','7','8'};

ax.YAxis.Label.String='a.u.'

ax.FontWeight='bold'
ax.Box='off'

    ax=axes('position',[.58,.15,.38,.2*aspect_ration])
    
    %text(ax.XLim(1),ax.XLim(1)-1,'C','HorizontalAlignment','right','FontSize',20,'FontWeight','bold')
    p=plot([1:8],a,'Marker','.','color',[.8,0,0]);
    p.MarkerSize=20;
    p.LineWidth=2
    ax.XLim=[0,9];
    ax.YLim=1.5*ax.YLim;
    ax.XAxis.Direction='normal';
    ax.XTick=[0:8];
    ax.XTickLabel={'Word:','1','2','3','4','5','6','7','8'};
ax.YAxis.Label.String='a.u.'
ax.FontSize=14;
ax.FontWeight='bold'
ax.Box='off'

    
   
    
    
        if ~exist(strcat(analysis_path,info.subject))
            mkdir(strcat(analysis_path,info.subject))
        end 
        %set(gcf,'PaperPosition',[.25 .25 8 6])
        print(gcf, '-djpeg', strcat(analysis_path,info.subject,'/',info.subject,'_representitive_fig_pmi_presentation')); 
    
end
% 



