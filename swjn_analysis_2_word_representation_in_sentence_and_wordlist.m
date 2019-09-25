%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 0: prepare the workspace 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all 
close all 
data_path='/Users/eghbalhosseiniasl1/MyData/ecog-sentence';
subject_id='AMC037';
d= dir([data_path,strcat('/**/',subject_id,'*_crunched.mat')]);
fprintf(' %d .mat files were found \n', length(d));
gamma_band_index=4;
save_path='/Users/eghbalhosseiniasl1/MyData/ECoG-sentence/crunched/';
analysis_path='/Users/eghbalhosseiniasl1/MyData/ECoG-sentence/analysis/swjn_analysis_2_word_representation_in_sentence_and_wordlist/';
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
    language_electrode=info.language_responsive_electrodes_hilbert_odd;
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
    hilbert_band_envelope_lang_elec_tensor=hilbert_band_envelope_tensor(find(language_electrode),:,:);

    session_sentence_hilbert_band_envelope_lang_elec_tensor=cat(2,session_sentence_hilbert_band_envelope_lang_elec_tensor,hilbert_band_envelope_lang_elec_tensor);
    % find sentence locations 
    example_sentence=cellfun(@(x) x(2:end),{sentences(:).trial_string},'UniformOutput',false);
    session_sentence_examples=[session_sentence_examples,example_sentence];
    example_sentence_locations=cellfun(@(x) (regexpi([node_cell{:,1}],x)),example_sentence,'UniformOutput',false);
    example_sentence_locations=cell2mat(cellfun(@(x) find_index(x), example_sentence_locations, 'UniformOutput',false));
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
    hilbert_band_envelope_lang_elec_tensor=hilbert_band_envelope_tensor(find(language_electrode),:,:);
    session_words_hilbert_band_envelope_lang_elec_tensor=cat(2,session_words_hilbert_band_envelope_lang_elec_tensor,hilbert_band_envelope_lang_elec_tensor);
    % 
    example_words=cellfun(@(x) x(2:end),{words(:).trial_string},'UniformOutput',false);
    session_wordlist_examples=[session_wordlist_examples,example_words];
    % 
    fprintf('added %s from %s \n',d(i).name, strcat(d(i).folder,'/',d(i).name));
    clear data subj
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  get common sentence and wordlists 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %
session_sentence_examples_split=transpose(cellfun(@(x) strsplit(x,' '),session_sentence_examples,'UniformOutput',false));
session_wordlist_examples_split=transpose(cellfun(@(x) strsplit(x,' '),session_wordlist_examples,'UniformOutput',false));
unique_session_sentence_examples_split=unique([session_sentence_examples_split{:}]);

% unique version 

unique_sentence_word_locations={};
unique_sentence_word_in_wordlist_locations={};
for k=1:length(unique_session_sentence_examples_split)
    sentence_patten_location=cell2mat(cellfun(@cell2mat,cellfun(@(z) cellfun(@(x) strcmp(x,unique_session_sentence_examples_split(k)),z,'UniformOutput',false),...
        session_sentence_examples_split,'UniformOutput',false),'UniformOutput',false));
    wordlist_patten_location=cell2mat(cellfun(@cell2mat,cellfun(@(z) cellfun(@(x) strcmp(x,unique_session_sentence_examples_split(k)),z,'UniformOutput',false),...
        session_wordlist_examples_split,'UniformOutput',false),'UniformOutput',false));
    unique_sentence_word_locations{k,1}=unique_session_sentence_examples_split{k};
    unique_sentence_word_locations{k,2}=sentence_patten_location;
    unique_sentence_word_locations{k,3}=sum(sentence_patten_location(:));
    %
    unique_sentence_word_in_wordlist_locations{k,1}=unique_session_sentence_examples_split{k};
    unique_sentence_word_in_wordlist_locations{k,2}=wordlist_patten_location;
    unique_sentence_word_in_wordlist_locations{k,3}=sum(wordlist_patten_location(:));
    
end


%% plot average for pre-post in word closing 
close all 
look_out_window=100; %200 ms
offset=30; % 100 ms
word_range=[1:8];
num_rows=1;
num_columns=8;
total_plots=num_rows*num_columns;
p=0;

for i=1:length(language_electrode_num)

    electrode_sentence_response=double(squeeze(session_sentence_hilbert_band_envelope_lang_elec_tensor(i,:,:)));
    electrode_response_sentence_cell=mat2cell(electrode_sentence_response,ones(1,size(electrode_sentence_response,1)),size(electrode_sentence_response,2)/8*ones(1,8));
    % 
    electrode_words_response=double(squeeze(session_words_hilbert_band_envelope_lang_elec_tensor(i,:,:)));
    electrode_response_words_cell=mat2cell(electrode_words_response,ones(1,size(electrode_words_response,1)),size(electrode_words_response,2)/8*ones(1,8));
    
    electrode_unique_sentence_pattern_response={};
    electrode_unique_wordlist_pattern_response={};
    for k=1:size(unique_sentence_word_locations,1)
        electrode_unique_sentence_pattern_response{k,2}=cell2mat(electrode_response_sentence_cell(unique_sentence_word_locations{k,2}));
        electrode_unique_sentence_pattern_response{k,1}=unique_sentence_word_locations{k,1};
        electrode_unique_wordlist_pattern_response{k,2}=cell2mat(electrode_response_words_cell(unique_sentence_word_in_wordlist_locations{k,2}));
        electrode_unique_wordlist_pattern_response{k,1}=unique_sentence_word_in_wordlist_locations{k,1};
    end
    shared_words=find([unique_sentence_word_in_wordlist_locations{:,3}]);
    electrode_unique_sentence_pattern_means=cellfun(@(x) nanmean(x(:,offset+[1:look_out_window]),2),electrode_unique_sentence_pattern_response(shared_words,2),'UniformOutput',false);
    electrode_unique_wordlist_pattern_means=cellfun(@(x) nanmean(x(:,offset+[1:look_out_window]),2),electrode_unique_wordlist_pattern_response(shared_words,2),'UniformOutput',false);
    colors=flipud(cbrewer('seq','YlOrRd',length(shared_words)+10));
    colors_words=flipud(cbrewer('seq','Blues',length(shared_words)+10));
    % 
    figure(fix((i-1)/total_plots)+1);
    set(gcf,'position',[1822,-156,1257,1269])
    sub_title=sprintf('subj: %s,\n electrodes: %d ', info.subject ,language_electrode_num(i));
    ax=subplot(num_rows,num_columns,i-total_plots*fix((i-1)/total_plots));
    
    bl=barh([1:length(shared_words)]+.2,cellfun(@nanmean,electrode_unique_sentence_pattern_means),'Facecolor','flat','BarWidth',.4);
    set(bl,'Linestyle','none');
    set(bl,'Displayname','sentences')
    bl.CData=colors(1:length(shared_words),:);
    hold on 
    e=arrayfun(@(x,y,z) errorbar(x,y,z,'horizontal','k.'),[cellfun(@nanmean,electrode_unique_sentence_pattern_means)]',bl.XData,[cellfun(@nanstd,electrode_unique_sentence_pattern_means)]');
    arrayfun(@(x) set(e,'Capsize',0),e);arrayfun(@(x) set(e,'Linewidth',1.5),e);arrayfun(@(x) set(e,'marker','none'),e);
    arrayfun(@(x,y) set(x,'Color',bl.CData(y,:)),e,[1:size(bl.CData,1)]);
     hAnnotation=arrayfun(@(x) get(x,'Annotation'),e);
    hLegendEntry = arrayfun(@(x) get(x,'LegendInformation'),hAnnotation);
    arrayfun(@(x) set(x,'IconDisplayStyle','off'),hLegendEntry);
    % 
    bl=barh([1:length(shared_words)]-.2,cellfun(@nanmean,electrode_unique_wordlist_pattern_means),'FaceColor','flat','BarWidth',.4);
    set(bl,'Linestyle','none');
    bl.CData=colors_words(1:length(shared_words),:);
    set(bl,'Displayname','words');
    e=arrayfun(@(x,y,z) errorbar(x,y,z,'horizontal','k.'),[cellfun(@nanmean,electrode_unique_wordlist_pattern_means)]',bl.XData,[cellfun(@nanstd,electrode_unique_wordlist_pattern_means)]');
    arrayfun(@(x) set(e,'Capsize',0),e);arrayfun(@(x) set(e,'Linewidth',1.5),e);arrayfun(@(x) set(e,'marker','none'),e);
   
    arrayfun(@(x,y) set(x,'Color',bl.CData(y,:)),e,[1:size(bl.CData,1)])
    hAnnotation=arrayfun(@(x) get(x,'Annotation'),e);
    hLegendEntry = arrayfun(@(x) get(x,'LegendInformation'),hAnnotation);
    arrayfun(@(x) set(x,'IconDisplayStyle','off'),hLegendEntry);
    hold on 

    ax.XAxis.Visible = 'off'; 
    set(ax,'ydir', 'normal','box','off');
    set(ax,'ytick',[1:length(shared_words)],'YTickLabel',[unique_sentence_word_in_wordlist_locations(shared_words,1)]);
   
    set(gca,'TickLength',[0.003, 0.003])
    title(sub_title);
    if ~mod(i,total_plots) | i==length(language_electrode_num)
        legend('show','Location','best')
        p=p+1;
        
        %print(gcf, '-depsc', strcat(analysis_path,info.subject,'/',info.subject,'_node_closing',num2str(p),'_window_',window,'.eps')); 
    end 
    
end

%%
emb = fastTextWordEmbedding;


%% 
close all 
look_out_window=100; %200 ms
offset=30; % 100 ms
word_range=[1:8];
num_rows=1;
num_columns=8;
total_plots=num_rows*num_columns;
p=0;

for i=1:length(language_electrode_num)

    electrode_sentence_response=double(squeeze(session_sentence_hilbert_band_envelope_lang_elec_tensor(i,:,:)));
    electrode_response_sentence_cell=mat2cell(electrode_sentence_response,ones(1,size(electrode_sentence_response,1)),size(electrode_sentence_response,2)/8*ones(1,8));
    % 
    electrode_words_response=double(squeeze(session_words_hilbert_band_envelope_lang_elec_tensor(i,:,:)));
    electrode_response_words_cell=mat2cell(electrode_words_response,ones(1,size(electrode_words_response,1)),size(electrode_words_response,2)/8*ones(1,8));
    
    electrode_unique_sentence_pattern_response={};
    electrode_unique_wordlist_pattern_response={};
    for k=1:size(unique_sentence_word_locations,1)
        electrode_unique_sentence_pattern_response{k,2}=cell2mat(electrode_response_sentence_cell(unique_sentence_word_locations{k,2}));
        electrode_unique_sentence_pattern_response{k,3}=word2vec(emb,unique_sentence_word_locations{k,1});
        electrode_unique_sentence_pattern_response{k,1}=unique_sentence_word_locations{k,1};
        electrode_unique_wordlist_pattern_response{k,2}=cell2mat(electrode_response_words_cell(unique_sentence_word_in_wordlist_locations{k,2}));
        electrode_unique_wordlist_pattern_response{k,1}=unique_sentence_word_in_wordlist_locations{k,1};
        electrode_unique_wordlist_pattern_response{k,3}=word2vec(emb,unique_sentence_word_in_wordlist_locations{k,1});
    end
    
    shared_words=find([unique_sentence_word_in_wordlist_locations{:,3}]);
    temp_word2vec=cell2mat(cellfun(@(x) x,electrode_unique_wordlist_pattern_response(shared_words,3),'UniformOutput',false));
    shared_words(isnan(sum(temp_word2vec,2)))=[];
    electrode_unique_sentence_pattern_means=cellfun(@(x) nanmean(x(:,offset+[1:look_out_window]),2),electrode_unique_sentence_pattern_response(shared_words,2),'UniformOutput',false);
    electrode_unique_wordlist_pattern_means=cellfun(@(x) nanmean(x(:,offset+[1:look_out_window]),2),electrode_unique_wordlist_pattern_response(shared_words,2),'UniformOutput',false);  
    % 
    electrode_unique_sentence_pattern_word2vec=cell2mat(cellfun(@(x) x,electrode_unique_wordlist_pattern_response(shared_words,3),'UniformOutput',false));
    electrode_sentence_pattern_word2vec=cell2mat(cellfun(@(x) x,electrode_unique_wordlist_pattern_response(:,3),'UniformOutput',false));

    % compute the cosine similarity in the word embeddings 
    % 
    Y = pdist(electrode_unique_sentence_pattern_word2vec,'cosine');
    Z = linkage(Y,'average'); 
    leafOrder=optimalleaforder(Z,Y);
    Y = pdist(electrode_unique_sentence_pattern_word2vec(leafOrder,:),'cosine');
    Z = linkage(Y,'average');
    c = cophenet(Z,Y);
    % 
    figure;
    set(gcf,'position',[1442 1 1275 1021])
    colors = cbrewer('div', 'RdBu', 128);
    colors=flipud(colors);
    colormap(colors);
    D=squareform(Y);
    subplot('position',[.07,.07,.7,.7])
    imagesc(D);
    tick_labels=[unique_sentence_word_in_wordlist_locations(shared_words,1)];
    set(gca,'box','off','ytick',[1:size(D,1)],'yticklabel',tick_labels(leafOrder),...
        'xtick',[1:size(D,1)],'xticklabel',tick_labels(leafOrder));
    xtickangle(90);
    main_fig=findobj(gcf,'type','figure');
    aspect_ratio=main_fig.Position(3)./main_fig.Position(4);
    main_fig = findobj(gca,'Type','axes');
    axpos = get(main_fig, 'Position');
    axpos(4)= axpos(4);
    axpos(3)=min([axpos(3),axpos(4)]);
    axpos(4)=min([axpos(3),axpos(4)]);
    axpos(3)=axpos(3)./aspect_ratio;
    set(main_fig, 'Position', axpos);
    % 
    colors=flipud(cbrewer('seq','YlOrRd',length(shared_words)+20));
    colors_words=flipud(cbrewer('seq','Blues',length(shared_words)+20));
    % get axis limits 
    all_data_points=[cell2mat(electrode_unique_sentence_pattern_means);cell2mat(electrode_unique_wordlist_pattern_means)];
    data_quantiles=quantile(all_data_points,50);
    % plot sentences 
    sub_title=sprintf('subj: %s,\n electrodes: %d ', info.subject ,language_electrode_num(i));
    ax=subplot('position',[main_fig.Position(1),main_fig.Position(2)+main_fig.Position(4),main_fig.Position(3),.10]);
    % 
    bl=bar([1:length(shared_words)],cellfun(@nanmean,electrode_unique_sentence_pattern_means(leafOrder)),'Facecolor','flat');
    xlim([1-bl.BarWidth/2,70+bl.BarWidth/2])
    set(bl,'Linestyle','none');
    set(bl,'Displayname','sentences')
    bl.CData=colors(1:length(shared_words),:);
    hold on 
    e=arrayfun(@(x,y,z) errorbar(x,y,z,'vertical','k.'),bl.XData,[cellfun(@nanmean,electrode_unique_sentence_pattern_means(leafOrder))]',[cellfun(@nanstd,electrode_unique_sentence_pattern_means(leafOrder))]');
    arrayfun(@(x) set(e,'Capsize',0),e);arrayfun(@(x) set(e,'Linewidth',2),e);arrayfun(@(x) set(e,'marker','none'),e);
    arrayfun(@(x,y) set(x,'Color',bl.CData(y,:)),e,[1:size(bl.CData,1)]);
     hAnnotation=arrayfun(@(x) get(x,'Annotation'),e);
    hLegendEntry = arrayfun(@(x) get(x,'LegendInformation'),hAnnotation);
    arrayfun(@(x) set(x,'IconDisplayStyle','off'),hLegendEntry);
    ax.XAxis.Visible = 'off'; 
    ax.Box='off';
    ax.YLim=[0+0*data_quantiles(1),data_quantiles(end)];
    legend('show','Location','best')
    %title(sub_title);
    ax=subplot('position',[main_fig.Position(1),main_fig.Position(2)+.10+main_fig.Position(4),main_fig.Position(3),.10]);
    plot(normalize(sum(D)),'-','linewidth',2,'color',[0.5,.5,.5]);
    hold on 
    plot(normalize([cellfun(@nanmean,electrode_unique_sentence_pattern_means(leafOrder))]),'r','LineWidth',3);
    ax.XAxis.Visible='off';
    ax.YAxis.Visible='off';
    text(ax.XLim(2),ax.YLim(1),strcat('correlation: ',num2str(corr(normalize([cellfun(@nanmean,electrode_unique_sentence_pattern_means(leafOrder))]),...
        normalize(sum(D))'))),'horizontalalignment','left','FontSize',16)
    % 
    set(ax,'ylim',[-4,4])
    ax=subplot('position',[main_fig.Position(1)+main_fig.Position(3),main_fig.Position(2),.10,main_fig.Position(4)]);
    bl=barh([1:length(shared_words)],cellfun(@nanmean,electrode_unique_wordlist_pattern_means(leafOrder)),'FaceColor','flat');
    ylim([1-bl.BarWidth/2,70+bl.BarWidth/2]);
    set(bl,'Linestyle','none');
    bl.CData=colors_words(1:length(shared_words),:);
    set(bl,'Displayname','words');
    hold on 
    e=arrayfun(@(x,y,z) errorbar(x,y,z,'horizontal','k.'),[cellfun(@nanmean,electrode_unique_wordlist_pattern_means(leafOrder))]',bl.XData,[cellfun(@nanstd,electrode_unique_wordlist_pattern_means(leafOrder))]');
    arrayfun(@(x) set(e,'Capsize',0),e);arrayfun(@(x) set(e,'Linewidth',2),e);arrayfun(@(x) set(e,'marker','none'),e);
    arrayfun(@(x,y) set(x,'Color',bl.CData(y,:)),e,[1:size(bl.CData,1)])
    hAnnotation=arrayfun(@(x) get(x,'Annotation'),e);
    hLegendEntry = arrayfun(@(x) get(x,'LegendInformation'),hAnnotation);
    arrayfun(@(x) set(x,'IconDisplayStyle','off'),hLegendEntry);
    ax.XAxis.Visible = 'on'; 
    ax.YAxis.Visible = 'off';
    set(ax,'ydir', 'reverse','box','off');
    ax.XLim=[0+0*data_quantiles(1),data_quantiles(end)];
    legend('show','Location','best');
    ax=subplot('position',[main_fig.Position(1)+main_fig.Position(3)+.1,main_fig.Position(2),.10,main_fig.Position(4)]);
    
    plot(normalize(sum(D)),1:length(sum(D)),'-','linewidth',2,'color',[0.5,.5,.5]);
    hold on 
    plot(normalize([cellfun(@nanmean,electrode_unique_wordlist_pattern_means(leafOrder))]),1:length(sum(D)),'color',[.1,.1,.8],'LineWidth',3);
    ax.XAxis.Visible='off';
    ax.YAxis.Visible='off';
    set(ax,'ydir', 'reverse','box','off');
    set(ax,'xlim',[-4,4])
    text(ax.XLim(2),ax.YLim(1),strcat('correlation: ',num2str(corr(normalize([cellfun(@nanmean,electrode_unique_wordlist_pattern_means(leafOrder))]),...
        normalize(sum(D))'))),'horizontalalignment','left','FontSize',16)
    
    
    
    if ~exist(strcat(analysis_path,info.subject))
            mkdir(strcat(analysis_path,info.subject))
    end   
    print(gcf, '-depsc', strcat(analysis_path,info.subject,'/',info.subject,'_ch_',num2str(language_electrode_num(i)),'_representation_sentence_vs_words.eps')); 
    close(gcf)
    
end

%% plot pre post for word closings 
close all 
look_out_window=word_length; %200 ms
offset=0; % 100 ms
word_range=[1:8];
num_rows=5;
num_columns=4;
total_plots=num_rows*num_columns;
for i=1:length(language_electrode_num)
    elec_response_to_node_closing_pattern=[];
    electrode_sentence_response=squeeze(session_sentence_hilbert_band_envelope_lang_elec_tensor(i,:,:));
    electrode_response_sentence_cell=mat2cell(electrode_sentence_response,ones(1,size(electrode_sentence_response,1)),size(electrode_sentence_response,2)/8*ones(1,8));
    node_closing={};
    node_closing_index_sentence={};
    for k=1:8
        closing_type=total_closing(k,:);
        node_closing{k,1}=[ electrode_response_sentence_cell(find(cell2mat([closing_type{2}]))),electrode_response_sentence_cell(find(cell2mat([closing_type{3}])))];
        node_closing_index_sentence{k,1}=k*ones(size([electrode_response_sentence_cell(find(cell2mat([closing_type{2}])))]))-1;
    end
    a=cellfun(@(x) transpose(cell2mat(x)),node_closing,'UniformOutput',false);
    b=cellfun(@(x) transpose(x),node_closing_index_sentence,'UniformOutput',false);
    a=[a{:}];
    b=[b{:}];
    a=a-repmat(nanmean(a,1),size(a,1),1);
    colors=flipud(cbrewer('div','RdBu',10));
    
    % 
    x=1000/info.downsample_sampling_rate*   repmat([-(word_length-1):word_length]',1,size(a,2));
    figure(fix((i-1)/total_plots)+1);
   set(gcf,'position',[1822,-156,1257,1269])
    sub_title=sprintf('subj: %s,\n electrodes: %d ', info.subject ,language_electrode_num(i));
    ax=subplot(num_rows,num_columns,i-total_plots*fix((i-1)/total_plots));
    colors=flipud(cbrewer('seq','OrRd',length(unique(b))+3));
    for k=unique(b)
        ind=find(b==k);
        bl = boundedline(double(nanmean(x(:,ind),2)'), double(nanmean(a(:,ind),2)'), double(nanstd(a(:,ind),0,2)')./sqrt(length(ind)), 'cmap', colors(k+1,:),'alpha',ax,'transparency',.0);
        bl.LineWidth=2;
        bl.DisplayName=num2str(k);    
    end 
    hold on 
    h1=plot(double(nanmean(x,2)),0*double(nanmean(x,2)),'k--');hAnnotation = get(h1,'Annotation');hLegendEntry = get(hAnnotation,'LegendInformation');set(hLegendEntry,'IconDisplayStyle','off');
    h1=plot([0,0],get(gca,'ylim'),'k--');hAnnotation = get(h1,'Annotation');hLegendEntry = get(hAnnotation,'LegendInformation');set(hLegendEntry,'IconDisplayStyle','off');
    xlim([min(x(:)),max(x(:))])
    set(ax,'ydir', 'normal','box','off');
    title(sub_title);
    if ~mod(i,total_plots) | i==length(language_electrode_num)
    legend('show','Location','northeastoutside')
    end 

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  do a regression for pre --> post signal change and node closing number 
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
    elec_response_to_node_closing_pattern=[];
    electrode_sentence_response=double(squeeze(session_sentence_hilbert_band_envelope_lang_elec_tensor(i,:,:)));
    electrode_response_sentence_cell=mat2cell(electrode_sentence_response,ones(1,size(electrode_sentence_response,1)),size(electrode_sentence_response,2)/8*ones(1,8));
    node_closing={};
    node_closing_index_sentence={};
    for k=1:8
        closing_type=total_closing(k,:);
        node_closing{k,1}=[ electrode_response_sentence_cell(find(cell2mat([closing_type{2}]))),electrode_response_sentence_cell(find(cell2mat([closing_type{3}])))];
        node_closing_index_sentence{k,1}=k*ones(size([electrode_response_sentence_cell(find(cell2mat([closing_type{2}])))]))-1;
    end
    
    closing_means_sentence=cellfun(@(x) mat2cell(transpose(repmat(nanmean(transpose(cell2mat(x)),1),size(x,2),1)),ones(1,size(x,1)),ones(1,size(x,2))),node_closing,'UniformOutput',false);
    normalize_node_closing_sentence=cellfun(@(x,y) cellfun(@(p,q) p-q,x,y,'uniformoutput',false),node_closing,closing_means_sentence,'UniformOutput',false);
    a=cellfun(@(x) cellfun(@(y) nanmean(y(offset+[1:look_out_window])),x),normalize_node_closing_sentence,'UniformOutput',false);
    normalize_node_closing_diff=cellfun(@(x) transpose(diff(x,1,2)),a,'UniformOutput',false);
    
    b=cellfun(@(x) transpose(x),node_closing_index_sentence,'UniformOutput',false);
    a=[normalize_node_closing_diff{:}]';
    b=[b{:}]';
    % remove 5
    remove_5_closing=find(b==5);
    b(remove_5_closing)=[];
    a(remove_5_closing)=[];
    X=[ones(size(b)),b];
    lm = fitlm(b,a);
    a_fit=lm.Fitted;
    colors=flipud(cbrewer('div','RdBu',10));
    % 
    figure(fix((i-1)/total_plots)+1);
    set(gcf,'position',[2533,21,705,950])
    sub_title=sprintf('subj: %s,\n electrodes: %d ', info.subject ,language_electrode_num(i));
    ax=subplot(num_rows,num_columns,i-total_plots*fix((i-1)/total_plots));
    colors=flipud(cbrewer('seq','OrRd',length(unique(b))+3));
    c=colors(b+1,:);
    b_jitter=beeswarm(b,a,'hex','random',.2);
    bl=scatter(b_jitter,a,20,c,'filled');
    
    hold on
    plot(b,a_fit,'k-','LineWidth',2)
    %xlim([min(x(:))-150,max(x(:))+250])
    h1=plot(get(gca,'xlim'),0*get(gca,'ylim'),'k--');hAnnotation = get(h1,'Annotation');hLegendEntry = get(hAnnotation,'LegendInformation');set(hLegendEntry,'IconDisplayStyle','off');
    ax.XAxis.Visible = 'on'; 
    y_quantile=quantile(a,50);
    set(ax,'ydir', 'normal','box','off','ylim',[y_quantile(1),y_quantile(end)]);
    text(max(get(gca,'xlim')),min(get(gca,'ylim')),{strcat('slope: ',num2str(lm.Coefficients.Estimate(2))),strcat('p_{value}: ',num2str(lm.Coefficients.pValue(2)))},...
        'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10,'fontweight','bold')
    title(sub_title);
    
    if ~mod(i,total_plots) | i==length(language_electrode_num)
        legend('show','Location','northeastoutside')
        p=p+1;
        if ~exist(strcat(analysis_path,info.subject))
            mkdir(strcat(analysis_path,info.subject))
        end 
        %window=sprintf('%d_%dms',floor(min(1000/info.downsample_sampling_rate*(offset+[1:look_out_window]))),...
        %    floor(max(1000/info.downsample_sampling_rate*(offset+[1:look_out_window]))));
        %print(gcf, '-depsc', strcat(analysis_path,info.subject,'/',info.subject,'_node_closing',num2str(p),'_window_',window,'.eps')); 
    end 
    
end

%% plot average for pre-post in word closing 
close all 
look_out_window=100; %200 ms
offset=30; % 100 ms
word_range=[1:8];
num_rows=5;
num_columns=2;
total_plots=num_rows*num_columns;
p=0;
for i=1:length(language_electrode_num)
    elec_response_to_node_closing_pattern=[];
    electrode_sentence_response=squeeze(session_sentence_hilbert_band_envelope_lang_elec_tensor(i,:,:));
    electrode_response_sentence_cell=mat2cell(electrode_sentence_response,ones(1,size(electrode_sentence_response,1)),size(electrode_sentence_response,2)/8*ones(1,8));
    node_closing={};
    node_closing_index_sentence={};
    for k=1:8
        closing_type=total_closing(k,:);
        node_closing{k,1}=[ electrode_response_sentence_cell(find(cell2mat([closing_type{2}]))),electrode_response_sentence_cell(find(cell2mat([closing_type{3}])))];
        node_closing_index_sentence{k,1}=k*ones(size([electrode_response_sentence_cell(find(cell2mat([closing_type{2}])))]))-1;
    end
    closing_means_sentence=cellfun(@(x) mat2cell(transpose(repmat(nanmean(transpose(cell2mat(x)),1),size(x,2),1)),ones(1,size(x,1)),ones(1,size(x,2))),node_closing,'UniformOutput',false);
    normalize_node_closing_sentence=cellfun(@(x,y) cellfun(@(p,q) p-q,x,y,'uniformoutput',false),node_closing,closing_means_sentence,'UniformOutput',false);
    a=cellfun(@(x) cellfun(@(y) nanmean(y(offset+[1:look_out_window])),x),normalize_node_closing_sentence,'UniformOutput',false);
    std_val=cellfun(@(x) cellfun(@(y) nanstd(y(offset+[1:look_out_window])),x),normalize_node_closing_sentence,'UniformOutput',false);
    %a=cellfun(@(x) cellfun(@nanmean ,x),normalize_node_closing,'UniformOutput',false);
    a=cellfun(@(x) transpose(x),a,'UniformOutput',false);
    std_val=cellfun(@(x) transpose(x),std_val,'UniformOutput',false);
    b=cellfun(@(x) transpose(x),node_closing_index_sentence,'UniformOutput',false);
    a=[a{:}];
    b=[b{:}];
    std_val=[std_val{:}];
    colors=flipud(cbrewer('div','RdBu',10));
    % 
    x=1000/info.downsample_sampling_rate*   repmat([-134/2,word_length/2]',1,size(a,2));
    figure(fix((i-1)/total_plots)+1);
   set(gcf,'position',[1822,-156,1257,1269])
    sub_title=sprintf('subj: %s,\n electrodes: %d ', info.subject ,language_electrode_num(i));
    ax=subplot(num_rows,num_columns,i-total_plots*fix((i-1)/total_plots));
    colors=flipud(cbrewer('seq','OrRd',length(unique(b))+3));
    means_sentence=[];
    stds_sentence=[];
    for k=unique(b)
        ind=find(b==k);
        means_sentence =[means_sentence;double(nanmean(a(:,ind),2)')];
        stds_sentence=[stds_sentence;double(nanstd(a(:,ind),0,2)')./sqrt(length(ind))];
    end 
    bl=bar(means_sentence');
    arrayfun(@(x,y) set(x,'FaceColor',colors(y,:)),bl,[1:size(bl,2)])
    arrayfun(@(x,y) set(x,'Displayname',num2str(y)),bl,[1:size(bl,2)]-1)
    arrayfun(@(x,y) set(x,'EdgeColor',[1,1,1]),bl);
    hold on
    e=errorbar(repmat(bl(1).XData,size(bl,2),1)+transpose(arrayfun(@(x) get(x,'XOffset'),bl)),means_sentence,stds_sentence,'k.');
    arrayfun(@(x) set(e,'Capsize',0),e);arrayfun(@(x) set(e,'Linewidth',2),e);
    hAnnotation=arrayfun(@(x) get(x,'Annotation'),e);
    hLegendEntry = arrayfun(@(x) get(x,'LegendInformation'),hAnnotation);
    arrayfun(@(x) set(x,'IconDisplayStyle','off'),hLegendEntry);
    hold on 
    %xlim([min(x(:))-150,max(x(:))+250])
    h1=plot(get(gca,'xlim'),0*get(gca,'ylim'),'k-');hAnnotation = get(h1,'Annotation');hLegendEntry = get(hAnnotation,'LegendInformation');set(hLegendEntry,'IconDisplayStyle','off');
    ax.XAxis.Visible = 'off'; 
    set(ax,'ydir', 'normal','box','off');
    set(ax,'xticklabel',{'pre','poset'});
    title(sub_title);
    if ~mod(i,total_plots) | i==length(language_electrode_num)
        legend('show','Location','northeastoutside')
        p=p+1;
        if ~exist(strcat(analysis_path,info.subject))
            mkdir(strcat(analysis_path,info.subject))
        end 
        window=sprintf('%d_%dms',floor(min(1000/info.downsample_sampling_rate*(offset+[1:look_out_window]))),...
            floor(max(1000/info.downsample_sampling_rate*(offset+[1:look_out_window]))));
        print(gcf, '-depsc', strcat(analysis_path,info.subject,'/',info.subject,'_node_closing',num2str(p),'_window_',window,'.eps')); 
    end 
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  compute correlation of actvitiy with the number of open nodes 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all 
anchor_points=1:word_length:size(session_sentence_hilbert_band_envelope_lang_elec_tensor,3);
look_out_window=100; %200 ms
offset=15; % 100 ms
p_val_threshold=0.05;
word_range=[4:8];
j=0;
[row,column]=ind2sub(size(ones(2,3)),find(ones(2,3)));
for i=1:size(session_sentence_hilbert_band_envelope_lang_elec_tensor,1)
    electrode_sentence_response=squeeze(session_sentence_hilbert_band_envelope_lang_elec_tensor(i,:,:));    
  
    electrode_response_sentence_cell=mat2cell(electrode_sentence_response,ones(1,size(electrode_sentence_response,1)),word_length*ones(1,8));
    electrode_post_word_sentence_mean=double(cell2mat(cellfun(@(x) nanmean(x(offset+(1:look_out_window))),electrode_response_sentence_cell,'UniformOutput',false)));
    word_position=repmat(1:size(electrode_post_word_sentence_mean,2),[size(electrode_post_word_sentence_mean,1),1]);
    % 
    electrode_response_sentence_cell=cellfun(@transpose,mat2cell(electrode_post_word_sentence_mean,ones(1,size(electrode_post_word_sentence_mean,1)),size(electrode_post_word_sentence_mean,2)),'UniformOutput',false);
    all_sentence_pattern_cell=cellfun(@transpose,mat2cell(all_sentence_pattern,ones(1,size(all_sentence_pattern,1)),size(all_sentence_pattern,2)),'uniformoutput',false);
    word_position_cell=cellfun(@transpose,mat2cell(word_position,ones(1,size(word_position,1)),size(word_position,2)),'uniformoutput',false);
    % select a subset of words 
    electrode_response_sentence_cell=cellfun(@(x) x(word_range),electrode_response_sentence_cell,'UniformOutput',false);
    all_sentence_pattern_cell=cellfun(@(x) x(word_range),all_sentence_pattern_cell,'UniformOutput',false);
    word_position_cell=cellfun(@(x) x(word_range),word_position_cell,'UniformOutput',false);
    
    [rho,p_val]=cellfun(@(x,y) corr(x,y,'type','Spearman'),electrode_response_sentence_cell,all_sentence_pattern_cell,'UniformOutput',false);
    elec_response_to_open_node_corr=[rho{:}];
    elec_response_to_open_node_corr_p_val=[p_val{:}];
    elec_response_to_open_node_corr_significant=elec_response_to_open_node_corr_p_val<p_val_threshold;
    
    % 
    [rho,p_val]=cellfun(@(x,y) corr(x,y,'type','Spearman'),electrode_response_sentence_cell,word_position_cell,'UniformOutput',false);
    
    elec_response_to_word_position_corr=[rho{:}];
    elec_response_to_word_position_corr_p_val=[p_val{:}];
    elec_response_to_word_position_corr_significant=elec_response_to_word_position_corr_p_val<p_val_threshold;
    % 
    significance_overlap= 1 |elec_response_to_open_node_corr_significant | elec_response_to_word_position_corr_significant;
    elec_response_to_open_node_corr=elec_response_to_open_node_corr(significance_overlap);
    elec_response_to_word_position_corr=elec_response_to_word_position_corr(significance_overlap);
    
    j=j+1;
    g(row(j),column(j))=gramm('x',elec_response_to_word_position_corr,'y',elec_response_to_open_node_corr);
    g(row(j),column(j)).geom_point('alpha','1');
    g(row(j),column(j)).stat_cornerhist('edges',-.75:0.025:.75,'aspect',.8,'location',1.1);
    g(row(j),column(j)).geom_abline();
    g(row(j),column(j)).axe_property('DataAspectRatio',[1 1 1]);
    g(row(j),column(j)).set_title(strcat(info.subject,' ch# ',num2str(language_electrode_num(i))));
    g(row(j),column(j)).set_names('x',' ','y',' ');
    g(1,1).set_names('x','correlation to word position','y','correlation to open nodes');
    if j==6 | i==length(language_electrode_num)
    g.set_title('correlation between gamma power and sentence');
    figure('Position',[100 100 800 550]);
    g.draw();
    j=0;
     % g.export('file_name',strcat(info.subject,'_channel_correlation_',num2str(ceil(i/6)),'_window_',num2str(ceil(offset*3.33))...
    %        ,'_',num2str(ceil((offset+look_out_window)*3.33))),'export_path',strcat(save_path),'file_type','pdf');
    end 
end 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  compute correlation based on unique patterns  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all 
unique_sentence_pattern=mat2cell(unique(all_sentence_pattern,'row'),ones(1,size(unique(all_sentence_pattern,'row'),1)),size(all_sentence_pattern,2));

uniq_sent_patt_loc_in_all_sentences=cellfun(@(x) find(ismember(all_sentence_pattern,x,'row')),unique_sentence_pattern,'UniformOutput',false);
look_out_window=100; %200 ms
offset=15; % 100 ms
word_range=[1:8];
Beta_mat_unique=[];
[row,column]=ind2sub(size(ones(2,3)),find(ones(2,3)));

j=0;
for i=1:length(language_electrode_num)
    elec_response_to_open_node_pattern=[];
    electrode_sentence_response=squeeze(session_sentence_hilbert_band_envelope_lang_elec_tensor(i,:,:));
    electrode_response_sentence_cell=mat2cell(electrode_sentence_response,ones(1,size(electrode_sentence_response,1)),word_length*ones(1,8));
    electrode_post_word_sentence_mean=cell2mat(cellfun(@(x) nanmean(x(offset+(1:look_out_window))),electrode_response_sentence_cell,'UniformOutput',false));
    elec_response_to_open_node_corr=[];
    elec_response_to_word_position_corr=[];
    for k=1:size(unique_sentence_pattern,1)
        pattern_response=mean(electrode_post_word_sentence_mean(uniq_sent_patt_loc_in_all_sentences{k},:),1);
        elec_response_to_open_node_pattern=[elec_response_to_open_node_pattern;pattern_response];
    end
    electrode_response_sentence_cell=cellfun(@transpose,mat2cell(elec_response_to_open_node_pattern,ones(1,size(elec_response_to_open_node_pattern,1)),size(elec_response_to_open_node_pattern,2)),'UniformOutput',false);
    unique_sentence_pattern_cell=cellfun(@transpose,unique_sentence_pattern,'uniformoutput',false);
    word_position=repmat(1:size(elec_response_to_open_node_pattern,2),[size(elec_response_to_open_node_pattern,1),1]);
    word_position_cell=cellfun(@transpose,mat2cell(word_position,ones(1,size(word_position,1)),size(word_position,2)),'uniformoutput',false);
    % select a subset of words
    electrode_response_sentence_cell=cellfun(@(x) x(word_range),electrode_response_sentence_cell,'UniformOutput',false);
    unique_sentence_pattern_cell=cellfun(@(x) x(word_range),unique_sentence_pattern_cell,'UniformOutput',false);
    word_position_cell=cellfun(@(x) x(word_range),word_position_cell,'UniformOutput',false);
    
    [rho,p_val]=cellfun(@(x,y) corr(x,y,'type','Spearman'),electrode_response_sentence_cell,unique_sentence_pattern_cell,'UniformOutput',false);
    elec_response_to_open_node_corr=[rho{:}];
    elec_response_to_open_node_corr_p_val=[p_val{:}];
    elec_response_to_open_node_corr_significant=elec_response_to_open_node_corr_p_val<p_val_threshold;
    
    %
    [rho,p_val]=cellfun(@(x,y) corr(x,y,'type','Spearman'),electrode_response_sentence_cell,word_position_cell,'UniformOutput',false);
    
    elec_response_to_word_position_corr=[rho{:}];
    elec_response_to_word_position_corr_p_val=[p_val{:}];
    elec_response_to_word_position_corr_significant=elec_response_to_word_position_corr_p_val<p_val_threshold;
    %
    significance_overlap= elec_response_to_open_node_corr_significant + 2*elec_response_to_word_position_corr_significant+1;
    elec_response_to_open_node_corr=elec_response_to_open_node_corr;
    elec_response_to_word_position_corr=elec_response_to_word_position_corr;
    conditions={'none','Node', 'Word', 'both'};
    c=conditions(significance_overlap);
    % do the do the correlation;
    g.set_color_options('map','d3_20b');
    j=j+1;
    g(row(j),column(j))=gramm('x',elec_response_to_word_position_corr,'y',elec_response_to_open_node_corr,'color',c);
    g(row(j),column(j)).geom_point('alpha','1');
    g(row(j),column(j)).stat_cornerhist('edges',-.75:0.025:.75,'aspect',.8,'location',1.1);
    g(row(j),column(j)).geom_abline();
    g(row(j),column(j)).axe_property('DataAspectRatio',[1 1 1]);
    g(row(j),column(j)).set_title(strcat(info.subject,' ch# ',num2str(language_electrode_num(i))));
    g(row(j),column(j)).set_names('x',' ','y',' ','color','pVal');
%     if j==1
%         g(row(j),column(j)).no_legend();
%     end 
    
    g(1,1).set_names('x','correlation to word position','y','correlation to open nodes');
    if j==6 | i==length(language_electrode_num)
        g.set_title('correlation between gamma power and unique open node vs. word position');
        figure('Position',[100 100 800 550]);
         g.set_color_options('map','lch');
        g.draw();
        j=0;
    end
    
   
end

%% plot pre post for word closings 
close all 
look_out_window=word_length; %200 ms
offset=0; % 100 ms
word_range=[1:8];
num_rows=5;
num_columns=4;
total_plots=num_rows*num_columns;
for i=1:length(language_electrode_num)
    elec_response_to_node_closing_pattern=[];
    electrode_sentence_response=squeeze(session_sentence_hilbert_band_envelope_lang_elec_tensor(i,:,:));
    electrode_response_sentence_cell=mat2cell(electrode_sentence_response,ones(1,size(electrode_sentence_response,1)),size(electrode_sentence_response,2)/8*ones(1,8));
    node_closing={};
    node_closing_index_sentence={};
    for k=1:8
        closing_type=total_closing(k,:);
        node_closing{k,1}=[ electrode_response_sentence_cell(find(cell2mat([closing_type{2}]))),electrode_response_sentence_cell(find(cell2mat([closing_type{3}])))];
        node_closing_index_sentence{k,1}=k*ones(size([electrode_response_sentence_cell(find(cell2mat([closing_type{2}])))]))-1;
    end
    a=cellfun(@(x) transpose(cell2mat(x)),node_closing,'UniformOutput',false);
    b=cellfun(@(x) transpose(x),node_closing_index_sentence,'UniformOutput',false);
    a=[a{:}];
    b=[b{:}];
    a=a-repmat(nanmean(a,1),size(a,1),1);
    colors=flipud(cbrewer('div','RdBu',10));
    
    % 
    x=1000/info.downsample_sampling_rate*   repmat([-(word_length-1):word_length]',1,size(a,2));
    figure(fix((i-1)/total_plots)+1);
   set(gcf,'position',[1822,-156,1257,1269])
    sub_title=sprintf('subj: %s,\n electrodes: %d ', info.subject ,language_electrode_num(i));
    ax=subplot(num_rows,num_columns,i-total_plots*fix((i-1)/total_plots));
    colors=flipud(cbrewer('seq','OrRd',length(unique(b))+3));
    for k=unique(b)
        ind=find(b==k);
        bl = boundedline(double(nanmean(x(:,ind),2)'), double(nanmean(a(:,ind),2)'), double(nanstd(a(:,ind),0,2)')./sqrt(length(ind)), 'cmap', colors(k+1,:),'alpha',ax,'transparency',.0);
        bl.LineWidth=2;
        bl.DisplayName=num2str(k);    
    end 
    hold on 
    h1=plot(double(nanmean(x,2)),0*double(nanmean(x,2)),'k--');hAnnotation = get(h1,'Annotation');hLegendEntry = get(hAnnotation,'LegendInformation');set(hLegendEntry,'IconDisplayStyle','off');
    h1=plot([0,0],get(gca,'ylim'),'k--');hAnnotation = get(h1,'Annotation');hLegendEntry = get(hAnnotation,'LegendInformation');set(hLegendEntry,'IconDisplayStyle','off');
    xlim([min(x(:)),max(x(:))])
    set(ax,'ydir', 'normal','box','off');
    title(sub_title);
    if ~mod(i,total_plots) | i==length(language_electrode_num)
    legend('show','Location','northeastoutside')
    end 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  compare closings for sentences and word lists 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all 
look_out_window=90; %200 ms
offset=30; % 100 ms
word_range=[1:8];
num_rows=5;
num_columns=3;
total_plots=num_rows*num_columns;
p=0;
colors=flipud(cbrewer('div','RdBu',10));
for i=1:length(language_electrode_num)
    % extract sentences 
    electrode_sentence_response=squeeze(session_sentence_hilbert_band_envelope_lang_elec_tensor(i,:,:));
    electrode_response_sentence_cell=mat2cell(electrode_sentence_response,ones(1,size(electrode_sentence_response,1)),size(electrode_sentence_response,2)/8*ones(1,8));
    % extract words 
    electrode_words_response=squeeze(session_words_hilbert_band_envelope_lang_elec_tensor(i,:,:));
    electrode_response_words_cell=mat2cell(electrode_words_response,ones(1,size(electrode_words_response,1)),size(electrode_words_response,2)/8*ones(1,8));
    %
    node_closing_sentence={};node_closing_index_sentence={};
    node_closing_words={};
    for k=1:8
        closing_type=total_closing(k,:);
        node_closing_sentence{k,1}=[ electrode_response_sentence_cell(find(cell2mat([closing_type{2}]))),electrode_response_sentence_cell(find(cell2mat([closing_type{3}])))];
        node_closing_words{k,1}=[ electrode_response_words_cell(find(cell2mat([closing_type{2}]))),electrode_response_words_cell(find(cell2mat([closing_type{3}])))];
        node_closing_index_sentence{k,1}=k*ones(size([electrode_response_sentence_cell(find(cell2mat([closing_type{2}])))]))-1;
    end
    % sentence
    closing_means_sentence=cellfun(@(x) mat2cell(transpose(repmat(nanmean(transpose(cell2mat(x)),1),size(x,2),1)),ones(1,size(x,1)),ones(1,size(x,2))),node_closing_sentence,'UniformOutput',false);
    normalize_node_closing_sentence=cellfun(@(x,y) cellfun(@(p,q) p-q,x,y,'uniformoutput',false),node_closing_sentence,closing_means_sentence,'UniformOutput',false);
    mean_val_sentence=cellfun(@(x) cellfun(@(y) nanmean(y(offset+[1:look_out_window])),x),normalize_node_closing_sentence,'UniformOutput',false);
    mean_val_sentence=cellfun(@(x) transpose(x),mean_val_sentence,'UniformOutput',false);
    closing_idx_val=cellfun(@(x) transpose(x),node_closing_index_sentence,'UniformOutput',false);
    mean_val_sentence=[mean_val_sentence{:}];closing_idx_val=[closing_idx_val{:}];
    % words 
    closing_means_words=cellfun(@(x) mat2cell(transpose(repmat(nanmean(transpose(cell2mat(x)),1),size(x,2),1)),ones(1,size(x,1)),ones(1,size(x,2))),node_closing_words,'UniformOutput',false);
    normalize_node_closing_words=cellfun(@(x,y) cellfun(@(p,q) p-q,x,y,'uniformoutput',false),node_closing_words,closing_means_words,'UniformOutput',false);
    mean_val_words=cellfun(@(x) cellfun(@(y) nanmean(y(offset+[1:look_out_window])),x),normalize_node_closing_words,'UniformOutput',false);
    mean_val_words=cellfun(@(x) transpose(x),mean_val_words,'UniformOutput',false);
    mean_val_words=[mean_val_words{:}];
    % 
    figure(fix((i-1)/total_plots)+1);
    set(gcf,'position',[ 1442,-129,1068,1269])
    sub_title=sprintf('subj: %s,\n electrodes: %d ', info.subject ,language_electrode_num(i));
    ax=subplot(num_rows,num_columns,i-total_plots*fix((i-1)/total_plots));
    colors=flipud(cbrewer('seq','OrRd',length(unique(b))+3));
    means_sentence=[];stds_sentence=[];
    means_words=[];stds_words=[];
    
    for k=unique(closing_idx_val)
        ind=find(closing_idx_val==k);
        means_sentence =[means_sentence;double(nanmean(mean_val_sentence(:,ind),2)')];
        stds_sentence=[stds_sentence;double(nanstd(mean_val_sentence(:,ind),0,2)')./sqrt(length(ind))];
        means_words =[means_words;double(nanmean(mean_val_words(:,ind),2)')];
        stds_words=[stds_words;double(nanstd(mean_val_words(:,ind),0,2)')./sqrt(length(ind))];
    end 
    % sentences 
    x=repmat([1,3],size(means_sentence,1),1);
    bl=bar(x',means_sentence');
    arrayfun(@(x,y) set(x,'FaceColor',colors(y,:)),bl,[1:size(bl,2)]);
    arrayfun(@(x,y) set(x,'Displayname',num2str(y)),bl,[1:size(bl,2)]-1);
    arrayfun(@(x,y) set(x,'EdgeColor',[0,0,0]),bl);
    arrayfun(@(x) set(x,'Barwidth',.4),bl);
    hold on
    bar_location=transpose(arrayfun(@(x) get(x,'XOffset'),bl));
    e=errorbar(repmat(bl(1).XData,size(bl,2),1)+bar_location,means_sentence,stds_sentence,'k.');
    arrayfun(@(x) set(e,'Capsize',0),e);arrayfun(@(x) set(e,'Linewidth',1.5),e);
    hAnnotation=arrayfun(@(x) get(x,'Annotation'),e);hLegendEntry = arrayfun(@(x) get(x,'LegendInformation'),hAnnotation);arrayfun(@(x) set(x,'IconDisplayStyle','off'),hLegendEntry);
    h1=plot(get(gca,'xlim'),0*get(gca,'ylim'),'k-');hAnnotation = get(h1,'Annotation');hLegendEntry = get(hAnnotation,'LegendInformation');set(hLegendEntry,'IconDisplayStyle','off');
    % words 
    ax_offset=unique(diff(arrayfun(@(x) get(x,'XOffset'),bl)));
    x=repmat([1,3]+ax_offset(1)/2,size(means_sentence,1),1);
    bl1=bar(x',means_words');
    arrayfun(@(x,y) set(x,'FaceColor',[1,1,1]),bl1,[1:size(bl1,2)]);
    arrayfun(@(x,y) set(x,'Displayname',num2str(y)),bl1,[1:size(bl1,2)]-1);
    arrayfun(@(x,y) set(x,'EdgeColor',colors(y,:)),bl1,[1:size(bl1,2)]);
    arrayfun(@(x) set(x,'Barwidth',.4),bl1);

    e1=errorbar(transpose(repmat(bl1(1).XData,size(bl1,2),1))+transpose(repmat(bar_location,1,2)),means_words',stds_words','.');
    arrayfun(@(x,y) set(x,'Color',colors(y,:)),e1,[1:size(bl1,2)]);
    arrayfun(@(x) set(e1,'Capsize',0),e1);arrayfun(@(x) set(e1,'Linewidth',1.5),e);
    hAnnotation=arrayfun(@(x) get(x,'Annotation'),e1);hLegendEntry = arrayfun(@(x) get(x,'LegendInformation'),hAnnotation);arrayfun(@(x) set(x,'IconDisplayStyle','off'),hLegendEntry);
    h1=plot(get(gca,'xlim'),0*get(gca,'ylim'),'k-');hAnnotation = get(h1,'Annotation');hLegendEntry = get(hAnnotation,'LegendInformation');set(hLegendEntry,'IconDisplayStyle','off');
    % 
    ax.XAxis.Visible = 'off'; 
    set(ax,'ydir', 'normal','box','off');
    set(ax,'xticklabel',{'pre','poset'});
    title(sub_title);
    if ~mod(i,total_plots) | i==length(language_electrode_num)
        legend('show','Location','northeastoutside')
        p=p+1;
        if ~exist(strcat(analysis_path,info.subject))
            mkdir(strcat(analysis_path,info.subject))
        end 
        window=sprintf('%d_%dms',floor(min(1000/info.downsample_sampling_rate*(offset+[1:look_out_window]))),...
            floor(max(1000/info.downsample_sampling_rate*(offset+[1:look_out_window]))));
        print(gcf, '-depsc', strcat(analysis_path,info.subject,'/',info.subject,'_node_closing_vs_words',num2str(p),'_window_',window,'.eps')); 
    end 
    
end





