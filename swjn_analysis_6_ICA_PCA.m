%% STEP 0: prepare the workspace
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
data_path='/Users/eghbalhosseiniasl1/MyData/ecog-sentence';
subject_id='AMC';
d= dir([data_path,strcat('/**/',subject_id,'*_crunched.mat')]);
fprintf(' %d .mat files were found \n', length(d));
gamma_band_index=4;

%
subject_ids=table2cell(unique(cell2table(cellfun(@(x) x(regexpi(x,'AMC')+[0:5]), {d.name},'UniformOutput',false)')));
%

save_path='/Users/eghbalhosseiniasl1/MyData/ECoG-sentence/crunched/';
analysis_path='/Users/eghbalhosseiniasl1/MyData/ECoG-sentence/analysis/swjn_analysis_6_ICA';
%
if 1
    fprintf('adding basic ecog tools to path \n');
    addpath('~/MyCodes/basic-ecog-tools/');
    addpath(genpath('~/MyCodes/ecog-sentence/'));
    addpath(genpath('~/MyCodes/basic-ecog-tools/activeBrain'));
    addpath(genpath('~/MyCodes/basic-ecog-tools/ecog-filters'));
    addpath(genpath('~/MyCodes/basic-ecog-tools/mex'));
end


find_index= @(input_cell) find(cell2mat(cellfun(@(x) ~isempty(x), input_cell,'UniformOutput',false)));
electrode_group='language';
do_print=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  STEP 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
all_sub_dat=[];

for k=1:length(subject_ids)-1
    session_sentence_examples={};
    session_wordlist_examples={};
    session_nonwords_examples={};
    session_jabberwocky_examples={};
    session_sentence_hilbert_band_envelope_tensor=[];
    session_words_hilbert_band_envelope_tensor=[];
    session_nonwords_hilbert_band_envelope_tensor=[];
    session_jabberwocky_hilbert_band_envelope_tensor=[];
    
    d_sub=(find(~cellfun(@isempty,cellfun(@(x) strfind(x,subject_ids{k}),{d.name},'UniformOutput',false))));
    for i=d_sub
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
        hilbert_band_envelope=cellfun(@(x) x(1:8),{sentences.signal_hilbert_zs_downsample_parsed},'UniformOutput',false);
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
    sub.sub_id=subject_ids{k};
    % signal
    sub.sess_sentence_hilb_tensor=session_sentence_hilbert_band_envelope_tensor;
    sub.sess_wordlist_hilb_tensor=session_words_hilbert_band_envelope_tensor;
    sub.sess_nonwords_hilb_tensor=session_nonwords_hilbert_band_envelope_tensor;
    sub.sess_jabberwocky_hilb_tensor=session_jabberwocky_hilbert_band_envelope_tensor;
    % stimuli
    sub.sess_sentence_stim=session_sentence_examples;
    sub.sess_wordlist_stim=session_wordlist_examples;
    sub.sess_nonwords_stim=session_nonwords_examples;
    sub.sess_jabberwocky_stim=session_jabberwocky_examples;
    sub.sess_stim_length=word_length;
    all_sub_dat=[all_sub_dat;sub];
    
    
    
end
%% find shared stimuli
stim_type={'sentence','wordlist','nonwords','jabberwocky'};
for k=1:size(stim_type,2)
    fprintf(['adding ',stim_type{k},'\n'])
    eval(strcat('A=all_sub_dat(1).sess_',stim_type{k},'_stim;')) ;
    for kk=1:size(all_sub_dat,1)
        eval(strcat('B=all_sub_dat(',num2str(kk),').sess_',stim_type{k},'_stim;')) ;
        A=intersect(A,B);
    end
    for kk=1:size(all_sub_dat,1)
        eval(strcat('all_sub_dat(kk).sess_',stim_type{k},'_stim_shared=A;')) ;
    end
end



%%
sub_ica_dat=struct;
find_first_index= @(input_cell,val) find(cell2mat(cellfun(@(x) ~isempty(x), input_cell,'UniformOutput',false)),val,'first');
stim_type={'sentence','wordlist','nonwords','jabberwocky'};
word_win=[1:135];% [0,450]
for k=1:size(stim_type,2)
    fprintf(['adding ',stim_type{k},' electrode:stim\n'])
    elec_stim_mat=[];
    elec_stim_str=[];
    elec_stim_wordpos=[];
    for kk=1:size(all_sub_dat,1)
        eval(strcat('A=all_sub_dat(',num2str(kk),').sess_',stim_type{k},'_stim;'));% get_stim
        eval(strcat('B=all_sub_dat(',num2str(kk),').sess_',stim_type{k},'_stim_shared;'));% get_stim_shared
        idx=cellfun(@(x) (regexpi(A,x)),B,'UniformOutput',false);
        idx_1st=cellfun(@(x) find_first_index(x,1), idx, 'UniformOutput',false);
        idx=cellfun(@(x) find_index(x), idx, 'UniformOutput',false);
        % get hilbert
        eval(strcat('Hilb=all_sub_dat(',num2str(kk),').sess_',stim_type{k},'_hilb_tensor;'));% get_stim_shared
        eval(strcat('Stim_L=all_sub_dat(',num2str(kk),').sess_stim_length;'));% get_word_length
        
        El=cellfun(@(x) double(squeeze(Hilb(:,x,:))),idx_1st','UniformOutput',false);
        
        El_demeaned=cellfun(@(x) transpose(bsxfun(@minus, x', mean(x'))),El,'UniformOutput',false);
        El_demeaned=El;
        El_cut_demeaned=[];
        El_cut=[];
        a=size(El_demeaned{1},1);
        for p=1:size(El_demeaned,2)
            El_cut=[El_cut,mat2cell(El{p},a,Stim_L*ones(1,8))];
            El_cut_demeaned=[El_cut_demeaned,mat2cell(El_demeaned{p},a,Stim_L*ones(1,8))];
        end
        El_cut_ave=cell2mat(cellfun(@(x) nanmean([x(:,word_win)],2),El_cut,'UniformOutput',false));
        %El_cut_ave=cell2mat(cellfun(@(x) nanmean(x,2),El_cut,'UniformOutput',false));
        %El_cut_demeaned_ave=cell2mat(cellfun(@(x) nanmean(x,2),El_cut_demeaned,'UniformOutput',false));
        El_cut_demeaned_ave=cell2mat(cellfun(@(x) nanmean([x(:,word_win)],2),El_cut_demeaned,'UniformOutput',false));
        El_string=cellfun(@(x) strsplit(A{x},' '),idx_1st','UniformOutput',false);
        El_wordpos=cellfun(@(x) 1:size(x,2),El_string,'UniformOutput',false);
        El_string=[El_string{:}];
        El_wordpos=[El_wordpos{:}];
        elec_stim_mat=[elec_stim_mat;El_cut_demeaned_ave];
        elec_stim_str=[elec_stim_str;El_string];
        elec_stim_wordpos=[elec_stim_wordpos;El_wordpos];
        
        %         figure
        %         for i=1:size(El_cut_ave,1)
        %         a=reshape(El_cut_demeaned_ave(i,:)',8,[]);
        %         subplot(7,7,i)
        %         plot(nanmean(a,2))
        %         hold on
        %         axis tight
        %         end
        %
    end
    eval(strcat('sub_ica_dat.sess_',stim_type{k},'_elec_stim_mat=elec_stim_mat;')) ;
    eval(strcat('sub_ica_dat.sess_',stim_type{k},'_elec_stim_str=elec_stim_str;')) ;
    eval(strcat('sub_ica_dat.sess_',stim_type{k},'_elec_stim_wordpos=elec_stim_wordpos;')) ;
end
%% 
subs_sentence_mat.sentence_matrix=sub_ica_dat.sess_sentence_elec_stim_mat;
subs_sentence_mat.sentence_str=sub_ica_dat.sess_sentence_elec_stim_str(1,:)
save(strcat(analysis_path,'/subj_sentence_mat.mat'),'subs_sentence_mat')
%%
condition='sentence';
switch condition
    case 'sentence'
        ica_mat=sub_ica_dat.sess_sentence_elec_stim_mat';
        stim_str=sub_ica_dat.sess_sentence_elec_stim_str(1,:);
    case 'nonwords'
        ica_mat=sub_ica_dat.sess_nonwords_elec_stim_mat';
        stim_str=sub_ica_dat.sess_nonwords_elec_stim_str(1,:);
    case 'wordlist'
        ica_mat=sub_ica_dat.sess_wordlist_elec_stim_mat';
        stim_str=sub_ica_dat.sess_wordlist_elec_stim_str(1,:);
    case 'jabberwocky'
        ica_mat=sub_ica_dat.sess_jabberwocky_elec_stim_mat';
        stim_str=sub_ica_dat.sess_jabberwocky_elec_stim_str(1,:);
end
close all
R_width=.55;

num_rows=5;
plot_dist=.03;
R_height=(.9-.05-plot_dist*num_rows)/num_rows;
R_start=R_height*(0:num_rows-1)+plot_dist*(1:num_rows)+.02;
num_columns=1;
total_plots=length(R_start);
p=0;
% dimensionality of data and components
M = 98; % number of features (e.g. sounds)
N = 416; % number of measures (e.g. fMRI voxels)
K = 10; % number of components
for K=5
    % create the data matrix
    % decomposition analysis
    N_RANDOM_INITS = 20;
    PLOT_FIGURES = 0;
    RAND_SEED = 1;
    [R_inferred, W_inferred,component_var] = nonparametric_ica(ica_mat, K, N_RANDOM_INITS, PLOT_FIGURES, RAND_SEED);
    fig1=figure;
    set(gcf,'position',[1376 353 988 992])
    axes('position',[.1,.6,.5,.3])
    colors = cbrewer('div', 'RdBu', 128);
    colors = flipud(colors); % puts red on top, blue at the bottom
    colormap(colors);
    imagesc(sub_ica_dat.sess_sentence_elec_stim_mat)
    set(gca, 'ydir', 'reverse','box','off');
    xlabel(condition)
    ylabel('Electrodes')
    title([condition,' Response'])
    %
    axes('position',[.1,.1,.5,.2])
    colors = cbrewer('div', 'RdBu', 128);
    colors = flipud(colors); % puts red on top, blue at the bottom
    colormap(colors);
    imagesc(W_inferred)
    set(gca, 'ydir', 'reverse','box','off');
    title('Weights')
    %
    axes('position',[.65,.1,.3,.4])
    colors = cbrewer('div', 'RdBu', 128);
    colors = flipud(colors); % puts red on top, blue at the bottom
    colormap(colors);
    imagesc(R_inferred)
    set(gca, 'ydir', 'reverse','box','off');
    title('Features');
    axes('position',[.65,.6,.3,.1])
    bl=bar(component_var,'Facecolor','flat');
    set(gca,'box','off')
    ylabel('%var explained');
    xlabel('component');
    x_lim=get(gca,'xlim');
    axes('position',[.65,.75,.3,.1])
    plot(cumsum(component_var),'k','Marker','o')
    set(gca,'xlim',x_lim);
    set(gca,'box','off');
    ylabel('%total var explained');
    if ~exist(strcat(analysis_path))
        mkdir(strcat(analysis_path))
    end
    coeff=450/135;
    if do_print==1
    print(fig1,'-bestfit','-painters', '-dpdf', strcat(analysis_path,'/','ECoG_ICA_',condition,'_num_components_',num2str(K),'_window_'...
        ,sprintf('%1d-%1d',coeff*(min(word_win)-1),coeff*max(word_win)),'.pdf'));
    % reorder based on wordpos
    end 
    word_pos=sub_ica_dat.sess_sentence_elec_stim_wordpos(1,:);
    sort_idx={};
    unique_word_pos=unique(word_pos);
    for i=1:length(unique(word_pos))
        sort_idx=[sort_idx,find(word_pos==unique_word_pos(i))];
    end
    ax_min=floor(min(R_inferred(:)));
    ax_max=ceil(max(R_inferred(:)));
    for i=1:K
        R_indx=i-num_rows*fix((i-1)/num_rows);
        f=figure(fix((i-1)/num_rows)+2);
        set(gcf,'position',[29,12,852,1270])
        colors = cbrewer('div', 'RdYlBu', 8);
        ax=axes('position',[.05,R_start(R_indx),R_width,R_height]);
        A=cellfun(@(x) transpose(R_inferred(x,i)),sort_idx,'UniformOutput',false);
        R_inferred_sort=cellfun(@(x) sort(x),A,'UniformOutput',false);
        R_sort_ave=cellfun(@mean,A);
        bl=bar(cell2mat(R_inferred_sort),'Facecolor','flat');
        bar_color=colors(word_pos(cell2mat(sort_idx)),:);
        set(bl,'Linestyle','none');
        set(bl,'Displayname','sentences')
        bl.CData=bar_color;
        set(gca,'box','off')
        ax.YLim=[ax_min,ax_max];
        ax.XAxis.Visible='off';
        if i==1
            xlabel([condition,': Feature dim sorted by word position']);
            ax.XAxis.Visible='on';
            ax.XAxis.FontWeight='bold'
        end
        ax.Title.String=[sprintf('component %d ',i),sprintf('; var explained: %%%1.1f',component_var(i))]
        %
        ax1=axes('position',[.1+R_width,R_start(R_indx),.2,R_height]);
        hold on;
        bl1=arrayfun(@(x,y) bar(x,y,'Facecolor','flat'),unique_word_pos,R_sort_ave);
        set(bl1,'Linestyle','none');
        arrayfun(@(x,y) set(x,'CData',colors(y,:)),bl1,unique_word_pos);
        arrayfun(@(x,y) set(x,'Displayname',num2str(y)),bl1,unique_word_pos);
        
        set(gca,'box','off')
        ax1.XAxis.Visible='off';
        if i==1
            ax1.Title.String={'average over word position'};
            legend('show','position',[[.35+R_width,R_start(R_indx),.05,R_height/1.5]])
        end
        set(gca,'box','off')
        
        if ~mod(i,total_plots) | i==K
            %legend('show','Location','northeastoutside')
            p=p+1;
            if ~exist(strcat(analysis_path))
                mkdir(strcat(analysis_path))
            end
            if do_print==1
                set(gcf,'PaperPosition',[.25 .25 8 6])
                print(f,'-bestfit', '-painters','-dpdf', strcat(analysis_path,'/','ECoG_ICA_',condition,'_num_components_',num2str(K),...
                    '_window_'...
                    ,sprintf('%1d-%1d',coeff*(min(word_win)-1),coeff*max(word_win)),'_fig_',num2str(fix((i-1)/num_rows)+2),'.pdf'));
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plot top 5 and bottom loadings for each components
    for i=1:K
        R_indx=i-num_rows*fix((i-1)/num_rows);
        f=figure(5+fix((i-1)/num_rows)+2);
        set(gcf,'position',[29,12,852,1270])
        colors = cbrewer('div', 'RdYlBu', size(stim_str,2));
        colors= flipud(colors);
        colors_top=cbrewer('seq','Oranges',size(stim_str,2));
        %colors_top= flipud(colors_top);
        colors_bottom=cbrewer('seq','Blues',size(stim_str,2));
        colors_bottom= flipud(colors_bottom);
        ax=axes('position',[R_start(R_indx),0.04,.1,.92]);
        
        [A_order,idx]=sort(R_inferred(:,i));
        idx_top=find(A_order>median(A_order));
        idx_bottom=find(A_order<=median(A_order));
        A_top=A_order(idx_top);
        stim_top=stim_str(idx(idx_top));
        word_idx_top=word_pos(idx(idx_top));
        A_bottom=A_order(idx_bottom);
        stim_bottom=stim_str(idx(idx_bottom));
        word_idx_bottom=word_pos(idx(idx_bottom));
        hold on
        % top
        offset=1;
        A=A_top;
        x=[1:length(A)]';
        bl=arrayfun(@(x,y) plot(offset+[0,x],[y,y],'k','LineWidth',2),A,x );
        XLim_top=[floor(min(A)),ceil(max(A))];
        arrayfun(@(x) set(bl(x),'color',colors_top(idx_top(x),:)),x);
        
        tx=arrayfun(@(x,y,z) text(offset+.3+x,y,stim_top(y)),A,x );
        arrayfun(@(x) set(tx(x),'fontsize',4),x);
        plot(0*x+offset,x,'k','LineWidth',.5);
        % bottom
        A=A_bottom;
        x=[1:length(A)]';
        bl1=arrayfun(@(x,y) plot(-offset+[0,x],[y,y],'k','LineWidth',2),A,x );
        arrayfun(@(x) set(bl1(x),'color',colors_bottom(idx_bottom(x),:)),x);
        tx=arrayfun(@(x,y,z) text(-offset-.3+x,y,stim_bottom(y),'HorizontalAlignment','right'),A,x );
        arrayfun(@(x) set(tx(x),'fontsize',3.5),x);
        plot(0*x-offset,x,'k','LineWidth',.5);
        XLim_bottom=[floor(min(A)),ceil(max(A))];
        ax.Title.String={sprintf('component %d ',i),sprintf('var explained: %%%1.1f',component_var(i))};
        ax.Title.FontSize=8;
        %
        ax.YLim=[-2,max([length(bl),length(bl1)])];
        
        % fix x axis
        ax.XAxis.Visible='off';
        ax.YAxis.Visible='off';
        set(gca,'box','off');
        % top
        plot(offset+XLim_top,[0,0],'k','LineWidth',.5)
        arrayfun(@(x) plot([x,x],[0,-.5],'k','LineWidth',.5),offset+XLim_top)
        arrayfun(@(x,y) text(x+offset,-.5,num2str(x),'VerticalAlignment','top','HorizontalAlignment','center','fontsize',7),XLim_top);
        text(mean(offset+XLim_top),-2,{'above', 'median'},'HorizontalAlignment','center','VerticalAlignment','top','fontsize',7)
        % bottom
        plot(-offset+XLim_bottom,[0,0],'k','LineWidth',.5)
        arrayfun(@(x) plot([x,x],[0,-.5],'k','LineWidth',.5),-offset+XLim_bottom)
        arrayfun(@(x,y) text(-offset+x,-.5,num2str(x),'VerticalAlignment','top','HorizontalAlignment','center','fontsize',7),XLim_bottom);
        text(mean(-offset+XLim_bottom),-2,{'below', 'median'},'HorizontalAlignment','center','VerticalAlignment','top','fontsize',7)
        
        %
        
        if ~mod(i,total_plots) | i==K
            %legend('show','Location','northeastoutside')
            p=p+1;
            if ~exist(strcat(analysis_path))
                mkdir(strcat(analysis_path))
            end
            if do_print==1
                set(gcf,'PaperPosition',[.25 .25 8 6])
                print(f,'-painters','-fillpage', '-dpdf', strcat(analysis_path,'/','ECoG_ICA_',condition,'_num_components_',num2str(K),...
                    '_component_loading_fig_',num2str(fix((i-1)/num_rows)+2),'_window_'...
                    ,sprintf('%1d-%1d',coeff*(min(word_win)-1),coeff*max(word_win)),'.pdf'));
            end
        end
    end
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% do the first top component and last 5
    for i=1:K
        R_indx=i-num_rows*fix((i-1)/num_rows);
        f=figure(10+fix((i-1)/num_rows)+2);
        set(gcf,'position',[29,12,852,1270])
        colors = cbrewer('div', 'RdYlBu', size(stim_str,2));
        colors= flipud(colors);
        colors_top=cbrewer('seq','Oranges',size(stim_str,2));
        %colors_top= flipud(colors_top);
        colors_bottom=cbrewer('seq','Blues',size(stim_str,2));
        colors_bottom= flipud(colors_bottom);
        ax=axes('position',[R_start(R_indx),0.04,.1,.92]);
        
        R_component=R_inferred(:,i);
        [A_order,idx]=sort(R_component);
        idx_sub=flipud(idx(end-4:end));
        a=arrayfun(@(x) [1:8]-x, word_pos(idx_sub)','UniformOutput',false);
        max_pos=find([a{:}]==0);
        sent_pos=cellfun(@(x,y) y+x,a,mat2cell(idx_sub,ones(size(idx_sub))),'UniformOutput',false);
        sent=cellfun(@(x) stim_str(x),sent_pos,'UniformOutput',false);
        val=cellfun(@(x) R_component(x)',sent_pos,'UniformOutput',false);
        y_val_1=arrayfun(@(x) [1:8]+8*x+2*x, [0:4]','UniformOutput',false);
        y_val=fliplr([y_val_1{:}]);
        val=([val{:}]);
        %
        hold on
        %offset=abs(min(A_order))+abs(min(A_order))+1;
        offset=max(abs(A_order));
        bl=arrayfun(@(x,y) plot(offset+[0,x],[y,y],'k','LineWidth',10),val,y_val );
        XLim_top=[floor(min(val)),ceil(max(val))];
        colors=cbrewer('seq','Oranges',size(a,2)+60);
        colors= flipud(colors);
        arrayfun(@(x) set(bl(x),'color',colors(x,:)),1:length(bl));
        a=[sent{:}];
        text_x=val;
        text_x(text_x<0)=0;
        tx=arrayfun(@(x,y,z) text(x+offset+.5,y,a(z),'HorizontalAlignment','left'),text_x,y_val,1:length(bl) );
        arrayfun(@(x) set(tx(x),'fontsize',4),1:length(bl));
        arrayfun(@(x) set(tx(x),'fontweight','bold','fontsize',5,'fontangle','italic','color',[.5,0,0]),max_pos);
        %plot(0*y_val+offset,y_val,'k--','LineWidth',.5);
        cellfun(@(x) plot(0*x+offset,x,'k--','LineWidth',.5),y_val_1);
        
        % bottom
        
        idx_sub=flipud(idx(1:5));
        a=arrayfun(@(x) [1:8]-x, word_pos(idx_sub)','UniformOutput',false);
        max_pos=find([a{:}]==0);
        sent_pos=cellfun(@(x,y) y+x,a,mat2cell(idx_sub,ones(size(idx_sub))),'UniformOutput',false);
        sent=cellfun(@(x) stim_str(x),sent_pos,'UniformOutput',false);
        val=cellfun(@(x) R_component(x)',sent_pos,'UniformOutput',false);
        y_val=arrayfun(@(x) [1:8]+8*x+2*x, [0:4]','UniformOutput',false);
        y_val=fliplr([y_val{:}]);
        val=([val{:}]);
        XLim_bottom=[floor(min(val)),ceil(max(val))];
        bl=arrayfun(@(x,y) plot(-offset+[0,x],[y,y],'k','LineWidth',10),val,y_val );
        
        colors=cbrewer('seq','Blues',size(a,2)+100);
        %colors= flipud(colors);
        
        arrayfun(@(x) set(bl(length(bl)-x+1),'color',colors(end-x,:)),1:length(bl));
        a=[sent{:}];
        text_x=val;
        text_x(text_x>0)=0;
        tx=arrayfun(@(x,y,z) text(x-offset-.5,y,a(z),'HorizontalAlignment','right'),text_x,y_val,1:length(bl) );
        arrayfun(@(x) set(tx(x),'fontsize',4),1:length(bl));
        arrayfun(@(x) set(tx(x),'fontweight','bold','fontsize',5,'fontangle','italic','color',[.5,0,0]),max_pos);
        cellfun(@(x) plot(0*x-offset,x,'k--','LineWidth',.5),y_val_1);
        %
        ax.Title.String=sprintf('component %d ; %%%1.1f ',i,component_var(i));
        ax.Title.Position(2)=49;
        ax.Title.FontSize=8;
        %
        ax.YLim=[-2,max(y_val)];
        
        % fix x axis
        ax.XAxis.Visible='off';
        ax.YAxis.Visible='off';
        set(gca,'box','off');
        % top
        plot(offset+XLim_top,[0,0],'k','LineWidth',.5);
        arrayfun(@(x) plot([x,x],[0,-.1],'k','LineWidth',.1),offset+XLim_top);
        arrayfun(@(x,y) text(x+offset,-.5,num2str(x),'VerticalAlignment','top','HorizontalAlignment','center','fontsize',7),XLim_top);
        text(mean(offset+XLim_top),-1,{'top'},'HorizontalAlignment','center','VerticalAlignment','top','fontsize',7);
        % bottom ;
        plot(-offset+XLim_bottom,[0,0],'k','LineWidth',.5);
        arrayfun(@(x) plot([x,x],[0,-.1],'k','LineWidth',.5),-offset+XLim_bottom);
        arrayfun(@(x,y) text(-offset+x,-.5,num2str(x),'VerticalAlignment','top','HorizontalAlignment','center','fontsize',7),XLim_bottom);
        text(mean(-offset+XLim_bottom),-1,{'bottom'},'HorizontalAlignment','center','VerticalAlignment','top','fontsize',7);
        %
        
        if ~mod(i,total_plots) | i==K
            %legend('show','Location','northeastoutside')
            p=p+1;
            if ~exist(strcat(analysis_path))
                mkdir(strcat(analysis_path));
            end
            if do_print==1
                set(gcf,'PaperPosition',[.25 .25 8 6])
                print(f,'-painters','-fillpage', '-dpdf', strcat(analysis_path,'/','ECoG_ICA_',condition,'_num_components_',num2str(K),...
                    '_component_top_bottom_loading_',num2str(fix((i-1)/num_rows)+2),'_window_'...
                    ,sprintf('%1d-%1d',coeff*(min(word_win)-1),coeff*max(word_win)),'.pdf'));
            end
        end
    end
    
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% display PCA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub_pca_dat=sub_ica_dat;
condition='sentence';
switch condition
    case 'sentence'
        pca_mat=sub_pca_dat.sess_sentence_elec_stim_mat';
        stim_str=sub_pca_dat.sess_sentence_elec_stim_str(1,:);
    case 'nonwords'
        pca_mat=sub_pca_dat.sess_nonwords_elec_stim_mat';
        stim_str=sub_pca_dat.sess_nonwords_elec_stim_str(1,:);
    case 'wordlist'
        pca_mat=sub_pca_dat.sess_wordlist_elec_stim_mat';
        stim_str=sub_pca_dat.sess_wordlist_elec_stim_str(1,:);
    case 'jabberwocky'
        pca_mat=sub_pca_dat.sess_jabberwocky_elec_stim_mat';
        stim_str=sub_pca_dat.sess_jabberwocky_elec_stim_str(1,:);
end
close all
R_width=.55;
num_rows=5;
plot_dist=.03;
R_height=(.9-.05-plot_dist*num_rows)/num_rows;
R_start=R_height*(0:num_rows-1)+plot_dist*(1:num_rows)+.02;
num_columns=1;
total_plots=length(R_start);
p=0;
% dimensionality of data and components
M = 98; % number of features (e.g. sounds)
N = 416; % number of measures (e.g. fMRI voxels)
K = 10; % number of components
for K=5
    % create the data matrix
    % decomposition analysis
    N_RANDOM_INITS = 20;
    PLOT_FIGURES = 0;
    RAND_SEED = 1;
    %[R_inferred, W_inferred,component_var] = nonparametric_pca(pca_mat, K, N_RANDOM_INITS, PLOT_FIGURES, RAND_SEED);
    [coeff,score,latent,tsquared,explained,mu] = pca(pca_mat') ;
    R_inferred=coeff(:,1:K);
    W_inferred=transpose(score(:,1:K));
    component_var=explained(1:K)';
    fig1=figure;
    set(gcf,'position',[1376 353 988 992]);
    axes('position',[.1,.6,.5,.3]);
    colors = cbrewer('div', 'RdBu', 128);
    colors = flipud(colors); % puts red on top, blue at the bottom
    colormap(colors);
    imagesc(sub_pca_dat.sess_sentence_elec_stim_mat)
    set(gca, 'ydir', 'reverse','box','off');
    xlabel(condition)
    ylabel('Electrodes')
    title([condition,' Response'])
    %
    axes('position',[.1,.1,.5,.2])
    colors = cbrewer('div', 'RdBu', 128);
    colors = flipud(colors); % puts red on top, blue at the bottom
    colormap(colors);
    imagesc(W_inferred)
    set(gca, 'ydir', 'reverse','box','off');
    title('Weights')
    %
    axes('position',[.65,.1,.3,.4])
    colors = cbrewer('div', 'RdBu', 128);
    colors = flipud(colors); % puts red on top, blue at the bottom
    colormap(colors);
    imagesc(R_inferred)
    set(gca, 'ydir', 'reverse','box','off');
    title('Features');
    axes('position',[.65,.6,.3,.1])
    bl=bar(component_var,'Facecolor','flat');
    set(gca,'box','off')
    ylabel('%var explained');
    xlabel('component');
    x_lim=get(gca,'xlim');
    axes('position',[.65,.75,.3,.1])
    plot(cumsum(component_var),'k','Marker','o');
    set(gca,'xlim',x_lim);
    set(gca,'box','off');
    ylabel('%total var explained');
    if ~exist(strcat(analysis_path))
        mkdir(strcat(analysis_path));
    end
    coeff=450/135;
    if do_print==1 
    print(fig1,'-bestfit','-painters', '-dpdf', strcat(analysis_path,'/','ECoG_pca_',condition,'_num_components_',num2str(K),'_window_'...
        ,sprintf('%1d-%1d',coeff*(min(word_win)-1),coeff*max(word_win)),'.pdf'));
    end 
    % reorder based on wordpos
    word_pos=sub_pca_dat.sess_sentence_elec_stim_wordpos(1,:);
    sort_idx={};
    unique_word_pos=unique(word_pos);
    for i=1:length(unique(word_pos))
        sort_idx=[sort_idx,find(word_pos==unique_word_pos(i))];
    end
    ax_min=floor(min(R_inferred(:)));
    ax_max=ceil(max(R_inferred(:)));
    for i=1:K
        R_indx=i-num_rows*fix((i-1)/num_rows);
        f=figure(fix((i-1)/num_rows)+2);
        set(gcf,'position',[29,12,852,1270]);
        colors = cbrewer('div', 'RdYlBu', 8);
        ax=axes('position',[.05,R_start(R_indx),R_width,R_height]);
        A=cellfun(@(x) transpose(R_inferred(x,i)),sort_idx,'UniformOutput',false);
        R_inferred_sort=cellfun(@(x) sort(x),A,'UniformOutput',false);
        R_sort_ave=cellfun(@mean,A);
        bl=bar(cell2mat(R_inferred_sort),'Facecolor','flat');
        bar_color=colors(word_pos(cell2mat(sort_idx)),:);
        set(bl,'Linestyle','none');
        
        bl.CData=bar_color;
        set(gca,'box','off')
        %ax.YLim=[ax_min,ax_max];
        ax.XAxis.Visible='off';
        if i==1
            xlabel([condition,': Feature dim sorted by word position']);
            ax.XAxis.Visible='on';
            ax.XAxis.FontWeight='bold';
        end
        ax.Title.String=[sprintf('component %d ',i),sprintf('; var explained: %%%1.1f',component_var(i))];
        %
        ax1=axes('position',[.1+R_width,R_start(R_indx),.2,R_height]);
        hold on;
        bl1=arrayfun(@(x,y) bar(x,y,'Facecolor','flat'),unique_word_pos,R_sort_ave);
        set(bl1,'Linestyle','none');
        arrayfun(@(x,y) set(x,'CData',colors(y,:)),bl1,unique_word_pos);
        arrayfun(@(x,y) set(x,'Displayname',num2str(y)),bl1,unique_word_pos);
        
        set(gca,'box','off')
        ax1.XAxis.Visible='off';
        if i==1
            ax1.Title.String={'average over word position'};
            legend('show','position',[[.35+R_width,R_start(R_indx),.05,R_height/1.5]]);
        end
        set(gca,'box','off')
        
        if ~mod(i,total_plots) | i==K
            %legend('show','Location','northeastoutside')
            p=p+1;
            if ~exist(strcat(analysis_path))
                mkdir(strcat(analysis_path));
            end
            if do_print==1
                set(gcf,'PaperPosition',[.25 .25 8 6])
                print(f,'-bestfit', '-painters','-dpdf', strcat(analysis_path,'/','ECoG_pca_',condition,'_num_components_',num2str(K),...
                    '_window_'...
                    ,sprintf('%1d-%1d',coeff*(min(word_win)-1),coeff*max(word_win)),'_fig_',num2str(fix((i-1)/num_rows)+2),'.pdf'));
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plot top 5 and bottom loadings for each components
    for i=1:K
        R_indx=i-num_rows*fix((i-1)/num_rows);
        f=figure(5+fix((i-1)/num_rows)+2);
        set(gcf,'position',[29,12,852,1270]);
        colors = cbrewer('div', 'RdYlBu', size(stim_str,2));
        colors= flipud(colors);
        colors_top=cbrewer('seq','Oranges',size(stim_str,2));
        %colors_top= flipud(colors_top);
        colors_bottom=cbrewer('seq','Blues',size(stim_str,2));
        colors_bottom= flipud(colors_bottom);
        ax=axes('position',[R_start(R_indx),0.04,.1,.92]);
        
        [A_order,idx]=sort(R_inferred(:,i));
        idx_top=find(A_order>median(A_order));
        idx_bottom=find(A_order<=median(A_order));
        A_top=A_order(idx_top);
        stim_top=stim_str(idx(idx_top));
        word_idx_top=word_pos(idx(idx_top));
        A_bottom=A_order(idx_bottom);
        stim_bottom=stim_str(idx(idx_bottom));
        word_idx_bottom=word_pos(idx(idx_bottom));
        hold on
        % top
        offset=mean(abs(A_order));
        A=A_top;
        x=[1:length(A)]';
        bl=arrayfun(@(x,y) plot(offset+[0,x],[y,y],'k','LineWidth',2),A,x );
        XLim_top=[min(A),(max(A))];
        arrayfun(@(x) set(bl(x),'color',colors_top(idx_top(x),:)),x);
        
        tx=arrayfun(@(x,y,z) text(offset+x,y,stim_top(y)),A,x );
        arrayfun(@(x) set(tx(x),'fontsize',4),x);
        plot(0*x+offset,x,'k','LineWidth',.5);
        % bottom
        A=A_bottom;
        x=[1:length(A)]';
        bl1=arrayfun(@(x,y) plot(-offset+[0,x],[y,y],'k','LineWidth',2),A,x );
        arrayfun(@(x) set(bl1(x),'color',colors_bottom(idx_bottom(x),:)),x);
        tx=arrayfun(@(x,y,z) text(-offset+x,y,stim_bottom(y),'HorizontalAlignment','right'),A,x );
        arrayfun(@(x) set(tx(x),'fontsize',3.5),x);
        plot(0*x-offset,x,'k','LineWidth',.5);
        XLim_bottom=[min(A),max(A)];
        ax.Title.String={sprintf('component %d ',i),sprintf('var explained: %%%1.1f',component_var(i))};
        ax.Title.FontSize=8;
        %
        ax.YLim=[-2,max([length(bl),length(bl1)])];
        
        % fix x axis
        ax.XAxis.Visible='off';
        ax.YAxis.Visible='off';
        set(gca,'box','off');
        % top
        plot(offset+XLim_top,[0,0],'k','LineWidth',.5);
        arrayfun(@(x) plot([x,x],[0,-.5],'k','LineWidth',.5),offset+XLim_top);
        arrayfun(@(x,y) text(x+offset,-.5,sprintf('%1.1f',x),'VerticalAlignment','top','HorizontalAlignment','center','fontsize',7),XLim_top);
        text(mean(offset+XLim_top),-2,{'above', 'median'},'HorizontalAlignment','center','VerticalAlignment','top','fontsize',7);
        % bottom
        plot(-offset+XLim_bottom,[0,0],'k','LineWidth',.5);
        arrayfun(@(x) plot([x,x],[0,-.5],'k','LineWidth',.5),-offset+XLim_bottom);
        arrayfun(@(x,y) text(-offset+x,-.5,sprintf('%1.1f',x),'VerticalAlignment','top','HorizontalAlignment','center','fontsize',7),XLim_bottom);
        text(mean(-offset+XLim_bottom),-2,{'below', 'median'},'HorizontalAlignment','center','VerticalAlignment','top','fontsize',7);
        %
        if ~mod(i,total_plots) | i==K
            %legend('show','Location','northeastoutside')
            p=p+1;
            if ~exist(strcat(analysis_path))
                mkdir(strcat(analysis_path));
            end
            if do_print==1
                set(gcf,'PaperPosition',[.25 .25 8 6]);
                print(f,'-painters','-fillpage', '-dpdf', strcat(analysis_path,'/','ECoG_pca_',condition,'_num_components_',num2str(K),...
                    '_component_loading_fig_',num2str(fix((i-1)/num_rows)+2),'_window_'...
                    ,sprintf('%1d-%1d',coeff*(min(word_win)-1),coeff*max(word_win)),'.pdf'));
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% do the first top component and last 5
    for i=1:K
        R_indx=i-num_rows*fix((i-1)/num_rows);
        f=figure(10+fix((i-1)/num_rows)+2);
        set(gcf,'position',[29,12,852,1270]);
        colors = cbrewer('div', 'RdYlBu', size(stim_str,2));
        colors= flipud(colors);
        colors_top=cbrewer('seq','Oranges',size(stim_str,2));
        %colors_top= flipud(colors_top);
        colors_bottom=cbrewer('seq','Blues',size(stim_str,2));
        colors_bottom= flipud(colors_bottom);
        ax=axes('position',[R_start(R_indx)+i*.02,0.04,.1,.92]);
        
        R_component=R_inferred(:,i);
        [A_order,idx]=sort(R_component);
        idx_sub=flipud(idx(end-4:end));
        a=arrayfun(@(x) [1:8]-x, word_pos(idx_sub)','UniformOutput',false);
        max_pos=find([a{:}]==0);
        sent_pos=cellfun(@(x,y) y+x,a,mat2cell(idx_sub,ones(size(idx_sub))),'UniformOutput',false);
        sent=cellfun(@(x) stim_str(x),sent_pos,'UniformOutput',false);
        val=cellfun(@(x) R_component(x)',sent_pos,'UniformOutput',false);
        y_val_1=arrayfun(@(x) [1:8]+8*x+2*x, [0:4]','UniformOutput',false);
        y_val=fliplr([y_val_1{:}]);
        val=([val{:}]);
        %
        hold on
        offset=max(abs(A_order));
        bl=arrayfun(@(x,y) plot(offset+[0,x],[y,y],'k','LineWidth',10),val,y_val );
        XLim_top=[(min(val)),(max(val))];
        colors=cbrewer('seq','Oranges',size(a,2)+60);
        colors= flipud(colors);
        arrayfun(@(x) set(bl(x),'color',colors(x,:)),1:length(bl));
        a=[sent{:}];
        text_x=val;
        text_x(text_x<0)=0;
        tx=arrayfun(@(x,y,z) text(x+offset+.05,y,a(z),'HorizontalAlignment','left'),text_x,y_val,1:length(bl) );
        arrayfun(@(x) set(tx(x),'fontsize',4),1:length(bl));
        arrayfun(@(x) set(tx(x),'fontweight','bold','fontsize',5,'fontangle','italic','color',[.5,0,0]),max_pos);
        %plot(0*y_val+offset,y_val,'k--','LineWidth',.5);
        cellfun(@(x) plot(0*x+offset,x,'k--','LineWidth',.5),y_val_1);
        
        % bottom
        
        idx_sub=flipud(idx(1:5));
        a=arrayfun(@(x) [1:8]-x, word_pos(idx_sub)','UniformOutput',false);
        max_pos=find([a{:}]==0);
        sent_pos=cellfun(@(x,y) y+x,a,mat2cell(idx_sub,ones(size(idx_sub))),'UniformOutput',false);
        sent=cellfun(@(x) stim_str(x),sent_pos,'UniformOutput',false);
        val=cellfun(@(x) R_component(x)',sent_pos,'UniformOutput',false);
        y_val=arrayfun(@(x) [1:8]+8*x+2*x, [0:4]','UniformOutput',false);
        y_val=fliplr([y_val{:}]);
        val=([val{:}]);
        XLim_bottom=[(min(val)),(max(val))];
        bl=arrayfun(@(x,y) plot(-offset+[0,x],[y,y],'k','LineWidth',10),val,y_val );
        
        colors=cbrewer('seq','Blues',size(a,2)+100);
        %colors= flipud(colors);
        
        arrayfun(@(x) set(bl(length(bl)-x+1),'color',colors(end-x,:)),1:length(bl));
        a=[sent{:}];
        text_x=val;
        text_x(text_x>0)=0;
        tx=arrayfun(@(x,y,z) text(x-offset-.05,y,a(z),'HorizontalAlignment','right'),text_x,y_val,1:length(bl) );
        arrayfun(@(x) set(tx(x),'fontsize',4),1:length(bl));
        arrayfun(@(x) set(tx(x),'fontweight','bold','fontsize',5,'fontangle','italic','color',[.5,0,0]),max_pos);
        cellfun(@(x) plot(0*x-offset,x,'k--','LineWidth',.5),y_val_1);
        %
        ax.Title.String=sprintf('component %d ; %%%1.1f ',i,component_var(i));
        ax.Title.Position(2)=49;
        ax.Title.FontSize=8;
        %
        ax.YLim=[-2,max(y_val)];
        
        % fix x axis
        ax.XAxis.Visible='off';
        ax.YAxis.Visible='off';
        set(gca,'box','off');
        % top
        plot(offset+XLim_top,[0,0],'k','LineWidth',.5);
        arrayfun(@(x) plot([x,x],[0,-.1],'k','LineWidth',.1),offset+XLim_top);
        arrayfun(@(x,y) text(x+offset,-.5,sprintf('%1.1f',x),'VerticalAlignment','top','HorizontalAlignment','center','fontsize',7),XLim_top);
        text(mean(offset+XLim_top),-1,{'top'},'HorizontalAlignment','center','VerticalAlignment','top','fontsize',7);
        % bottom
        plot(-offset+XLim_bottom,[0,0],'k','LineWidth',.5);
        arrayfun(@(x) plot([x,x],[0,-.1],'k','LineWidth',.5),-offset+XLim_bottom);
        arrayfun(@(x,y) text(-offset+x,-.2,sprintf('%1.1f',x),'VerticalAlignment','top','HorizontalAlignment','center','fontsize',7),XLim_bottom);
        text(mean(-offset+XLim_bottom),-1,{'bottom'},'HorizontalAlignment','center','VerticalAlignment','top','fontsize',7);
        %
        
        if ~mod(i,total_plots) | i==K
            %legend('show','Location','northeastoutside')
            p=p+1;
            if ~exist(strcat(analysis_path))
                mkdir(strcat(analysis_path));
            end
            if do_print==1
                set(gcf,'PaperPosition',[.25 .25 8 6])
                print(f,'-painters','-bestfit', '-dpdf', strcat(analysis_path,'/','ECoG_pca_',condition,'_num_components_',num2str(K),...
                    '_component_top_bottom_loading_',num2str(fix((i-1)/num_rows)+2),'_window_'...
                    ,sprintf('%1d-%1d',coeff*(min(word_win)-1),coeff*max(word_win)),'.pdf'));
            end
        end
    end
    
    
end
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% schematic for presentation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub_pca_dat=sub_ica_dat;
condition='sentence';
switch condition
    case 'sentence'
        pca_mat=sub_pca_dat.sess_sentence_elec_stim_mat';
        stim_str=sub_pca_dat.sess_sentence_elec_stim_str(1,:);
    case 'nonwords'
        pca_mat=sub_pca_dat.sess_nonwords_elec_stim_mat';
        stim_str=sub_pca_dat.sess_nonwords_elec_stim_str(1,:);
    case 'wordlist'
        pca_mat=sub_pca_dat.sess_wordlist_elec_stim_mat';
        stim_str=sub_pca_dat.sess_wordlist_elec_stim_str(1,:);
    case 'jabberwocky'
        pca_mat=sub_pca_dat.sess_jabberwocky_elec_stim_mat';
        stim_str=sub_pca_dat.sess_jabberwocky_elec_stim_str(1,:);
end
close all
R_width=.55;
num_rows=5;
plot_dist=.03;
R_height=(.9-.05-plot_dist*num_rows)/num_rows;
R_start=R_height*(0:num_rows-1)+plot_dist*(1:num_rows)+.02;
num_columns=1;
total_plots=length(R_start);
p=0;
% dimensionality of data and components
M = 98; % number of features (e.g. sounds)
N = 416; % number of measures (e.g. fMRI voxels)
K = 10; % number of components
for K=5
    % create the data matrix
    % decomposition analysis
    N_RANDOM_INITS = 20;
    PLOT_FIGURES = 0;
    RAND_SEED = 1;
    %[R_inferred, W_inferred,component_var] = nonparametric_pca(pca_mat, K, N_RANDOM_INITS, PLOT_FIGURES, RAND_SEED);
    [coeff,score,latent,tsquared,explained,mu] = pca(pca_mat') ;
    R_inferred=coeff(:,1:K);
    W_inferred=transpose(score(:,1:K));
    component_var=explained(1:K)';
  
    fig1=figure;
    aspect_ration=9.32./4.13;
    y=500;
    set(fig1,'position',[591 455 aspect_ration*y y]);
    ax=axes('position',[.1,.6,.3,.3]);
    colors = cbrewer('div', 'RdBu', 128);
    colors = flipud(colors); % puts red on top, blue at the bottom
    colormap(colors);
    imagesc(sub_pca_dat.sess_sentence_elec_stim_mat)
    set(gca, 'ydir', 'reverse','box','on');
    xlabel(condition)
    ylabel('Voxels (Electrodes)')
    title([condition,' Responses'])
    ax.XTick=[];
    ax.YTick=[];
    hold on 
    ax.FontSize=12;
    %arrayfun(@(x) plot([x,x],ax.YLim,'w-'),[1:8:size(pca_mat,1)])
    %
    ax=axes('position',[.1,.1,.3,.2])
    colors = cbrewer('div', 'RdBu', 128);
    colors = flipud(colors); % puts red on top, blue at the bottom
    colormap(colors);
    imagesc(W_inferred)
    set(gca, 'ydir', 'reverse','box','on');
    title('Weights')
    ax.XLabel.String='Voxels (Electrodes)';
    ax.YLabel.String='Component';
    ax.XTick=[];
    ax.YTick=[];
    ax.FontSize=12;
    %
    ax=axes('position',[.45,.1,.2./aspect_ration,.3*aspect_ration])
    colors = cbrewer('div', 'RdBu', 128);
    colors = flipud(colors); % puts red on top, blue at the bottom
    colormap(colors);
    imagesc(R_inferred)
    set(gca, 'ydir', 'reverse','box','on');
    title('Features');
    ax.YLabel.String='sentence';
    ax.XLabel.String='Component';
    ax.XTick=[];
    ax.YTick=[];
    ax.FontSize=12;
    
    % 
    ax=axes('position',[.6,.1,(.2./aspect_ration)/5,.3*aspect_ration])
    colors = cbrewer('div', 'RdBu', 128);
    colors = flipud(colors); % puts red on top, blue at the bottom
    colormap(colors);
    imagesc(R_inferred(:,1),[min(R_inferred(:)),max(R_inferred(:))])
    
    set(gca, 'ydir', 'reverse','box','on');
    title('Feature 1');
    ax.XTick=[];
    ax.YTick=[];
    ax.FontSize=12;
    
    % 
    
     ax=axes('position',[.67,.1,(.2./aspect_ration)/5,.3*aspect_ration])
    colors = cbrewer('div', 'RdBu', 128);
    colors = flipud(colors); % puts red on top, blue at the bottom
    colormap(colors);
   
    
    F1 = min(R_inferred(:)) + (max(R_inferred(:))-min(R_inferred(:)))*rand(size(R_inferred,1),1);
    imagesc(F1,[min(R_inferred(:)),max(R_inferred(:))]);
    title('F_1');
    ax.FontSize=12;
    set(gca, 'ydir', 'reverse','box','on');
    ax.XTick=[];
    ax.YTick=[];
    
    % 
    analysis_path='/Users/eghbalhosseiniasl1/MyData/ECoG-sentence/analysis/swjn_Aim_2_graphics';
    if ~exist(strcat(analysis_path))
            mkdir(strcat(analysis_path))
    end 
    print(fig1, '-djpeg', strcat(analysis_path,'/ICA_analysis.jpeg'));
end

