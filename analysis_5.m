
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
%% 
find_index= @(input_cell) find(cell2mat(cellfun(@(x) ~isempty(x), input_cell,'UniformOutput',false)));
%electrode_pattern_section_tensor=nan*ones(sum(language_electrode),size(unique_open_node_patterns,1),8*135);
electrode_pattern_section_cell= cell(size(unique_open_node_patterns,1),8);
all_example_sentence_locations=[];

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
    sentences=[data{sentence_trial_index}];
    sentence_gamma_band_envelope=cellfun(@(x) x(1:8),{sentences.signal_hilbert_zs_downsample_parsed},'UniformOutput',false);
    % creat a cell with wordposition(row)*time in trial(column) structure
    sentence_gamma_band_envelope=[sentence_gamma_band_envelope{:,:}];
    sentence_lang_elec_gamma_band_envelope=cellfun(@(x) x(language_electrode_num,:),sentence_gamma_band_envelope,'UniformOutput',false);
    
    % find sentence locations 
    example_sentence={sentences(:).trial_string};
    example_sentence=cellfun(@(x) x(2:end),example_sentence,'UniformOutput',false);
    % find trial location in node cell 
    example_sentence_locations=cellfun(@(x) (regexpi([node_cell{:,1}],x)),example_sentence,'UniformOutput',false);
    example_sentence_locations=cell2mat(cellfun(@(x) find_index(x), example_sentence_locations, 'UniformOutput',false));
    all_example_sentence_locations=[all_example_sentence_locations,example_sentence_locations];
    
    for k=1:length(example_sentence_locations)
        
        sentence_elec_gamma_bands=sentence_lang_elec_gamma_band_envelope(:,k);
        sentence_pattern_membership=all_pattern_locations(example_sentence_locations(k),2:end);
        total_iterations=sum(cell2mat(cellfun(@length,sentence_pattern_membership,'UniformOutput',false)));
        fprintf('doing sentence %d total iterations %d \n',example_sentence_locations(k),total_iterations);
        for kk=1:size(sentence_elec_gamma_bands,1)
            section_membership=sentence_pattern_membership{kk};
            fprintf('doing section %d \n',kk);
            for kkk=1:length(section_membership)
                current_cell_state=electrode_pattern_section_cell{section_membership(kkk),kk};
                C = cat(3,current_cell_state,sentence_elec_gamma_bands{kk});
                
                electrode_pattern_section_cell{section_membership(kkk),kk}=C;
                % assign the channel response to the s
            end 
        end
        
    end    
    % convert the gamma into a tensor of shape : channel*time in trial*
    % trial number
    %sentence_lang_elec_gamma_band_envelope_tensor=cell2mat(permute(sentence_lang_elec_gamma_band_envelope,[3,1,2]));

    fprintf('added %s from %s \n',d(i).name, strcat(d(i).folder,'/',d(i).name));
    
end


%% plot example electrode for all conditions close all;
pattern_num=length(unique_open_node_patterns);

colors = cbrewer('div', 'RdBu', 128);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);

num_rows=4;
num_columns=2;
total_plots=num_rows*num_columns;

for i=1:pattern_num
    pattern_response=electrode_pattern_section_cell(i,:);
    pattern_mean_ch_response=cellfun(@(x) mean(x,3), pattern_response,'UniformOutput',false);
    template_mat=pattern_mean_ch_response{1,1}*nan;
    % find empty secions
    empty_cell=find(cell2mat(cellfun(@isempty,pattern_mean_ch_response,'UniformOutput',false)));
    for nn=empty_cell
        pattern_mean_ch_response{1,nn}=template_mat;
    end
    pattern_mean_ch_response=cell2mat(pattern_mean_ch_response);
    
    figure(fix((i-1)/total_plots)+1);
    colors = cbrewer('div', 'RdBu', 128);
    colors = flipud(colors); % puts red on top, blue at the bottom
    colormap(colors);
    %set(gcf,'Position',[-882 449 877 1325]);
    subplot(num_rows,num_columns,i-total_plots*fix((i-1)/total_plots));
    imagesc(squeeze(pattern_mean_ch_response));
    ylim=get(gca,'ylim');
    text(0.0:135:size(pattern_mean_ch_response,2)-1,ylim(2)+zeros(1,8),strsplit(num2str(all_pattern_locations{i,1}),' '),'FontSize',12,...
        'VerticalAlignment','bottom');
    set(gca, 'ydir', 'normal','box','off','xtick',0.0:135:size(pattern_mean_ch_response,2)-1,'xticklabel',[0:8]);
    handles = colorbar;
    xlabel('word position'); ylabel('Ch');    
    

end  
%% find patterns with highest correlations 
pattern_corr_mat=corr(unique_open_node_patterns');
[sort_pattern,ind]=sort(pattern_corr_mat,2);
sort_pattern=fliplr(sort_pattern);
ind=fliplr(ind);
num_of_pattern_to_plot_together=3;
%% plot each electrode response for different pattern 
% add nan to the electrodes 
empty_cell=find(cell2mat(cellfun(@isempty,electrode_pattern_section_cell,'UniformOutput',false)));
template_mat=electrode_pattern_section_cell{1,1,1}*nan;
for nn=empty_cell'
    electrode_pattern_section_cell{nn}=template_mat;
end
%% 
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
%% 
num_rows=5;
num_columns=7;
total_plots=num_rows*num_columns;

for i=1:length(language_electrode_num)
    electrode_response=cellfun(@(x) x(i,:,:),electrode_pattern_section_cell_modified,'UniformOutput',false);
    electrode_mean_response=cellfun(@(x) nanmean(x,3), electrode_response,'UniformOutput',false);
    % get max size of electrode cell 
    max_size=max(max(cell2mat(cellfun(@(x) size(x,3),electrode_response,'UniformOutput',false))));
    % add nan to fill the matrix 
    % find non_max cells 
    non_max_cell=find(cell2mat(cellfun(@(x) size(x,3)~=max_size, electrode_response,'UniformOutput',false))) ;
    for nn=non_max_cell'
        cell_response=[];
        cell_response=electrode_response{nn};
        cell_dummy_response=nan*ones(size(cell_response,1),size(cell_response,2),max_size-size(cell_response,3));
        electrode_response{nn}=cat(3,cell_response,cell_dummy_response);
    end
    electrode_response=cellfun(@(x) squeeze(x),electrode_response,'UniformOutput',false);
    colors = cbrewer('qual', 'Set2', size(unique_open_node_patterns,1));
    %figure(fix((i-1)/total_plots)+1);
    %close all;
    figure('position',[-2154 536 2111 1207]);
    hold on;
    y_lims=[];
    for ii=1:size(unique_open_node_patterns,1)
        ax{ii}=subplot(num_rows,num_columns,ii);
        electrode_unique_pattern_response=double(transpose(cell2mat(transpose(electrode_response(ii,:)))));
        time=1:size(electrode_unique_pattern_response,2);
        mean_electrode_pattern_response=nanmean(electrode_unique_pattern_response,1);
        std_electrode_pattern_response=nanstd(electrode_unique_pattern_response,1)./sqrt( ...
            sum(~isnan(electrode_unique_pattern_response ),1));
        non_nan_times=~isnan(mean_electrode_pattern_response);
        bl = boundedline(time(non_nan_times), mean_electrode_pattern_response(non_nan_times), std_electrode_pattern_response(non_nan_times), ...
            'cmap', [1,0,0],'alpha');
        % make the sudo line for open node pregression
        pattern_prgression=reshape(time,[],length(all_pattern_locations{ii,1}));
        diff_pattern=diff(all_pattern_locations{ii,1});
        
        pattern_progression=repmat(diff_pattern,size(pattern_prgression,1),1);
        pattern_pregression=cumsum(pattern_progression(:))./max(cumsum(pattern_progression(:)));
        bl = boundedline(time(1:length(pattern_pregression)), pattern_pregression,nan*pattern_pregression,'--', ...
            'cmap', [.5,.5,.5],'alpha');
        % plot 2 other highly correlated pattern
        %     electrode_corr1_pattern_response=double(transpose(cell2mat(transpose(electrode_response(ind(ii,2),:)))));
        %     mean_electrode_corr1_pattern_response=nanmean(electrode_corr1_pattern_response,1);
        %         std_electrode_corr1_pattern_response=nanstd(electrode_corr1_pattern_response,1)./sqrt( ...
        %             sum(~isnan(electrode_corr1_pattern_response ),1));
        %         non_nan_times=~isnan(mean_electrode_corr1_pattern_response);
        %         bl = boundedline(time(non_nan_times), mean_electrode_corr1_pattern_response(non_nan_times), std_electrode_corr1_pattern_response(non_nan_times), ...
        %         'cmap', [.7,.7,1],'alpha');
        %
        ylim=get(gca,'ylim');
        y_lims=[y_lims;get(gca,'ylim')];
        set(gca,'xlim',[0,max(time)]);
        %text(0.0:135:max(time)-1,ylim(2)+zeros(1,8),strsplit(num2str(all_pattern_locations{ii,1}),' '),'FontSize',8,...
        %    'VerticalAlignment','top');
        set(gca, 'ydir', 'normal','box','off','xtick',0.0:135:max(time)-1,'xticklabel',strsplit(num2str(all_pattern_locations{ii,1}),' '));
        if ii==size(unique_open_node_patterns,1)
           
             set(gca, 'ydir', 'normal','box','off','xtick',0.0:135:max(time)-1,'xticklabel',strsplit(num2str(all_pattern_locations{ii,1}),' '));
             xlabel('open node pattern');ylabel('z-score');
             text(max(time),ylim(1),strcat(' ',info.subject,'- electrode ',num2str(language_electrode_num(i))),'FontSize',12,...
            'HorizontalAlignment','Left');
             
              
        end
        
    end
    new_bounds=quantile(y_lims(:),3);
       for l=1:size(unique_open_node_patterns,1)
           ax{l}.YLim=[new_bounds(1),new_bounds(end)];
       end 
    %pause(5)
    %print(gcf, '-depsc', strcat('/Users/eghbalhosseiniasl1/MyData/ecog-sentence/syntax-response/',...
    %        info.subject,'_channel_',num2str(language_electrode_num(i)),'_response_syntax_comparison','.eps'));
end 

%% 

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

%% sort trials based on value of last merge 
id_1=10;
id_2=30;

P1=squeeze(double(nanmean(cat(2,electrode_pattern_section_cell_modified{id_1,:}),1)));
P2=squeeze(double(nanmean(cat(2,electrode_pattern_section_cell_modified{id_2,:}),1)));

time=1:size(P1,1);
mean_electrode_pattern_response=nanmean(P1,2);
std_electrode_pattern_response=nanstd(P1')./sqrt(sum(~isnan(P1),2));
non_nan_times=~isnan(mean_electrode_pattern_response);

figure;hold on; 
bl = boundedline(time(non_nan_times), mean_electrode_pattern_response(non_nan_times), std_electrode_pattern_response(non_nan_times),'cmap', [1,0,0],'alpha');
ylim=get(gca,'ylim');set(gca,'xlim',[0,max(time)]);
text(0.0:135:max(time)-1,nanmean(reshape(mean_electrode_pattern_response,135,[]),1)+.05,strsplit(num2str(unique_open_node_patterns(id_1,:)),' '),...
    'FontSize',12,'color',[.6,0,0],'VerticalAlignment','bottom','FontWeight','bold');

% 
time=1:size(P2,1);
mean_electrode_pattern_response=nanmean(P2,2);
std_electrode_pattern_response=nanstd(P2')./sqrt(sum(~isnan(P2),2));
non_nan_times=~isnan(mean_electrode_pattern_response);

bl = boundedline(time(non_nan_times), mean_electrode_pattern_response(non_nan_times), std_electrode_pattern_response(non_nan_times),'cmap', [0,0,1],'alpha');
text(0.0:135:max(time)-1,nanmean(reshape(mean_electrode_pattern_response,135,[]),1)+.05,strsplit(num2str(unique_open_node_patterns(id_2,:)),' '),...
    'FontSize',12,'color',[0,0,.6],'VerticalAlignment','top','FontWeight','bold');


set(gca, 'ydir', 'normal','box','off','xtick',0.0:135:max(time)-1,'xticklabel',[1:8]);
xlabel('word position');ylabel('z-score');
