clear all 
close all 
analysis_path='/Users/eghbalhosseiniasl1/MyData/ECoG-sentence/analysis/swjn_Aim_1_graphics';
data_path='/Users/eghbalhosseiniasl1/MyData/ecog-sentence';
if 1
    fprintf('adding basic ecog tools to path \n');
    addpath('~/MyCodes/basic-ecog-tools/');
    addpath(genpath('~/MyCodes/ecog-sentence/'));
    addpath(genpath('~/MyCodes/basic-ecog-tools/activeBrain'));
    addpath(genpath('~/MyCodes/basic-ecog-tools/ecog-filters'));
    addpath(genpath('~/MyCodes/basic-ecog-tools/mex'));
end
d_pmi= dir([data_path,'/**/*pmi_sentences.csv']);
fprintf(' %d pmi files were found \n', length(d_pmi));
i=1;
[pmi_cell]=generate_pmi_table_eh(strcat(d_pmi(i).folder,'/',d_pmi(i).name));
subject_id='AMC026';
d= dir([data_path,strcat('/**/',subject_id,'*_crunched.mat')]);
fprintf(' %d .mat files were found \n', length(d));

d_nmerge= dir([data_path,'/**/*nmerge.txt']);
fprintf(' %d nmerge files were found \n', length(d_nmerge));
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

%%  
AMC026_elec_fsaverage=load(['/Users/eghbalhosseiniasl1/MyData/ecog-sentence/subjects_raw/AMC026/IMAGING/','AMC026_elec_fsavg.mat']);
AMC026_elec_fsaverage=AMC026_elec_fsaverage.elec_fsavg_frs;
AMC029_elec_fsaverage=load(['/Users/eghbalhosseiniasl1/MyData/ecog-sentence/subjects_raw/AMC029/IMAGING/','AMC029_elec_fsavg.mat']);
AMC029_elec_fsaverage=AMC029_elec_fsaverage.elec_fsavg_frs;

%% get 

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
    language_electrode=info.ramp_electrodes_hilbert_odd;
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
    words_trial_index=~cellfun(@isempty,cellfun(@ (x) strfind(x,'N'),info.word_type,'UniformOutput',false));
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
%% 
idx=17;
electrode_id=language_electrode_num(idx);

sent_ramp_reponse=squeeze(session_sentence_hilbert_band_envelope_lang_elec_tensor(idx,:,:));
words_ramp_reponse=squeeze(session_words_hilbert_band_envelope_lang_elec_tensor(idx,:,:));

%% 
colors=cbrewer('qual','Paired',10);
close all 
f=figure;
set(f,'position',[-1133 382 847 847]);
ax=axes('position',[.1,.4,.5,.5]);
fspial_lh = ft_read_headshape('/Applications/freesurfer/subjects/fsaverage/surf/lh.pial');
fspial_lh.coordsys = 'fsaverage';
h=ft_plot_mesh(fspial_lh);
h.DisplayName='Average';
h_1=ft_plot_sens(AMC026_elec_fsaverage,'elecshape','point','elecsize',35,'facecolor',colors(4,:),'edgecolor',[0,0,0]);
h_1.DisplayName='AMC026';
h_2=ft_plot_sens(AMC029_elec_fsaverage,'elecshape','point','elecsize',35,'facecolor',colors(10,:),'edgecolor',[0,0,0]);
h_2.DisplayName='AMC029';
view([-100 10]);
material dull;
lighting gouraud;
camlight; 
legend('position',[.16,.45,.04,0.05],'fontsize',18)

a=annotation(f,'textbox',[.1 .74 .1 .1],'String','A','FontSize',20,'FontWeight','bold');
a.LineStyle='none';

colors=cbrewer('qual','Set1',10);
ax=axes('position',[.20,.32,.3,.10]);
y1=double(mean(sent_ramp_reponse,1));
x=1:length(y1);
e1=double(nanstd(sent_ramp_reponse,1)./sqrt(sum(~isnan(sent_ramp_reponse),1)));
[l,p] = boundedline(x, y1, e1, '-b');
hAnnotation=arrayfun(@(x) get(x,'Annotation'),p);
hLegendEntry = arrayfun(@(x) get(x,'LegendInformation'),hAnnotation);
arrayfun(@(x) set(x,'IconDisplayStyle','off'),hLegendEntry);
l.DisplayName='S';
l.LineWidth=2;
hold on 
y1=double(mean(words_ramp_reponse,1));
x=1:length(y1);
e1=double(nanstd(words_ramp_reponse,1)./sqrt(sum(~isnan(words_ramp_reponse),1)));
[l,p] = boundedline(x, y1, e1, '-r');
l.LineWidth=2;
hAnnotation=arrayfun(@(x) get(x,'Annotation'),p);
hLegendEntry = arrayfun(@(x) get(x,'LegendInformation'),hAnnotation);
arrayfun(@(x) set(x,'IconDisplayStyle','off'),hLegendEntry);
l.DisplayName='N';
leg=legend('position',[.53,.4,.03,0.05],'fontsize',18);
ax.YTick=[];
ax.XLim=[-150,1080];
ax.XTick=[-150,1:135:1080];
ax.XTickLabel={'Word:','1','2','3','4','5','6','7','8'};
ax.XTickLabelMode
ax.YAxis.Label.String='\gamma'
ax.FontSize=18
ax.XLim=[-150,1080];
a=annotation(f,'textbox',[.1 .32 .1 .1],'String','B','FontSize',20,'FontWeight','bold');
a.LineStyle='none';
string=['\bf','Figure 1. ', ' \rm\color{black}','(A) Example of two subjects overlayed on an average brain,','(B)','\rm ','Response of a sample electrode to S and N conditions']
%a=annotation(f,'textbox',[.1 .17 .55 .1],'String',...
%    'Figure.1 (A) Example of two subjects overlayed on an average brain, (B) Response of a sample electrode to S and N conditions ',...
%    'FontSize',25);

a=annotation(f,'textbox',[.1 .17 .55 .1],'String',...
    string,...
    'FontSize',25);
a.LineStyle='none';
 
if ~exist(strcat(analysis_path))
    mkdir(strcat(analysis_path))
end

print(f, '-djpeg', strcat(analysis_path,'/average_anatomicals.jpeg'));
%% for presentation 

colors=cbrewer('qual','Set1',10);
close all 
f=figure;
f.Units = 'inches';
set(f,'position',[-1133 382 847 847]);

set(f, 'Position', [9.9167 , 7.0, 9.32,4.13]);
ax=axes('position',[.1,.4,.5,.5]);
fspial_lh = ft_read_headshape('/Applications/freesurfer/subjects/fsaverage/surf/lh.pial');
fspial_lh.coordsys = 'fsaverage';
h=ft_plot_mesh(fspial_lh);
h.DisplayName='Average';
h_1=ft_plot_sens(AMC026_elec_fsaverage,'elecshape','point','elecsize',5,'facecolor',colors(1,:),'edgecolor',[0,0,0]);
h_1.DisplayName='AMC026';
h_2=ft_plot_sens(AMC029_elec_fsaverage,'elecshape','point','elecsize',5,'facecolor',colors(2,:),'edgecolor',[0,0,0]);
h_2.DisplayName='AMC029';
view([-100 10]);
material dull;
lighting gouraud;
camlight; 
legend('position',[.16,.45,.04,0.05],'fontsize',18)

a=annotation(f,'textbox',[.1 .74 .1 .1],'String','A','FontSize',20,'FontWeight','bold');
a.LineStyle='none';
%
k1=70;y1=double(sent_ramp_reponse(k1,:));
k2=45;y2=double(sent_ramp_reponse(k2,:));
% 
k3=70;y3=double(words_ramp_reponse(k3,:));
k4=72;y4=double(words_ramp_reponse(k4,:));
% 
ax_max=max([y1,y2,y3,y4]);
ax_min=min([y1,y2,y3,y4]);

%
ax=axes('position',[.65,.92,.3,0.05]);
x=1:length(y1);
l=plot(x,y1);
hold on
l.LineWidth=2;
ax.XLim=[-100,1080];
ax.YLim=[ax_min,1.2*ax_max];
ax.XTick=[1:135:1080];
arrayfun(@(x)plot([x,x],ax.YLim,'k--'),ax.XTick);
ticks=strsplit(session_sentence_examples{k1},' ');
str=arrayfun(@(x) text(ax.XTick(x),ax.YLim(2),ticks(x),...
    'VerticalAlignment','top','HorizontalAlignment','left','FontSize',6,'FontWeight','bold','Rotation',45),1:8);
ax.YAxis.Visible='off';
ax.XAxis.Visible='off';

ax=axes('position',[.65,.82,.3,0.05]);
l=plot(x,y2);
hold on
l.LineWidth=2;
ax.XLim=[-100,1080];
ax.YLim=[ax_min,1.2*ax_max];
ax.XTick=[1:135:1080];
arrayfun(@(x)plot([x,x],ax.YLim,'k--'),ax.XTick);

ticks=strsplit(session_sentence_examples{k2},' ');
str=arrayfun(@(x) text(ax.XTick(x),ax.YLim(2),ticks(x),...
    'VerticalAlignment','top','HorizontalAlignment','left','FontSize',6,'FontWeight','bold','Rotation',45),1:8)
ax.YAxis.Visible='off';
ax.XAxis.Visible='off';
box off

% words 

ax=axes('position',[.65,.65,.3,0.05]);
l=plot(x,y3,'r');
hold on
l.LineWidth=2;
ax.XLim=[-100,1080];
ax.YLim=[ax_min,1.2*ax_max];
ax.XTick=[1:135:1080];
arrayfun(@(x)plot([x,x],ax.YLim,'k--'),ax.XTick)
ticks=strsplit(session_wordlist_examples{k3},' ')

str=arrayfun(@(x) text(ax.XTick(x),ax.YLim(2),ticks(x),'VerticalAlignment','top','HorizontalAlignment','left','FontSize',10,'FontWeight','bold','Rotation',45),1:8)
ax.YAxis.Visible='off';
ax.XAxis.Visible='off';
box off

ax=axes('position',[.65,.55,.3,0.05]);
l=plot(x,y4,'r');
hold on
l.LineWidth=2;
ax.XLim=[-100,1080];
ax.YLim=[ax_min,1.2*ax_max];
ax.XTick=[1:135:1080];
arrayfun(@(x)plot([x,x],ax.YLim,'k--'),ax.XTick)

ticks=strsplit(session_wordlist_examples{k4},' ')
str=arrayfun(@(x) text(ax.XTick(x),ax.YLim(2),ticks(x),'VerticalAlignment','top','HorizontalAlignment','left','FontSize',10,'FontWeight','bold','Rotation',45),1:8)
ax.YAxis.Visible='on';
ax.XAxis.Visible='off';
box off

% 
%
hold on 
y1=double(sent_ramp_reponse(75,:));
x=1:length(y1);
a=plot3(x,x*0+1,y1);
% 
ax=axes('position',[.20,.32,.3,.10]);
y1=double(mean(sent_ramp_reponse,1));
x=1:length(y1);
e1=double(nanstd(sent_ramp_reponse,1)./sqrt(sum(~isnan(sent_ramp_reponse),1)));
[l,p] = boundedline(x, y1, e1, '-b');
hAnnotation=arrayfun(@(x) get(x,'Annotation'),p);
hLegendEntry = arrayfun(@(x) get(x,'LegendInformation'),hAnnotation);
arrayfun(@(x) set(x,'IconDisplayStyle','off'),hLegendEntry);
l.DisplayName='S';
l.LineWidth=2;
hold on 
y1=double(mean(words_ramp_reponse,1));
x=1:length(y1);
e1=double(nanstd(words_ramp_reponse,1)./sqrt(sum(~isnan(words_ramp_reponse),1)));
[l,p] = boundedline(x, y1, e1, '-r');
l.LineWidth=2;
hAnnotation=arrayfun(@(x) get(x,'Annotation'),p);
hLegendEntry = arrayfun(@(x) get(x,'LegendInformation'),hAnnotation);
arrayfun(@(x) set(x,'IconDisplayStyle','off'),hLegendEntry);
l.DisplayName='N';
leg=legend('position',[.53,.4,.03,0.05],'fontsize',18);
ax.YTick=[];
ax.XLim=[-150,1080];
ax.XTick=[-150,1:135:1080];
ax.XTickLabel={'Word:','1','2','3','4','5','6','7','8'};
ax.XTickLabelMode
ax.YAxis.Label.String='\gamma'
ax.FontSize=18
ax.XLim=[-150,1080];
a=annotation(f,'textbox',[.1 .32 .1 .1],'String','B','FontSize',20,'FontWeight','bold');
a.LineStyle='none';
a=annotation(f,'textbox',[.1 .17 .55 .1],'String',...
    'Figure.1 (A) Example of two subjects overlayed on an average brain, (B) Response of a sample electrode to S and N conditions ',...
    'FontSize',25);
a.LineStyle='none';
 
if ~exist(strcat(analysis_path))
    mkdir(strcat(analysis_path))
end

print(f, '-djpeg', strcat(analysis_path,'/average_anatomicals.jpeg'));