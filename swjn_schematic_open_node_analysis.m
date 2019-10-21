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
analysis_path='/Users/eghbalhosseiniasl1/MyData/ECoG-sentence/analysis/swjn_schematic_open_node_analysis/';
%
if 1
    fprintf('adding basic ecog tools to path \n');
    addpath('~/MyCodes/basic-ecog-tools/');
    addpath('~/MyCodes/ecog-sentence/');
    addpath(genpath('~/MyCodes/basic-ecog-tools/activeBrain'));
    addpath(genpath('~/MyCodes/basic-ecog-tools/ecog-filters'));
    addpath(genpath('~/MyCodes/basic-ecog-tools/mex'));
end 
print_opt=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find the node information 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d_nmerge= dir([data_path,'/**/*nmerge.txt']);
fprintf(' %d nmerge files were found \n', length(d_nmerge));
i=2;
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
node_table(1,:)
sentence=node_cell{1,1};
open_nodes=node_cell{1,2};
closed_nodes=node_table.closed_nodes{1};

extra_branch=sum(closed_nodes);
node_x=1:8;
node_y=ones(1,8);
n=8;
tape=[];

x=(1:8);
y=ones(1,8);
node_con=[];
for i=1:length(closed_nodes)
    cur=closed_nodes(i)
    tape=[tape,i]
    
    while cur~=0
        cur=cur-1;
        new_node=n+1;
        
        merge=tape(end-1:end);
        new_node_x=mean(x(merge));
        new_node_y=max(y(merge))+1;
        x=[x,new_node_x];
        y=[y,new_node_y]
        for k=1:length(merge)
        node_con=[node_con;[merge(k),new_node] ]   
        end 
        tape=tape(1:end-1);
        tape(end)=new_node;
        n=n+1;
    end 
end 

%node_mat=zeros(1,n)
%node_mat(node_con(:,1))=(n+1)-node_con(:,2)
%tree_vec=fliplr(node_mat);
%treeplot(tree_vec)
%a=get(gca)
%a.XDir='reverse'
%shg
% plotting 
f=figure;
ax=axes('Position',[.2,.3,.6,.6]);
hold on 
bl=arrayfun(@(p) plot(x(p{:}),y(p{:}),'k'), num2cell(node_con', 1),'UniformOutput',false);
bl1=plot(x,y,'k.')
bl1.MarkerSize=20

ax.XTick=[1:8];
sentence_split=strsplit(sentence);
ax.XTickLabel=sentence_split(1:8);
ax.XAxis.Visible='off'
ax.YAxis.Visible='off'
daspect([1,2,1]);

arrayfun(@(p) text(x(p),y(p)-.8,sentence_split{p},'verticalalignment','bottom','horizontalalignment','center','fontsize',10),1:8)
arrayfun(@(p) text(x(p),y(p)-1.2,sprintf('%1.0d',open_nodes(p)),'verticalalignment','bottom','horizontalalignment','center','fontsize',8),1:8)
arrayfun(@(p) text(x(p),y(p)-1.6,sprintf('%1.1d',closed_nodes(p)),'verticalalignment','bottom','horizontalalignment','center','fontsize',8),1:8)
arrayfun(@(p) text(x(p),y(p)-2,sprintf('%1.1d',(p)),'verticalalignment','bottom','horizontalalignment','center','fontsize',8),1:8)
text(x(1)-1.2,y(1)-1.2,'open nodes:','verticalalignment','bottom','horizontalalignment','left','fontsize',8);
text(x(1)-1.2,y(1)-1.6,'node closing:','verticalalignment','bottom','horizontalalignment','left','fontsize',8);
text(x(1)-1.2,y(1)-2,'word position:','verticalalignment','bottom','horizontalalignment','left','fontsize',8);



   if ~exist(strcat(analysis_path))
                mkdir(strcat(analysis_path));
            end
            if print_opt==1
                set(gcf,'PaperPosition',[.25 .25 8 6])
                print(f,'-painters', '-djpeg', strcat(analysis_path,'/','swjn_schematic_open_node.jpeg'));
            end