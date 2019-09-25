function tree_branch=make_hierarchical_tree(data)
% data is a list of opening and closings 

data_idx=1:length(data);
data_in=data;
tree_branch={};
tree_plot_pairs={};
i=1;
f_k_means=[];
while ~isempty(data_idx)
   f=[data_idx(find(diff(data_in)==-2));data_idx(find(diff(data_in)==-2)+1)]; 
   tree_branch{i}=f;
   data_temp=data_in;
   data_idx([(find(diff(data_temp)==-2));(find(diff(data_temp)==-2)+1)])=[];
   data_temp([(find(diff(data_temp)==-2));(find(diff(data_temp)==-2)+1)])=[];
   data_in=data_temp;   
   i=i+1;
end 
% make points for creating the tree
%figure
%hold on 
%cellfun(@(x) plot(x(1,:),x(2,:),'k-'), tree_plot_pairs,'UniformOutput',false)
end 


%for k=1:size(f,2)
%        if i==1
%         f_k_means_x=mean(f,1);
%         f_k_means_y=i+0*mean(f,1);
%         tree_plot_x=[f(1,k),mean(f(:,k)),f(2,k)];
%         tree_plot_y=[i-1,i,i-1];
%         tree_plot_pairs=[tree_plot_pairs,[tree_plot_x;tree_plot_y]];
%        else 
%           f_k_mean_in_range=f_k_means_x>f(1,k) & f_k_means_x<f(2,k);
%           f_k_mean_in_range_x=f_k_means_x(f_k_mean_in_range);
%           f_k_mean_in_range_y=f_k_means_y(f_k_mean_in_range);
%           if ~isempty(f_k_mean_in_range_x)
%           tree_plot_x=[f_k_mean_in_range_x(1),mean(f_k_mean_in_range_x),f_k_mean_in_range_x(2)];
%           tree_plot_y=[f_k_mean_in_range_y(1),sum(f_k_mean_in_range_y),f_k_mean_in_range_y(2)];
%           tree_plot_pairs=[tree_plot_pairs,[tree_plot_x;tree_plot_y]];
%           % remove previous point 
%           f_k_means_x(f_k_mean_in_range)=[];
%           f_k_means_y(f_k_mean_in_range)=[];
%           % add new points 
%           f_k_means_x=[f_k_means_x,mean(f_k_mean_in_range_x)];
%           f_k_means_y=[f_k_means_y,sum(f_k_mean_in_range_y)];
%           end 
%        end      
%end 

% for k=1:size(f,2)
%        if i==1
%         f_k_means_x=mean(f,1);
%         tree_plot_x=[f(1,k),mean(f(:,k)),f(2,k)];
%         tree_plot_pairs=[tree_plot_pairs,[tree_plot_x]];
%        else 
%           f_k_mean_in_range=f_k_means_x>f(1,k) & f_k_means_x<f(2,k);
%           f_k_mean_in_range_x=f_k_means_x(f_k_mean_in_range);
%           while ~isempty(f_k_mean_in_range_x)
%           tree_plot_x=[mean([f_k_mean_in_range_x(1),f_k_mean_in_range_x(2)])];
%           tree_plot_pairs=[tree_plot_pairs,[tree_plot_x]];  
%           f_k_means_x=([f_k_means_x,mean(f_k_mean_in_range_x)]);
%           % remove previous point 
%           f_k_means_x(f_k_mean_in_range)=[];
%           % add new points 
%           
%           end 
%  