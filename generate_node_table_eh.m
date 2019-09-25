%% Generate Node Table %%
% by eghbal hosseini
%Creates a struct with open and closed node information

%% load tree %%
function [output,output_table]=generate_node_table_eh(data_path)
sentences_ptb = readtable(data_path);

sentence_trial_index=~cellfun(@isempty,cellfun(@ (x) strfind(x,'.'),sentences_ptb.word,'UniformOutput',false));
%% get matrices of nodes for all sentences%%
output={};
sentence_trial_bounds=[0;find(sentence_trial_index)];
for sentence_idx = 1:length(sentence_trial_bounds)-1 % get node information sentence by sentence
    % get sentence
    sentence_range=[sentence_trial_bounds(sentence_idx)+1:sentence_trial_bounds(sentence_idx+1)-1];
    sentence=cellfun(@(x) [x,' '],sentences_ptb.word(sentence_range),'UniformOutput',false);
    sentence=[sentence{:}];
    try 
        sentence = strrep(sentence," '","'");
        
    end 
    try 
        sentence = strrep(sentence,"can not","cannot");
        
    end 
    % get node closing 
    node_closing=sentences_ptb.nmerge(sentence_range);
    node_open=[];
    for word = 1:length(node_closing) %get node information for each word
        if (word == 1 || word == 2) 
        node_open(word) = word;
        else
        node_open(word) = node_open(word-1) - node_closing(word-1) + 1; %open node formula from Nelson et al. 2017
        end
    end
    output=[output;{sentence,node_open,node_closing'}];
end
output_table = cell2table(output,...
    'VariableNames',{'sentence' 'open_nodes' 'closed_nodes'});


end    
  
