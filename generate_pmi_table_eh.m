%% Generate pmi Table %%
% by eghbal hosseini
%Creates a struct with pmi information

%% load tree %%
function [output]=generate_pmi_table_eh(data_path)
sentences_ptb = readtable(data_path);
sentences=table2cell(unique(cell2table(sentences_ptb.Var1)));
sentences_pmi=[];
% 
for i=1:length(sentences)
    sentence=sentences(i);
    sentence_index=find(ismember(sentences_ptb.Var1,sentence));
    sentence_table=sentences_ptb(sentence_index,:);
    sentence_first_word=sentence_table.Var4;
    sentence_second_word=sentence_table.Var5;
    sentence_pmi=sentence_table.Var6;
    pmi_mat=zeros(max(sentence_second_word)+1);
    for k=1:length(sentence_first_word)
    pmi_mat(sentence_first_word(k)+1,sentence_second_word(k)+1)=sentence_pmi(k);
    end 
    sentences_pmi=[sentences_pmi;sentence,pmi_mat];
end 
output=sentences_pmi;
end    
  
