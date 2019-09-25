cd ~/MyData/ecog-sentence/exp_stimuli/
d_stimuli= dir([cd,'/*words.XLS']);
fprintf(' %d stimuli files were found \n', length(d_stimuli));
word_list = readtable(d_stimuli(1).name);
trial_id=table2array(word_list(:,4));
unique_trial_id=unique(trial_id);
word_location_in_trial=transpose(reshape(1:size(trial_id,1),[],size(unique_trial_id,1)));
table2cell(word_list(word_location_in_trial(1,:),2))
all_words=[];

fileID = fopen('words.txt','w');
for i=1:length(unique_trial_id)
    trial_words=table2cell(word_list(word_location_in_trial(i,:),2));
    trial_words=join(trial_words);
    all_words=[all_words;trial_words];
    fprintf(fileID,'%s\n',trial_words{:});
end 
fclose(fileID);
