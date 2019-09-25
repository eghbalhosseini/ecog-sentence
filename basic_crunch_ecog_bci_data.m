function basic_crunch_ecog_bci_data(data_path)
    
    fprintf('extracting from %s \n',data_path);
    [signal_broadband,signal_bandpass,signal_envelope,states,parameters]=filter_channels_using_schalk({data_path});
    subject_name=data_path(strfind(data_path,'AMC')+[0:5]);
    session_name=data_path(strfind(data_path,'ECOG'):length(data_path)-4);
    %[ signal, states, parameters ] = load_bcidat( strcat(d(i).folder,'/',d(i).name));
    % calculate the signal at different frequency scales: 
    
    % start with an empty strcuture for data and info 
    list_var_to_get={''};
    stim_types={'S','W','N','J'};
    dat={};
    info=struct;
    
    % step 1: find start and end of trials 
    info.sample_rate=parameters.SamplingRate.NumericValue
    stimuli_squence=parameters.Sequence.NumericValue;
    trials_value=parameters.Stimuli.NumericValue;
    stimuli_value=parameters.Stimuli.Value;
    %
    trials_indx=cell2mat(cellfun(@(x) contains(x,'TrialNumber'),parameters.Stimuli.RowLabels,'UniformOutput',false));
    caption_indx=cell2mat(cellfun(@(x) contains(x,'caption'),parameters.Stimuli.RowLabels,'UniformOutput',false));
    wordtype_indx=cell2mat(cellfun(@(x) contains(x,'WordType'),parameters.Stimuli.RowLabels,'UniformOutput',false));
    StimType_indx=cell2mat(cellfun(@(x) contains(x,'StimType'),parameters.Stimuli.RowLabels,'UniformOutput',false));
    IsRight_indx=cell2mat(cellfun(@(x) contains(x,'IsRight'),parameters.Stimuli.RowLabels,'UniformOutput',false));
    %
    trial_for_stimuli_seq=trials_value(trials_indx,:);
    trials=unique(trial_for_stimuli_seq);
    fprintf('%d trials were found \n',length(trials));
    % 
    trial_caption_seq=trials_value(caption_indx,:);
    trial_seq_cell={};
    for ii=1:length(trials)
        trial_seq_cell{ii,1}=stimuli_squence(find(trial_for_stimuli_seq==trials(ii)));
    end
    % extracting data per trial for subject
    trial_reponse=[];
    for k=1:length(trials)
        trial_indx=trial_seq_cell{k};
        % find trial type 
        wordtype=trials_value(find(wordtype_indx),trial_for_stimuli_seq==k);
        wordtype(isnan(wordtype))=[];
        if ~isempty(wordtype)
            info.word_type{k,1}=stim_types{unique(wordtype)};
        else
            info.word_type{k,1}='N/A';
        end 
        trial=struct;
        trial_index=[];
        trial_string=[];
        trial_type=[];
        stimuli_range=[];
        stimuli_type={};
        stimuli_string={};
        fprintf('adding trial %d  \n',(k))
        for kk=1:length(trial_indx)
            stimulus_index=find(states.StimulusCode==trial_indx(kk));
            stimuli_type{kk,1}=stimuli_value{caption_indx,trial_indx(kk)};
            stimuli_string{kk,1}=stimuli_value{StimType_indx,trial_indx(kk)};
            trial_index=[trial_index;stimulus_index];
            stimuli_range=[stimuli_range;[min(stimulus_index),max(stimulus_index)]];
            trial_string=[trial_string,' ',stimuli_value{caption_indx,trial_indx(kk)}];
            trial_type=[trial_type,' ',stimuli_value{StimType_indx,trial_indx(kk)}];
            %trial.(strcat('stim_',num2str(kk),'_string'))=stimuli_value{caption_indx,trial_indx(kk)};
            %trial.(strcat('stim_',num2str(kk),'_type'))=stimuli_value{StimType_indx,trial_indx(kk)};
            %trial.(strcat('stim_',num2str(kk),'_broadband'))=signal_broadband(stimulus_index,:);
            %trial.(strcat('stim_',num2str(kk),'_bandpass'))=signal_bandpass(stimulus_index,:);
            %trial.(strcat('stim_',num2str(kk),'_envelope'))=signal_envelope(stimulus_index,:);
            %trial.(strcat('stim_',num2str(kk),'_LeftEyeGazeX'))=states.EyetrackerLeftEyeGazeX(stimulus_index);
            %trial.(strcat('stim_',num2str(kk),'_LeftEyeGazeY'))=states.EyetrackerLeftEyeGazeY(stimulus_index);
            %trial.(strcat('stim_',num2str(kk),'_LeftEyePosX'))=states.EyetrackerLeftEyePosX(stimulus_index);
            %trial.(strcat('stim_',num2str(kk),'_LeftEyePosY'))=states.EyetrackerLeftEyePosY(stimulus_index);
            %trial.(strcat('stim_',num2str(kk),'_LeftPupilSize'))=states.EyetrackerLeftPupilSize(stimulus_index);

        end 
        trial.(strcat('tiral','_broadband'))=signal_broadband(trial_index,:);
        trial.(strcat('tiral','_bandpass'))=signal_bandpass(trial_index,:);
        trial.(strcat('tiral','_envelope'))=signal_envelope(trial_index,:);
        trial.(strcat('trial','_LeftEyeGazeX'))=states.EyetrackerLeftEyeGazeX(trial_index);
        trial.(strcat('trial','_LeftEyeGazeY'))=states.EyetrackerLeftEyeGazeY(trial_index);
        trial.(strcat('trial','_LeftEyePosX'))=states.EyetrackerLeftEyePosX(trial_index);
        trial.(strcat('trial','_LeftEyePosY'))=states.EyetrackerLeftEyePosY(trial_index);
        trial.(strcat('trial','_LeftPupilSize'))=states.EyetrackerLeftPupilSize(trial_index);
        trial.(strcat('trial','_string'))=trial_string;
        trial.keydown=states.KeyDown(trial_index);
        trial.keyup=states.KeyUp(trial_index);
        trial.isRight=states.IsRight(trial_index);
        trial.stimuli_range=stimuli_range-min(stimuli_range(:))+1;
        trial.stimuli_type=stimuli_type;
        trial.stimuli_string=stimuli_string;
        % find subject response: 
        index_isright_start = trial_index(find(diff(double(trial.isRight > 0)) == 1)+1);
        index_isright_stop  = trial_index(find(diff(double(trial.isRight > 0)) == -1));
        buffer_before = info.sample_rate * 1; % 1 sec 
        buffer_after  = info.sample_rate * 2; % 2 sec
        KeyDown = unique(states.KeyDown((index_isright_start-buffer_before):(index_isright_stop+buffer_after)));
        KeyDown = intersect(KeyDown,[67,77,99,109]);
        if length(KeyDown) ~= 1                    % too many key's pressed or incorrect response
           TrialResponse = 0;
        elseif KeyDown == 67 || KeyDown == 99      % response is yes (1)
           TrialResponse = 1; 
        elseif KeyDown == 77 || KeyDown == 109     % response is no  (2) 
           TrialResponse = 2; 
        else                                       % incorrect response 
           TrialResponse = 0;
        end
        % 
        trial.trial_response=TrialResponse;
        info.trial_response(k,1)=TrialResponse;
        
        dat{k,1}=trial;
        
        if ~contains(trial_type,'word')
            info.trial_type{k,1}='fixation';
        else
            info.trial_type{k,1}='word';
        end
        
    end
    info.subject=subject_name;
    info.session_name=session_name;
    eval(strcat(subject_name,'_',session_name,'.data=dat')) 
    eval(strcat(subject_name,'_',session_name,'.info=info'))
    save(strcat(subject_name,'_',session_name,'_crunched.mat'),strcat(subject_name,'_',session_name))
end

