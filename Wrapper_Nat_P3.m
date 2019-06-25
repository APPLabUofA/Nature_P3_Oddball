%% PREPROCESSING AND ANALYSIS WRAPPER
%clear and close everything
clear all;
close all;
clc;

%% Settings for loading the raw data
%Datafiles must be in the format exp_participant, e.g. EEGexp_001.vhdr
exp.name = 'p3_nature';
exp.participants = {'012';'013';'014';'015';'016';'017';'018';'019';'020';'021';'022';'023';'024';'025';'026';'027';'028';'029'};
exp.pathname = 'M:\Data\Nature_P3\eeg';

% Choose how to organize your datasets.
% The sets are different trial types that contain comparable events.
% Examples of different sets: 10Hz/12Hz stimulation on different trials; emotional vs non-emotional images on different trials.
% Each row is a different set.
% The events are stimuli within the trial.
% Examples of different events. Missing vs present targets, valid vs invalid targets, detected vs undetected responses.
% Each column is a different event.
% You can collect multiple triggers into one event with square brackets [].

triggers = 2; %%1 for auditory tones; 2 for visual scenes
preprocessing = 1; %%%1 to only run preprocessing, 0 to only run analysis

if triggers == 1
    %%%%%Triggers to tone onset%%%%%
    
    % Each item in the exp.events matrix will become a seperate dataset, including only those epochs referenced by the events in that item.
    %e.g. 3 rows x 4 columns == 12 datasets/participant
    exp.events = {[11],[12];...
        [21],[22];...
        [31],[32]};    %must be matrix (sets x events)
    
    exp.setname = {'Baseline';'Nature';'Urban'}; %name the rows
    exp.event_names = {'Low_Tones','High_Tones'}; %name the columns
    
    % The settings will be saved as a new folder. It lets you save multiple datasets with different preprocessing parameters.
    % exp.settings = 'P3_Nature';%%%Tone Onset
    
    exp.settings = 'P3_Nature_DELETE';%%%Picture Onset
    
    exp.epochslims = [-1 1.2]; %Tone Onset in seconds; epoched trigger is 0 e.g. [-1 2] USED [-2 4]
    
    % % Blink Correction
    %the Blink Correction wants dissimilar events (different erps) seperated by commas and similar events (similar erps) seperated with spaces. See 'help gratton_emcp'
    exp.selection_cards = {'11','21','31','12','22','32'};
    
elseif triggers == 2
    %%%%Triggers to picture onset (grouped)%%%%%
    
    % Each item in the exp.events matrix will become a seperate dataset, including only those epochs referenced by the events in that item.
    %e.g. 3 rows x 4 columns == 12 datasets/participant
    exp.events = {[999];...
        [101:150];...
        [201:250]};    %must be matrix (sets x events)
    
    exp.setname = {'Baseline';'Nature';'Urban'}; %name the rows
    exp.event_names = {'Pictures'}; %name the columns
    
    % The settings will be saved as a new folder. It lets you save multiple datasets with different preprocessing parameters.
    exp.settings = 'P3_Nature_Pictures_DELETE';%%%Picture Onset
    
    exp.epochslims = [-1 2.5]; %Picture Onset in seconds; epoched trigger is 0 e.g. [-1 2]
    
    % % Blink Correction
    %the Blink Correction wants dissimilar events (different erps) seperated by commas and similar events (similar erps) seperated with spaces. See 'help gratton_emcp'
    exp.selection_cards = {'999',num2str(101:150),num2str(201:250)};
    
end
%% Preprocessing Settings
%segments settings
exp.segments = 'on'; %Do you want to make new epoched datasets? Set to "off" if you are only changing the tf settings.
%Where are your electrodes? (.ced file)
exp.electrode_locs = 'M:\Analysis\ElectrodeLocs\EOG-electrode-locs-32.ced';

%% Filter the data?
exp.filter = 'on';
exp.lowpass = 50;
exp.highpass = 0;

%% Re-referencing the data
exp.refelec = 1; %which electrode do you want to re-reference to?
exp.brainelecs = [2:32]; %list of every electrode collecting brain data (exclude mastoid reference, EOGs, HR, EMG, etc.

%% Epoching the data
%Choose what to epoch to. The default [] uses every event listed above.
%Alternatively, you can epoch to something else in the format {'TRIG'}. Use triggers which are at a consistent time point in every trial.
exp.epochs = [];
exp.epochbaseline = [-200 0]; %remove the for each epoched set, in ms. e.g. [-200 0]

%% Artifact rejection.
% Choose the threshold to reject trials. More lenient threshold followed by an (optional) stricter threshold
exp.preocularthresh = [-1000 1000]; %First happens before the ocular correction.
exp.postocularthresh = [-500 500]; %Second happens after. Leave blank [] to skip

%% Time-Frequency settings
%Do you want to run time-frequency analyses? (on/off)
exp.tf = 'off';
%Do you want to save the single-trial data? (on/off) (Memory intensive!!!)
exp.singletrials = 'off';
%Do you want to use all the electrodes or just a few? Leave blank [] for all (will use same as exp.brainelecs)
exp.tfelecs = [];
%Saving the single trial data is memory intensive. Just use the electrodes you need.
exp.singletrialselecs = [3];

%% Wavelet settings
%how long is your window going to be? (Longer window == BETTER frequency resolution & WORSE time resolution)
exp.winsize = 512; %in ms; use numbers that are 2^x, e.g. 2^10 == 1024ms
%baseline will be subtracted from the power variable. It is relative to your window size.
exp.erspbaseline = [-1000 0]; %e.g., [-200 0] will use [-200-exp.winsize/2 0-exp.winsize/2]; Can use just NaN for no baseline
%Instead of choosing a windowsize, you can choose a number of cycles per frequency. See "help popnewtimef"
exp.cycles = [0]; %leave it at 0 to use a consistent time window
exp.freqrange = [1 40]; % what frequencies to consider? default is [1 50]
%%%%

%% Save your pipeline settings
save([exp.settings '_Settings'],'exp') %save these settings as a .mat file. This will help you remember what settings were used in each dataset

if preprocessing == 1
    % % Run Preprocessing
    Preprocessing_Nat_P3(exp) %comment out if you're only doing analysis
    
elseif preprocessing == 0
    % % Run Analysis
    %Don't want to change all the above settings? Load the settings from the saved .mat file.
    
    %choose the data types to load into memory (on/off)
    anal.segments = 'on'; %load the EEG segments?
    anal.tf = 'off'; %load the time-frequency data?
    
    anal.singletrials = 'off'; %load the single trial data?
    anal.entrainer_freqs = [20; 15; 12; 8.5; 4]; %Single trial data is loaded at the event time, and at the chosen frequency.
    
    anal.tfelecs = []; %load all the electodes, or just a few? Leave blank [] for all.
    anal.singletrialselecs = [2 3 4 6];
    
    Analysis_Nat_P3(exp,anal) % The Analysis primarily loads the processed data. It will attempt to make some figures, but most analysis will need to be done in seperate scripts.
end
