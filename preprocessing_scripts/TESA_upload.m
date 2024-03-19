%% **** TMS–EEG signal analyser (TESA) _ 2018 / 09 ****
% Revisited 2019 / 08
% Giacomo Bertazzoli, giacomo.bertazzoli@unitn.it
% Tested on MATLAB 2017b and EEGLAB14

%% TESA PAPER:
% Rogasch NC, Sullivan C, Thomson RH, Rose NS, Bailey NW, Fitzgerald PB, Farzan F,
% Hernandez-Pavon JC. Analysing concurrent transcranial magnetic stimulation and
% electroencephalographic data: a review and introduction to the open-source TESA
% software. NeuroImage. 2017; 147: 934-951.

%% Toolbox and packages needed:
% EEGLAB with TESA extension

%% Dataset description
% Data are resting-state TMS-EEG. TMS was delivered at 4 different areas, 2
% parietal and 2 frontal (left-right Dorsolateral prefrontal cortex (lPF -
% rPF) and left-right inferior parietal lobule (lP - rP). Each area was
% stimulated 120 times with single monophasic stimuli (random ISI between
% 2-10 seconds). Total number of subjects (with test and re-test session):

%% METHOD TO REMOVE LARGE TMS-evoked MUSCLE ARTIFACT
%This script employes a FastICA decomposition with semi-automatic component
%selection to remove the TMS-evoked MUSCLE ARTIFACT (one of the alternative
%offered by TESA toolbox) Rogasch et al. 2017

%% TRIAL AND CHANNEL CHECK
% trial and channel are visually checked prior to the pipeline (as
% suggested in Rogash paper. Bad trial and Bad channel are annoted and
% saved in a structure and then removed during to the pipeline.

%% Define directories and subjects

clear   	% prepare a clean workspace
close all   % close all open events
code_dir=extractBefore(mfilename( 'fullpath' ), ['\' mfilename]); %define the folder in which the script is running
main_dir=extractBefore(code_dir, '\code'); %define the main folder

%% Add plug-in % creare cartella generale che contenga plug in con path relative
addpath([code_dir '\PLUG-INs\natsortfiles']); %natsort
addpath([code_dir '\PLUG-INs\TESA_v1.0.1']); %TESA
addpath([code_dir '\PLUG-INs\eeglab_current\eeglab2020_0']); %eeglab
addpath([code_dir '\PLUG-INs\FastICA_25']); %fastICA

% initialize eeglab;
eeglab nogui %initialize eeg lab

%% Import BIDS data with EEGlab (it transform each brain-vison dataset in a .set dataset. It takes time).
% pop_importbids(main_dir);

%% Variables
activity_prob=5; %[max] absolute thresold or activity probability limit(s) (in std. dev.) for automatic channel rejection
epoching_long=[-1 1]; %epoching in sec
demeaning_interval=[-1000 999]; %interval in ms for demeaning
downsample=1000; %Hz to reach for the downsampling
tr_rej_stnddev_first=5; %number of stnd dev to consider when computing joint probability and remove trial (trial rejection)
baseline_long=[-1000 -2]; %full period used for baseline correction
trigger_original= 16; %triggers (if more than 1, es. [132; 122; 1]
interp_interval=[-1 6]; %interval where to cut and interpolate TMS pulse
low_pass_filt=90; %low pass filter limit in Hz
hp_filt=1; %high pass filter cut-off fr in Hz
notch_filt=[48 52]; %notch filter limit in Hz
final_ref=[]; %reref for the last image. [] for avg reref

%% DEFINE DATA FOLDER
datafolder_parent = main_dir;  % parent data folder. It includes in it and all the subdir all the dataset to be processed
mkdir(main_dir, '\derivatives'); % build output folder

%% choose datasets to process
all_datasets_list = dir([datafolder_parent '\**\*.vhdr']);
all_datasets_list = all_datasets_list(~endsWith({all_datasets_list.name}, 'EI_eeg.vhdr')); % remove rest eye open initial
all_datasets_list = all_datasets_list(~endsWith({all_datasets_list.name}, 'EF_eeg.vhdr')); % remove rest eye open final

%% order the datasets and isolate names of the subjects
all_datasets_list = table2struct(sortrows(struct2table(all_datasets_list), 'name')); %alphabetic order
all_datasets_name =  extractBefore({all_datasets_list.name}, '.'); %make a list of just the name of the dataset in alphabetic order

all_subjects_name =  extractBefore({all_datasets_list.name}, 7); %make a list of just the subjects name
all_subjects_name =  unique(all_subjects_name,'sorted'); %remove string with the same name

%% Main loop
for maincounter=1:length(all_datasets_list) %loops through all datasets
    
    %build save name of the dataset
    current_datasets = all_datasets_name{maincounter}; %current dataset name
    curreent_session= extractBetween(current_datasets, 9,10); %current dataset session (T0, T1... etc)
    current_datasets_savename= all_datasets_name{maincounter}; %name with the bids structure
    
    %generate output folder in BIDS structure
    mkdir([main_dir '\derivatives\' extractBefore(current_datasets_savename, '_ses') '\'  'ses-' curreent_session{1} '\eeg\TESA\' current_datasets_savename]); %folder where to save preprocessing intermediates;
    
    %save path for intermediates and figures
    current_output_folder= [main_dir '\derivatives\' extractBefore(current_datasets_savename, '_ses') '\'  'ses-' curreent_session{1} '\eeg\TESA\' current_datasets_savename];
    
    %start saving workspace
    diary([current_output_folder '\' current_datasets_savename '_workspace_basic_prep']);
    
    %print which datasets it's being processed
    fprintf('\n******\n Processing %s\n******\n\n', ...
        current_datasets);
    
    %% load dataset
    %check the name of the current dataset in the list of .vhdr. Load the correct .vhdr from the corrisponding folder in all_datasets_list
    EEG = pop_loadbv(all_datasets_list(maincounter).folder,...
        [current_datasets '.vhdr']);
    
    %% Create a copy in the dataset in the intermediate folder
    pop_saveset(...
        EEG,...
        'filename', [current_datasets_savename '_TESA.set'],...
        'filepath', current_output_folder);
    
    %% Read in the channel location file and get channel location
    EEG = eeg_checkset( EEG );
    EEG = pop_chanedit(EEG, 'lookup','standard-10-5-cap385.elp'); % Add the channel locations. Automatic eeglab channel setup
    
    %% remove EOG electrods
    EEG=pop_select(EEG, 'nochannel', {'HEOG', 'VEOG'});
    
    %%  Remove bad electrodes
    % Save the original EEG locations for use in interpolation later
    EEG.allchan = EEG.chanlocs;
    EEG = pop_rejchan(EEG, 'elec', 1:size(EEG.data,1) ,'threshold',activity_prob,'norm','on','measure','kurt');
    
    %change set name and save this step on disk
    EEG.setname = [current_datasets_savename '_ChanRem'];
    EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);
    
    
    %% Epoch data
    % epoch data around the trigger S127 (TMS pulse) , -1000ms
    % +1000ms
    EEG = pop_epoch( EEG, {  ['S ',num2str(trigger_original)]  }, epoching_long , 'epochinfo', 'yes');
    
    %change set name and save this step on disk
    EEG.setname = [current_datasets_savename '_epochs'];
    EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);
    
    %% Demeaning Correction
    %apply demeaning, subtraction of the mean voltage of the whole
    %epoch to each point in the epoch. This steps increase the
    %reliability of the ICA
    EEG = pop_rmbase( EEG, demeaning_interval);
    
    %change set name and save this step on disk
    EEG.setname = [current_datasets_savename '_Demeaning'];
    EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);
    
    
    %% Remove TMS pulse artifact
    % remove data around the TMS pulse
    EEG = pop_tesa_removedata( EEG, interp_interval);
    
    %change set name and save this step on disk
    EEG.setname = [current_datasets_savename '_TMSpulseREM'];
    EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);
    
    %% Interpolate missing data around TMS pulse
    % interpolation of missing data aroun TMS pulse using a cubic
    % interpolation and 1 ms before an after the removed window
    EEG = pop_tesa_interpdata( EEG, 'cubic', [1,1] );
    
    %change set name and save this step on disk
    EEG.setname = [current_datasets_savename '_TMSpInt'];
    EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);
    
    %% Downsample data
    % now that the TMS pulse artifact is removed it is possible to
    % demeaning without inserting aliasing artifacts
    EEG = pop_resample( EEG, downsample);
    
    EEG.setname = [current_datasets_savename '_DownSamp'];
    EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);
    
    %% Remove bad trials
    EEG = pop_jointprob(EEG,1, 1:size(EEG.data,1),tr_rej_stnddev_first,tr_rej_stnddev_first,0,0); %automatic selection of bad epochs
    disp(['number of rejected epochs: ' find(EEG.reject.rejjp)]) %print numebr of rejected epochs
    bad_epochs=find(EEG.reject.rejjp); %create a new variables containing the indeces of bad epochs and store it in the next step
    save([current_output_folder '\bad_epochs_' current_datasets_savename], 'bad_epochs');
    EEG = pop_rejepoch( EEG, EEG.reject.rejjp, 0); %reject marked epochs
    
    %% Replace interpolated data around TMS pulse with constant amplitude data (-1 to 6 ms)
    % it is because interpolated data should not be present when
    % performing an ICA (because it is like you are adding
    % information to the dataset)
    EEG = pop_tesa_removedata( EEG, interp_interval );
    
    EEG.setname = [current_datasets_savename '_TMSpulse0'];
    EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);
    
    %% Remove TMS-evoked muscle activity (using FastICA and semi-auto component selection)
    EEG = pop_tesa_fastica( EEG, 'approach', 'symm', 'g', 'tanh', 'stabilization', 'off' );
    
    EEG.setname = [current_datasets_savename '_ICA-TMSmuscle'];
    EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);
    
    
    %parameters for automatic component selection
    
    EEG = pop_tesa_compselect( EEG, ...
        'compCheck','on',...
        'comps', 15, ...
        'figSize','small',...
        'plotTimeX',[-100 399],...
        'plotFreqX',[1 100],...
        'tmsMuscle','on',...
        'tmsMuscleThresh',8,...
        'tmsMuscleWin',[11 30],...
        'tmsMuscleFeedback','off',...
        'blink','off',...
        'blinkThresh',2.5,...
        'blinkElecs',{'Fp1','Fp2'},...
        'blinkFeedback','off',...
        'move','off',...
        'moveThresh',2,...
        'moveElecs',{'F7','F8'},...
        'moveFeedback','off',...
        'muscle','off',...
        'muscleThresh',0.6,...
        'muscleFreqWin',[30 100],...
        'muscleFeedback','off',...
        'elecNoise','off',...
        'elecNoiseThresh',4,...
        'elecNoiseFeedback','off' );
    
    
    EEG.setname = [current_datasets_savename '_CompSel1'];
    EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);
    
    
    %% load dataset afeter first component selection
    
    EEG = pop_loadset(...
        'filename', [current_datasets_savename  '_CompSel1.set'],...
        'filepath',  current_output_folder);
    
    %%  data removal (-1 to 6 ms)
    % this step is suggested from the TESA script example, even if
    % i don't understand why it is needed. I suppose is to be sure
    % there are high frequency artifact at the edge of the removed
    % window
    EEG = pop_tesa_removedata( EEG, interp_interval );
    
    EEG.setname = [current_datasets_savename '_TMSp0ext'];
    EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);
    
    
    %% Interpolate missing data around TMS pulse
    %interpolate back remove data that were substituted with
    %costant values. This is because for filtering it is better to
    %do not have sharp edges in the epoch.
    EEG = pop_tesa_interpdata( EEG, 'cubic', [5,5] );
    
    EEG.setname = [current_datasets_savename '_TMSp0Int'];
    EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);
    
    
    %% Bandpass (1-90 Hz) and bandstop (48-52 Hz) filter data
    EEG = pop_tesa_filtbutter( EEG, hp_filt, low_pass_filt, 4, 'bandpass' );
    EEG = pop_tesa_filtbutter( EEG, notch_filt(1), notch_filt(2), 4, 'bandstop' );
    
    EEG.setname = [current_datasets_savename '_notch&90_1_BandpassFilter'];
    EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);
    
    %% Replace interpolated data around TMS pulse with constant amplitude data (-1 to 6 ms)
    EEG = pop_tesa_removedata( EEG, interp_interval );
    
    EEG.setname = [current_datasets_savename '_TMSp0'];
    EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);
    
    %% Remove all other artifact (using FastICA and semi-auto component selection)
    EEG = pop_tesa_fastica( EEG, 'approach', 'symm', 'g', 'tanh', 'stabilization', 'off' );
    
    EEG.setname = [current_datasets_savename '_ICA2'];
    EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);
    
    EEG = pop_tesa_compselect( EEG,...
        'compCheck','on',...
        'comps',[],...
        'figSize','medium',...
        'plotTimeX',[-100 399],...
        'plotFreqX',[1 100],...
        'tmsMuscle','on',...
        'tmsMuscleThresh',8,...
        'tmsMuscleWin',[11 30],...
        'tmsMuscleFeedback','off',...
        'blink','on',...
        'blinkThresh',2.5,...
        'blinkElecs',{'Fp1','Fp2'},...
        'blinkFeedback','off',...
        'move','on',...
        'moveThresh',2,...
        'moveElecs',{'F7','F8'},...
        'moveFeedback','off',...
        'muscle','on',...
        'muscleThresh',0.6,...
        'muscleFreqWin',[30 100],...
        'muscleFeedback','off',...
        'elecNoise','on',...
        'elecNoiseThresh',4,...
        'elecNoiseFeedback','off' );
    
    
    EEG.setname = [current_datasets_savename '_CompSel2'];
    EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);
    
    
    %% load dataset afeter second component selection
    EEG = pop_loadset(...
        'filename', [current_datasets_savename  '_CompSel2.set'],...
        'filepath', current_output_folder);
    
    %% Interpolate missing data around TMS pulse
    EEG = pop_tesa_interpdata( EEG, 'cubic', [5,5] );
    
    EEG.setname = [current_datasets_savename '_After2ndCompSel_TMSpInt'];
    EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);
    
    %% Interpolate missing channels
    %channel removed because were bad are now interpolated
    EEG = pop_interp(EEG, EEG.allchan, 'spherical');
    
    EEG.setname = [current_datasets_savename '_ChanInt'];
    EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);
    
    %% Re-reference to average of all electrodes
    %Data were recorded with online reference in Tp9.
    EEG = pop_reref( EEG, final_ref);
    
    EEG.setname = [current_datasets_savename '_AvgReref'];
    EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);
    
    %% Baseline Correction (-1000 ms to -2 ms)
    EEG = pop_rmbase( EEG, baseline_long);
    
    EEG.setname = [current_datasets_savename '_baselineCorr'];
    EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);
    
    
    %% Save point
    EEG.setname= [current_datasets_savename 'TESA_Processed'];
    EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', current_output_folder);
    
    %% Plot the results
    % final TEP  scaled
    
    figure
    plot(EEG.times,mean(EEG.data, 3)); %plot avg trial TEP
    hold on
    xlim([-100 400]);%restrict epoch
    ylim([-20 30]); % restrict amplitude range
    title('Final scaled avg reref');
    
    hold off
    %save figure
    saveas(gcf,...  %save .fig
        [current_output_folder '\' current_datasets_savename '_Final_scaled_avg_ref'])
    saveas(gcf,... %save .png
        [current_output_folder '\' current_datasets_savename '_Final_scaled_avg_ref'],...
        'png')
    close all
    diary off %stop saving workspace
    
end