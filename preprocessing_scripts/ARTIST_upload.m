%% ARTIST %%
%% ****Automated aRTIfact rejection for Single-pulse TMS-EEG Data (ARTIST)_
% 2018 / 09 ****
% Revisited 2019 / 08
% Giacomo Bertazzoli, giacomo.bertazzoli@unitn.it
% Tested on MATLAB 2017b and EEGLAB14

%% ARTIST PAPER:
%Wu W, Keller C, Rogasch N, Longwell P, Shpigel E, Rolle C, Etkin A.
%ARTIST: a fully automated artifact rejection algorithm for single-pulse
%TMS-EEG data. Human Brain Mapping, 2018, DOI: 10.1002/hbm.23938.

%% Toolbox and packages needed:
% Wavelet Toolbox
% EEGLAB
% ARTIST pakage

%% Dataset description
% Data are resting-state TMS-EEG. TMS was delivered at 4 different areas, 2
% parietal and 2 frontal (left-right Dorsolateral prefrontal cortex (lPF -
% rPF) and left-right inferior parietal lobule (lP - rP). Each area was
% stimulated 120 times with single monophasic stimuli (random ISI between
% 2-10 seconds). Total number of subjects (with test and re-test session):

%% TRIAL AND CHANNEL CHECK
%  ARTIST is fully automated

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
addpath([code_dir '\PLUG-INs\ARTISTforDistribution_ExternalLabs']); %ARTIST

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
    mkdir([main_dir '\derivatives\' extractBefore(current_datasets_savename, '_ses') '\'  'ses-' curreent_session{1} '\eeg\ARTIST\' current_datasets_savename]); %folder where to save preprocessing intermediates;
    
    %save path for intermediates and figures
    current_output_folder= [main_dir '\derivatives\' extractBefore(current_datasets_savename, '_ses') '\'  'ses-' curreent_session{1} '\eeg\ARTIST\' current_datasets_savename];
    
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
        'filename', [current_datasets_savename '_ARTIST.set'],...
        'filepath', current_output_folder);
    
    %% Read in the channel location file and get channel location
    EEG = eeg_checkset( EEG );
    EEG = pop_chanedit(EEG, 'lookup','standard-10-5-cap385.elp'); % Add the channel locations. Automatic eeglab channel setup
    
    %% remove EOG electrods
    EEG=pop_select(EEG, 'nochannel', {'HEOG', 'VEOG'});
    
    
    %% cfg parameters
    cfg.EventCode= ['S ' num2str(trigger_original)]; % - Event code for the TMS pulses in cfg.EventCode (needs to be a string or a numeric value)
    cfg.TMSEEGrootFolder = current_output_folder; % - Path of the folder to store results
    % Optional fields:
    cfg.TrialStart= epoching_long(1)*1000; % - Start of each epoch in ms (default: -1000)
    cfg.TrialEnd= epoching_long(2)*1000; %- End of each epoch in ms (default: 1500)
    cfg.PulseLen= interp_interval(2); % - Length of the TMS pulse to be removed and interpolated in ms (default: 10)
    cfg.PulseShift= interp_interval(1); %- Time (ms) from zero where you want to start removing the TMS pulse (make sure it's before the rise of the pulse. Default: -2)
    cfg.BaseLine= baseline_long; % - Baseline time window in ms (default: [-300, -100])
    cfg.plottimes= [15 30:20:300]; % - Times (ms) at which topo plots are shown in the butterfly plots (default: [15,25,40,60,75,100,150,200,300])
    cfg.NameProject='Analyses of TMS-EEG signal'; %- Name of the project (default: '')
    cfg.NameCond = ''; %- Name of the condition (default: '')
    cfg.NameSub= current_datasets_savename; % - Name of the subject (default: '')
    
    
    %% ARTIST main function
    % the main function has some chenges to make it match to the other
    % pipelines of this study. Check the main function to see the changes
    % the function savefig in the ARTIST package was renamed savefig_ARTIST
    % because it interferes with matlab savefig. Origninal ARTIST function
    % is in the same folder
    EEG = ARTISTMain_GIACOMO(EEG,cfg);
    
    %% Save point
    %% save dataset after being loaded
    EEG.setname= [current_datasets_savename '_ARTIST'];
    
    pop_saveset(EEG, 'filename',[current_datasets_savename '_ARTIST'],'filepath', current_output_folder);
    
    
    %% save rejections
    %form a strcuture containing all the removed chan trials and ICcomp
    
    ARTIST_badstuff.name=current_datasets_savename;
    ARTIST_badstuff.artchan=EEG.artchan;
    ARTIST_badstuff.numartchan=length(EEG.artchan);
    ARTIST_badstuff.arttrial=EEG.arttrial;
    ARTIST_badstuff.artcomp1=EEG.artcomp1;
    ARTIST_badstuff.artcomp2=EEG.artcomp2;
    
    
    save([current_output_folder '\ARTIST_badstuff'], 'ARTIST_badstuff');
    
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

disp('All preprocessing ended')
