%% **** TMSEEG _ 2018 / 09 ****
% Revisited 2019 / 08
% Giacomo Bertazzoli, giacomo.bertazzoli@unitn.it
% Tested on MATLAB 2017b and EEGLAB14

%% TOOLBOX PRESENTED IN THIS PAPER:
%Atluri S, Frehlich M, Mei Y, GarciaDominguez L, Rogasch NC, Wong
%W,Daskalakis ZJ and Farzan F (2016)TMSEEG: A MATLAB-BasedGraphical User
%Interface forProcessing ElectrophysiologicalSignals during Transcranial
%MagneticStimulation.Front. Neural Circuits 10:78.

%% Toolbox and packages needed:
% EEGLAB
% TMSEEG app for matlab
% natsort
% TESA
% fastICA

%% BUGs
% it happened once that the trigger was misplaced by the function, causing
% a shift in the trials. Probably an error due to the too advanced version
% of matalb employd in the beginning.
% Sometimes it gived errors in the remove channel/trial pahse to to
% interaction with other function in the path (in my case it was fixed by
% removing fieldtrip from the matlab path).

%% Dataset description
% Data are resting-state TMS-EEG. TMS was delivered at 4 different areas, 2
% parietal and 2 frontal (left-right Dorsolateral prefrontal cortex (lPF -
% rPF) and left-right inferior parietal lobule (lP - rP). Each area was
% stimulated 120 times with single monophasic stimuli (random ISI between
% 2-10 seconds). Total number of subjects (with test and re-test session):

%% TRIAL AND CHANNEL CHECK
% fully manual

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
addpath([code_dir '\PLUG-INs\TMSEEG-4.0']); %TMSEEG-4.0
addpath([code_dir '\PLUG-INs\EEGDataPro-1.0']); %EEGDataPro-1.0

% initialize eeglab;
eeglab nogui %initialize eeg lab

%% Import BIDS data with EEGlab (it transform each brain-vison dataset in a .set dataset. It takes time).
% pop_importbids(main_dir);

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
    curreent_session= extractBetween(current_datasets, 15,16); %current dataset session (T0, T1... etc)
    current_datasets_savename= all_datasets_name{maincounter}; %name with the bids structure
    
    %generate output folder in BIDS structure
    mkdir([main_dir '\derivatives\' extractBefore(current_datasets_savename, '_ses') '\'  'ses-' curreent_session{1} '\eeg\TMSEEG\' current_datasets_savename]); %folder where to save preprocessing intermediates;
    
    %save path for intermediates and figures
    current_output_folder= [main_dir '\derivatives\' extractBefore(current_datasets_savename, '_ses') '\'  'ses-' curreent_session{1} '\eeg\TMSEEG\' current_datasets_savename];
    
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
        'filename', [current_datasets_savename '_TMSEEG.set'],...
        'filepath', current_output_folder);
    
    
    %% Open tmseeg interface
    %open tmseeg interface and block the script
    cd(current_output_folder)
    tmseeg_main
    keyboard %block the script while doing the cleaning 
    %press continue on the top >> to continue the script

    %% load cleaned set
    %after tmseeg toolbox cleaning a final .set file is generated.
    %load that file and ultimate the cleaning.
    EEG = pop_loadset(...
        'filename',[current_datasets_savename '_TMSEEG_1_2_3_4_5_6_7_8_9_10.set'],...
        'filepath', current_output_folder);
    
    %% deals with data gap left from TMS pulse removal
    %tmseeg fills the TMS pulse removed data with NaN. Using TESA
    %function to remove the NaN and replace it with interpolate
    %data.
    
    %find indicies of interval with NaNs
    NaN_interval=find(isnan(EEG.data(1,:,1)));
    interval_to_replace=[min(EEG.times(NaN_interval)) round(max(EEG.times(NaN_interval)))]; %find time in ms where are the NaNs
    save([current_output_folder '\interp_interval' ], 'NaN_interval' , 'interval_to_replace')
    EEG = pop_tesa_removedata( EEG, interval_to_replace ); %remove data. Needed to make use of pop_tesa_interpdata
    EEG = pop_tesa_interpdata( EEG, 'cubic', [5,5] );
    
    
    %% Baseline Correction (-1000 ms to -2 ms)
    %the baseline correction during the guided procedure in tmseeg
    %is used as a demeaning (remove the average of the whole epoch
    %to each epoch point) to improve ICA reliability. This baseline
    %correction is for visualization porpouses.
    EEG = pop_rmbase( EEG, [-1000  -2]);
    
    %% Downsample to 1000Hz
    % during the procedure the downsalple is done BEFORE the
    % removal of the TMS-pulse artifact, hance bringing a ringing
    % artifact before and after the pulse. To avoid that I have
    % skipped that step during the guided procedure and I am doing
    % it now.
    EEG = pop_resample(EEG, 1000);
    %% Save point
    EEG.setname= [current_datasets_savename '_TMSEEG_Processed'];
    EEG = pop_saveset(EEG,...
        'filename', [EEG.setname '.set'],...
        'filepath', current_output_folder);
    
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