%% PREPREOCESSING raw TMS-EEG data with SOUND and SSP-SIR

% Origninal description from Tuomas Mutanen:
% This is an EEGLAB-compatible MATLAB script that demonstrates how DDWiener
% and SOUND can be used to clean a partially noisy EEG dataset.
%
% Required toolbox to run this script: EEEGlab, BioSig, TESA, bva-io (plugin to import BrainVision data)
% 18 September 2017 : Tuomas Mutanen and Marta Bortoletto.
% Department of Neuroscience and Biomedical Engineering (NBE), School of Science, Aalto University

%original papers:
% Mutanen, T. P., Metsomaa, J., Liljander, S., & Ilmoniemi, R. J. (2018).
% Automatic and robust noise suppression in EEG and MEG: The SOUND
% algorithm. Neuroimage, 166, 135-151.

% Mutanen, T. P., Kukkonen, M., Nieminen, J. O., Stenroos, M., Sarvas, J.,
% & Ilmoniemi, R. J. (2016). Recovering TMS-evoked EEG responses masked by
% muscle artifacts. Neuroimage, 139, 157-166.

% .........................................................................

% UPDATE Giacomo Bertazzoli 14th october 2019
% UPDATE Giacomo Bertazzoli November 2019
% UPDATE Giacomo Bertazzoli December 2019
% Needs Fuzzy Logic Toolbox
% UPDATE Giacomo Bertazzoli, Marta Bortoletto, Guido Barchiesi, December 2020
% UPDATE Giacomo Bertazzoli, Marta Bortoletto, Guido Barchiesi, Jenuary 2021 MATLAB version R2020b
% UPDATE Giacomo Bertazzoli, Marta Bortoletto, April 2021 MATLAB version R2020b
% UPDATE Giacomo Bertazzoli, Marta Bortoletto, June 2021 MATLAB version R2020b

%% Define directories and subjects

clear   	% prepare a clean workspace
close all   % close all open events
code_dir=extractBefore(mfilename( 'fullpath' ), ['\' mfilename]); %define the folder in which the script is running
main_dir=extractBefore(code_dir, '\code'); %define the main folder

%% Add plug-in % creare cartella generale che contenga plug in con path relative
addpath([code_dir '\PLUG-INs\SOUND']);%SOUND
addpath([code_dir '\PLUG-INs\natsortfiles']); %natsort
addpath([code_dir '\PLUG-INs\TESA1.1.1']); %TESA
addpath([code_dir '\PLUG-INs\eeglab_current\eeglab2020_0']); %eeglab

% initialize eeglab;
eeglab nogui %initialize eeg lab

%% Import BIDS data with EEGlab (it transform each brain-vison dataset in a .set dataset. It takes time).
% pop_importbids(main_dir);

%% Variables

epoching_long=[-1 1]; %epoching in sec
epoching_short=[-0.1 0.350]; %final epoching for short version of TEPs
downsample=1000; %Hz to reach for the downsampling 
baseline_long=[-1000 -2]; %full period used for baseline correction 
trigger_original= 16; %triggers (if more than 1, es. [132; 122; 1]
interp_interval=[-1 6]; %interval where to cut and interpolate TMS pulse
ssp_sir_timerange=[-1 50]; % interval in which sspsir estimates the muscolar artifact
low_pass_filt=90; %low pass filter limit in Hz
hp_filt=1; %high pass filter cut-off fr in Hz
notch_filt=[48 52]; %notch filter limit in Hz 
ref_channel= 'TP9'; %ref for data and spheric lead field before SOUND 
final_ref=[]; %reref for the last image. [] for avg reref 
max_SSP_SIR_PC_to_remove=6; % number of PC removed sequentially during SSP-SIR step. 6 means that SSP-SIR will output 6 datasets in which SSP-SIR removed 0, 1,2,... 5 PCs. From this pool of dataset you will choose the correct number of PC to remove to cancel the muscle artifact

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

%% 1. BASIC PREPROCESSING
for maincounter=1:length(all_datasets_list) %loops through all datasets
    
    %build save name of the dataset
    current_datasets = all_datasets_name{maincounter}; %current dataset name
    curreent_session= extractBetween(current_datasets, 15,16); %current dataset session (T0, T1... etc)
    current_datasets_savename= all_datasets_name{maincounter}; %name with the bids structure
    
    %generate output folder in BIDS structure
    mkdir([main_dir '\derivatives\' extractBefore(current_datasets_savename, '_ses') '\'  'ses-' curreent_session{1} '\eeg\SOUND_SSP-SIR\' current_datasets_savename]); %folder where to save preprocessing intermediates;
    
    %save path for intermediates and figures
    current_output_folder= [main_dir '\derivatives\' extractBefore(current_datasets_savename, '_ses') '\'  'ses-' curreent_session{1} '\eeg\SOUND_SSP-SIR\' current_datasets_savename];
    
    %start saving workspace
    diary([current_output_folder '\' current_datasets_savename '_workspace_basic_prep']);
    
    %print which datasets it's being processed
    fprintf('\n******\n Processing %s\n******\n\n', ...
        current_datasets);
    
    % load dataset
    %check the name of the current dataset in the list of .vhdr. Load the correct .vhdr from the corrisponding folder in all_datasets_list
    EEG = pop_loadbv(all_datasets_list(maincounter).folder,...
        [current_datasets '.vhdr']);
    
    % Read in the channel location file and get channel location
    EEG = eeg_checkset( EEG );
    EEG = pop_chanedit(EEG, 'lookup','standard-10-5-cap385.elp'); % Add the channel locations. Automatic eeglab channel setup
    
    % remove EOG electrods
    EEG=pop_select(EEG, 'nochannel', {'HEOG', 'VEOG'});
    
    % define  EEG electrodes and reference
    [EEG.chanlocs(1:end).type] = deal('EEG');
    
    % Create lead-field matrix which allows to convert signal between sensors and sources
    % Let's make a dummy EEGLAB set for constructing a leadfield that has the correct reference:
    EEG_dummy = EEG;
    EEG_dummy.data(63,length(EEG_dummy.data)) = 0; %adds one 0 microV data channel (ref TP9)
    EEG_dummy.nbchan = 63; %change nbchan to keep consistency in number of data channels
    EEG_dummy.chanlocs(63).type = 'EEG';
    EEG_dummy.chanlocs(63).labels = ref_channel;
    EEG_dummy = pop_chanedit(EEG_dummy,  'lookup', 'standard-10-5-cap385.elp');
    [LFM_sphere] = construct_spherical_lead_field(EEG_dummy, [],[],[],[],[],[],[],1000);
    
    % save EEG orig, dataset and figure
    EEG.setname=current_datasets_savename;
    EEG = pop_saveset( EEG, 'filename',[EEG.setname '_orig'],...
        'filepath', current_output_folder );
    
    % High pass filter data
    EEG = pop_eegfiltnew(EEG, 'locutoff', hp_filt,'plotfreqz',1);
    close
    
    % epoching
    EEG = pop_epoch( EEG, {  ['S ',num2str(trigger_original)] }, epoching_long, 'epochinfo', 'yes');
    
    %   interp TMS-pulse artifact
    EEG_before_interpolation = EEG; %store dataset before interpolation
    EEG = pop_tesa_removedata( EEG, interp_interval); %remove TMS-pulse and replace with 0s
    
    % interpolation of missing data aroun TMS pulse using a cubic
    % interpolation and 1 ms before an after the removed window
    EEG = pop_tesa_interpdata( EEG, 'cubic', [1,1] );
    
    % rm baseline
    EEG = pop_rmbase(EEG, baseline_long);
    
    % Re-reference lead field to ref
    ref_chan_num=find(strcmp({EEG_dummy.chanlocs.labels}, ref_channel)); %chan num of TP9
    LFM_sphere_ref= repmat(LFM_sphere(ref_chan_num,:), size(LFM_sphere, 1)-1, 1); % for lead field reref matrix
    LFM_sphere = LFM_sphere([1:ref_chan_num-1 ref_chan_num+1:end],:) - LFM_sphere_ref; % Adjust lead-field matrix on the new reference
    
    % visualize TEP after basic preprocessing
    EEG_epoch=EEG;
    EEG_epoch= pop_rmbase(EEG_epoch, [-750 0]);
    TEP_dummy=mean(EEG_epoch.data, 3); %average across trials.
    figure
    plot(EEG_epoch.times, TEP_dummy(:,:)); %plot avg trial TEP
    hold on
    xlim([-100 400]);%restrict epoch
    ylim([-50 50]); % restrict amplitude range
    title('After  basic pre-processing');
    hold off
    saveas(gcf,...  %save .fig
        [current_output_folder '\' current_datasets_savename '_basic_prep'])
    saveas(gcf,... %save .png
        [current_output_folder '\' current_datasets_savename '_basic_prep'],...
        'png')
    close
    
    % save dataset
    EEG.setname= [current_datasets_savename '_basic_prep'];
    pop_saveset(EEG, 'filename',[current_datasets_savename '_basic_prep'],...
        'filepath', current_output_folder);
    
    % save leadfield
    save([current_output_folder '\LFM_sphere'], 'LFM_sphere'); % save leadfield
    
    disp('Basic preprocessing ended');
    
    diary off %stop saving workspace
end


%% 2. Use SOUND to clean the channel-specific noise

for maincounter=1:length(all_datasets_list) %loops through all datasets  %loops through all datasets
    
    %build save name of the dataset
    current_datasets = all_datasets_name{maincounter}; %current dataset name
    curreent_session= extractBetween(current_datasets, 15,16); %current dataset session (T0, T1... etc)
    current_datasets_savename= all_datasets_name{maincounter}; %name with the bids structure
    
    %generate output folder in BIDS structure
    mkdir([main_dir '\derivatives\' extractBefore(current_datasets_savename, '_ses') '\'  'ses-' curreent_session{1} '\eeg\SOUND_SSP-SIR\' current_datasets_savename]); %folder where to save preprocessing intermediates;
    
    %save path for intermediates and figures
    current_output_folder= [main_dir '\derivatives\' extractBefore(current_datasets_savename, '_ses') '\'  'ses-' curreent_session{1} '\eeg\SOUND_SSP-SIR\' current_datasets_savename];
    
    %start saving workspace for SOUND
    diary([current_output_folder '\' current_datasets_savename '_workspace_SOUND']);
    
    % load dataset (from previous shell)
    EEG = pop_loadset('filename',[current_datasets_savename '_basic_prep.set'],'filepath',  current_output_folder);
    
    close all
    
    % load leadfield position for each subject
    load([current_output_folder '\LFM_sphere'],  'LFM_sphere');
    
    % Set the reqularization parameter of the minimum-norm estimate when
    % finding the noise-free current estimates. This can be adjusted if clear
    % under- or overcorrection is observed. lambda_value = 0.1 was used in the original
    % article.
    lambda_value = 0.1;
    
    % Set the number of iterations when estimating the channel-specific noise
    % estimates. Iter = 5 was found to be sufficient in the datasets studied in
    % the original article.
    iter = 5;
    
    EEG_clean = EEG;
    
    chanN = size(EEG.data,1);
    
    % For loop for estimating the noise at each time sample separately:
    
    % profile on;
    
    all_sigmas = zeros(size(EEG.data,1),size(EEG.data,2));
    
    tic
    
    % parpool()
    parfor SOUND_maincounter = 1:(size(EEG.data,2))% 1:ms_samples:(size(EEG.data,2) - ms_samples + 1 )
        
        % Find the number of trials:
        N_tr= size(EEG.data,3);
        
        % Run the SOUND algorithm:
        
        % Inform about the progress of SOUND:
        disp( num2str( EEG.times(SOUND_maincounter)))
        
        % Estimate the noise for each sample
        
        tmp_data = reshape(EEG.data(:, SOUND_maincounter ,:), [EEG.nbchan, N_tr] );
        
        %if i ==1
        %   sigmas = ones(size(tmp_data,1),1);
        [y_solved, sigmas] = DDWiener(tmp_data);
        
        % end
        
        [~, ~, sigmas] = SOUND_alt(tmp_data, LFM_sphere, iter, lambda_value,sigmas,1);
        
        all_sigmas(:,SOUND_maincounter) = sigmas;
        %        disp(['processing ' all_datasets{maincounter} ' _ dataset N° ' num2str(maincounter)])
        %  profile viewer
        
    end
    
    
    % Now we clean each sample by smoothly estimating the sample-specific
    % noises with a smooth time-window
    
    scaling_w = 2;
    smoothness = 1;
    all_sigmas_ext = [repmat(all_sigmas(:,1),[1  scaling_w]), all_sigmas,  repmat(all_sigmas(:,end),[1  scaling_w])];
    step =  scaling_w + 1;
    
    smoothing_func = gaussmf(1:(2*scaling_w + 1),[smoothness  scaling_w + 1]);
    smoothing_func = smoothing_func/sum(smoothing_func); figure; plot(smoothing_func)
    smoothing_func = repmat(  smoothing_func, [size(EEG.data,1), 1]);
    
    
    for SOUND_maincounter = 1:size(EEG.data,2)% 1:ms_samples:(size(EEG.data,2) - ms_samples + 1 )
        
        % Find the number of trials:
        N_tr= size(EEG.data,3);
        
        % Make the final estimate for the noise distribution for this
        % sample, by using the noise estimate and its neigbourhood:
        
        sigmas =  sum(all_sigmas_ext(:,(step-scaling_w):(step +scaling_w) ).*smoothing_func,2);
        
        % Final cleaning step:
        
        W = diag(1./sigmas);
        WL = W*LFM_sphere;
        WLLW = WL*WL';
        
        data = reshape(EEG.data(:, SOUND_maincounter , :), [chanN, N_tr] );
        
        x = WL'*((WLLW + lambda_value*trace(WLLW)/chanN*eye(chanN))\(W*data));
        corrected_data = LFM_sphere*x;
        EEG_clean.data(:, SOUND_maincounter, :) = reshape(corrected_data,[size(corrected_data,1), 1, N_tr]);
        
        step = step + 1;
        
    end
    
    disp(['ELAPSED SOUND TIME: min ' num2str(toc/60)]);
    
    %  figure after SOUND
    EEG_epoch=EEG_clean;
    EEG_epoch= pop_rmbase(EEG_epoch, [-750 0]);
    TEP_dummy=mean(EEG_epoch.data, 3); %average across trials.
    figure
    plot(EEG_epoch.times, TEP_dummy(:,:)); %plot avg trial TEP
    hold on
    xlim([-100 400]);%restrict epoch
    ylim([-50 50]); % restrict amplitude range
    title('After SOUND');
    hold off
    saveas(gcf,...  %save .fig
        [current_output_folder '\' current_datasets_savename '_afterSOUND'])
    saveas(gcf,... %save .png
        [current_output_folder '\' current_datasets_savename '_afterSOUND'],...
        'png')
    close all
    
    
    %   save data after SOUND:
    EEG=EEG_clean; %sound saved the clean session in EEG_clean
    EEG.setname= [current_datasets_savename '_afterSOUND'];
    pop_saveset(EEG, 'filename',[current_datasets_savename '_afterSOUND'],'filepath', current_output_folder);
    
    disp('SOUND ended')
    diary off %stop saving workspace
    
end



%% 3. Artifact rejection (MANUAL)

for maincounter=1:length(all_datasets_list) %loops through all datasets
    
    %build save name of the dataset
    current_datasets = all_datasets_name{maincounter}; %current dataset name
    curreent_session= extractBetween(current_datasets, 15,16); %current dataset session (T0, T1... etc)
    current_datasets_savename= all_datasets_name{maincounter}; %name with the bids structure
    
    %generate output folder in BIDS structure
    mkdir([main_dir '\derivatives\' extractBefore(current_datasets_savename, '_ses') '\'  'ses-' curreent_session{1} '\eeg\SOUND_SSP-SIR\' current_datasets_savename]); %folder where to save preprocessing intermediates;
    
    %save path for intermediates and figures
    current_output_folder= [main_dir '\derivatives\' extractBefore(current_datasets_savename, '_ses') '\'  'ses-' curreent_session{1} '\eeg\SOUND_SSP-SIR\' current_datasets_savename];
    
    %start saving workspace for SOUND
    diary([current_output_folder '\' current_datasets_savename '_workspace_TrailRem']);
    
    
    % load dataset (from previous shell)
    EEG = pop_loadset('filename',[current_datasets_savename '_afterSOUND.set'],'filepath',  current_output_folder);
    
    % Artifact rejection to manually reject epochs
    % (a window pops up) scroll through epochs and select trials to be rejected by
    % clicking on the bad epoch (they turn into yellow). At the end press
    % 'reject' button.
    %mark bad epochs without rejection. I do not reject epochs here because
    %i need to save the indecies of bad epochs. Direct rejection apparently
    %does not allow to to that
    
    eeglab redraw
    close all
    pop_eegplot( EEG,1,1,0); %plot data scroll
    waitfor( findobj('parent', gcf, 'string', 'UPDATE MARKS'), 'userdata'); %select trial to remove and wait for the user to press UPDATE MARKS
    close all
    disp(['index  of rejected epochs: ' num2str(find(EEG.reject.rejmanual))]) %print numebr of rejected epochs
    rejected_epochs=find(EEG.reject.rejmanual); %create a new variables containing the indeces of bad epochs and store it in the next step
    save([current_output_folder '\rejected_epochs'], 'rejected_epochs');
    EEG = pop_rejepoch( EEG, EEG.reject.rejmanual, 0); %reject marked epochs
    
    % visualize TEP after after rej
    EEG_epoch=EEG;
    EEG_epoch= pop_rmbase(EEG_epoch, [-750 0]);
    TEP_dummy(:,:)=mean(EEG_epoch.data, 3); %average across trials.
    figure
    plot(EEG_epoch.times, TEP_dummy(:,:)); %plot avg trial TEP
    hold on
    xlim([-100 400]);%restrict epoch
    ylim([-50 50]); % restrict amplitude range
    title('After trialrej');
    hold off
    saveas(gcf,...  %save .fig
        [current_output_folder '\' current_datasets_savename '_trialrej'])
    saveas(gcf,... %save .png
        [current_output_folder '\' current_datasets_savename '_trialrej'],...
        'png')
    
    close all
    
    % save data after trials rejection:
    EEG.setname= [current_datasets_savename '_trialrej'];
    pop_saveset(EEG, 'filename',[current_datasets_savename '_trialrej'],'filepath', current_output_folder);
    
    
    disp('Trialrej ended')
    diary off %stop saving workspace
end

%% 4. ICA

for maincounter=1:length(all_datasets_list) %loops through all datasets
    
    %build save name of the dataset
    current_datasets = all_datasets_name{maincounter}; %current dataset name
    curreent_session= extractBetween(current_datasets, 15,16); %current dataset session (T0, T1... etc)
    current_datasets_savename= all_datasets_name{maincounter}; %name with the bids structure
    
    %generate output folder in BIDS structure
    mkdir([main_dir '\derivatives\' extractBefore(current_datasets_savename, '_ses') '\'  'ses-' curreent_session{1} '\eeg\SOUND_SSP-SIR\' current_datasets_savename]); %folder where to save preprocessing intermediates;
    
    %save path for intermediates and figures
    current_output_folder= [main_dir '\derivatives\' extractBefore(current_datasets_savename, '_ses') '\'  'ses-' curreent_session{1} '\eeg\SOUND_SSP-SIR\' current_datasets_savename];
    
    %start saving workspace for SOUND
    diary([current_output_folder '\' current_datasets_savename '_workspace_ICA']);
    
    
    % load dataset (from previous shell)
    EEG = pop_loadset('filename',[current_datasets_savename '_trialrej.set'],'filepath',  current_output_folder);
    
    % ICA decomposition
    EEG = pop_runica(EEG, 'chanind', 1:length(EEG.chanlocs), 'interupt','on'); % to exclude APBs
    
    
    % save set
    EEG.setname= [current_datasets_savename '_After_ICA'];
    pop_saveset(EEG, 'filename',[current_datasets_savename '_After_ICA'],'filepath', current_output_folder);
    
    
    disp('ICA decomposition ended')
    diary off %stop saving workspace
    
end




%% 5. Remove ocular artifat with ICA

for maincounter=1:length(all_datasets_list) %loops through all datasets  %loops through all datasets
    
    %build save name of the dataset
    current_datasets = all_datasets_name{maincounter}; %current dataset name
    curreent_session= extractBetween(current_datasets, 15,16); %current dataset session (T0, T1... etc)
    current_datasets_savename= all_datasets_name{maincounter}; %name with the bids structure
    
    %generate output folder in BIDS structure
    mkdir([main_dir '\derivatives\' extractBefore(current_datasets_savename, '_ses') '\'  'ses-' curreent_session{1} '\eeg\SOUND_SSP-SIR\' current_datasets_savename]); %folder where to save preprocessing intermediates;
    
    %save path for intermediates and figures
    current_output_folder= [main_dir '\derivatives\' extractBefore(current_datasets_savename, '_ses') '\'  'ses-' curreent_session{1} '\eeg\SOUND_SSP-SIR\' current_datasets_savename];
    
    %start saving workspace for SOUND
    diary([current_output_folder '\' current_datasets_savename '_workspace_ICAremoval']);
    
    EEG = pop_loadset('filename',[current_datasets_savename '_After_ICA.set'],'filepath',  current_output_folder);
    eeglab redraw % otherwise sometimes it makes error during ICA selection
    close all
    
    
    %to select ocular IC
    EEG = pop_selectcomps(EEG, 1:size(EEG.icaweights,1) );%plot interactive topography
    
    waitfor( findobj('parent', figure(2), 'string', 'OK'), 'userdata'); %to pause the script until ocular IC have been marked
    close all% close all the windows after you pressed OK in one figure
    
    %  save bad ICA comp on HD
    EEG.BadIC=find(EEG.reject.gcompreject==1);
    bad_ICAcomp=EEG.BadIC;
    save([current_output_folder '\bad_ICAcomp_' current_datasets_savename], 'bad_ICAcomp');
    
    %to remove marked ICA components
    if isempty(EEG.BadIC)==0 %if EEG.BadIC is not empty
        EEG = pop_subcomp( EEG, EEG.BadIC, 0); %remove the selected components
    end
    
    % visualize TEP
    EEG_epoch=EEG;
    EEG_epoch= pop_rmbase(EEG_epoch, [-750 0]);
    TEP_dummy(:,:)=mean(EEG_epoch.data, 3); %average across trials.
    figure
    plot(EEG_epoch.times, TEP_dummy(:,:)); %plot avg trial TEP
    hold on
    xlim([-100 400]);%restrict epoch
    ylim([-50 50]); % restrict amplitude range
    title('After ICAremoval');
    hold off
    saveas(gcf,...  %save .fig
        [current_output_folder '\' current_datasets_savename '_ICAremoval'])
    saveas(gcf,... %save .png
        [current_output_folder '\' current_datasets_savename '_ICAremoval'],...
        'png')
    close all
    
    
    %  save data after ICA rejection:
    EEG.setname= [current_datasets_savename '_ICAremoval'];
    pop_saveset(EEG, 'filename',[current_datasets_savename '_ICAremoval'],'filepath', current_output_folder);
    
    
    disp('ICArej ended')
    diary off %stop saving workspace
    
end



%% 6. Use SSP-SIR to clean muscle artifacts



for maincounter=1:length(all_datasets_list) %loops through all datasets  %loops through all datasets
    
    %build save name of the dataset
    current_datasets = all_datasets_name{maincounter}; %current dataset name
    curreent_session= extractBetween(current_datasets, 15,16); %current dataset session (T0, T1... etc)
    current_datasets_savename= all_datasets_name{maincounter}; %name with the bids structure
    
    %generate output folder in BIDS structure
    mkdir([main_dir '\derivatives\' extractBefore(current_datasets_savename, '_ses') '\'  'ses-' curreent_session{1} '\eeg\SOUND_SSP-SIR\' current_datasets_savename]); %folder where to save preprocessing intermediates;
    
    %save path for intermediates and figures
    current_output_folder= [main_dir '\derivatives\' extractBefore(current_datasets_savename, '_ses') '\'  'ses-' curreent_session{1} '\eeg\SOUND_SSP-SIR\' current_datasets_savename];
    
    %start saving workspace for SOUND
    diary([current_output_folder '\' current_datasets_savename '_workspace_SSP-SIR']);
    
    % load dataset (from previous shell)
    EEG = pop_loadset('filename',[current_datasets_savename '_ICAremoval.set'],'filepath',  current_output_folder);
    
    
    % replace the interp time with constant 0 values
    EEG = pop_tesa_removedata( EEG, interp_interval);
    
    % Use SSP-SIR to clean muscle artifacts
    
    load([current_output_folder '\LFM_sphere'],  'LFM_sphere'); %load lead field
    
    EEG_clean_muscle = EEG;
    EEG_clean_ica = EEG;
    
    dimN  =rank(mean(EEG_clean_ica.data(1:end,:,:),3));  %GIACOMO changed 1:61 in 1:end
    
    %to create TEPs with incremental number of removed components, up to 5
    for iii=1:max_SSP_SIR_PC_to_remove
        
        % artifact_topographies has been added to save this variable in the EEG
        % dataset
        % old version in which components to remove are selected manually
        %[~, artifact_topographies, ~, ~, filt_ker,SSP_SIR_operator] = SSP_SIR(mean(EEG_clean_ica.data(1:66,:,:),3), [], LFM_sphere,  dimN - PC, [], 0 ,EEG_clean_ica.srate,EEG.times,[0,50]); %[0, 50] is the time window on which SSP-SIR is applied
        [~, artifact_topographies, ~, ~, filt_ker,SSP_SIR_operator, PCs_variance] = SSP_SIR_GB(...
            mean(EEG_clean_ica.data(1:end,:,:),3),...  %GIACOMO changed 1:61 in 1:end
            iii,...
            LFM_sphere,...
            [],...
            [],...
            0,...
            EEG_clean_ica.srate,...
            EEG.times,...
            ssp_sir_timerange); % is the time window in ms on which SSP-SIR is applied. it should start at least at the begninning of the interpolated time. If it starts before, it will mess up a bit the baseline, if it starts after it might not caputure the component to remove.
        
        list_PCs_variance(maincounter, iii)=PCs_variance; %store variance of SSP_SIR PCs
        
        % to save figures of PC created by SSP-SIR
        saveas(gcf,[current_output_folder '\' current_datasets_savename '_PC' num2str(iii-1)],'fig')
        saveas(gcf,[current_output_folder '\' current_datasets_savename  '_PC' num2str(iii-1)],'jpg')
        
        
        close all
        
        for i = 1:size(EEG_clean_ica.data,3)
            EEG_clean_muscle.data(1:end,:,i) = filt_ker.*(SSP_SIR_operator*EEG_clean_ica.data(1:end,:,i)) +....
                EEG_clean_ica.data(1:end,:,i) - filt_ker.*EEG_clean_ica.data(1:end,:,i);
        end
        
        
        % to save topographies of removed components. This allows to count how many
        % topographies have been removed
        EEG_clean_muscle.SSPSIR_topographies = artifact_topographies;
        
        %save dataset with -N PC
        EEG = EEG_clean_muscle;
        EEG.setname= 'After SSP-SIR';
        EEG = pop_saveset( EEG, 'filename',[current_datasets_savename '_After_SSP-SIR_PC_' num2str(iii-1) '.set'],...
            'filepath', current_output_folder);
        
        % visualize TEP
        EEG_epoch=EEG;
        TEP_dummy(:,:)=mean(EEG_epoch.data, 3); %average across trials. Save all subjects in order to plot the grand-average artifact in the end of the shell
        figure
        plot(EEG_epoch.times, TEP_dummy(:,:)); %plot avg trial TEP
        hold on
        xlim([-100 400]);%restrict epoch
        ylim([-50 50]); % restrict amplitude range
        title(['After SSP-SIR check effect on the whole epoch PC removed: ' num2str(iii-1)]);
        
        hold off
        %save figure
        saveas(gcf,...  %save .fig
            [current_output_folder '\' current_datasets_savename '_afterSSP-SIR_PCs_' num2str(iii-1)])
        saveas(gcf,... %save .png
            [current_output_folder '\' current_datasets_savename '_afterSSP-SIR_PCs_' num2str(iii-1)],...
            'png')
        
        close all
        
        
        % prepare dataset to topoplot from GFP
        EEG_vis=pop_reref(EEG,[]);
        EEG_vis=pop_rmbase(EEG_vis, [-100 -2]);
        
        TEPs=mean(EEG_vis.data,3); %make TEPs (need for plot)
        
        % calculate Global Field Power
        gfp = std(mean(EEG_vis.data,3)).^2; %^2 becasue std does a sqrt on variance
        gfp_baseline.mean= mean(gfp(:,4500:4995),2);
        gfp_baseline.std= std(gfp(:,4500:4995));
        
        % find peaks
        [PKS,LOCS]=findpeaks(gfp,... %all gfp data
            5000,...%SR
            'MinPeakHeight', gfp_baseline.mean+(3*gfp_baseline.std),... %minimum peak height, should be 3 SD from baseline (mean + 3* SD baseline) %do it on GFP
            'MinPeakProminence', 0.75);
        LOCS=LOCS*1000; %transform latencies in ms from the TMS pulse
        LOCS=LOCS-1000; % move 0 to point 5000 (the epochs are of 10 000 point, 0 is point 5000. SR= 5000Hz)
        peaks_after_pulse=find(LOCS>-5 & LOCS<50);%take latencies and aplitude of peaks between -5 and 50 ms (it starts 4 ms before the interp interval -1 to 6, to be sure it will caputre peaks near 0)
        LOCS=LOCS(peaks_after_pulse);
        PKS=PKS(peaks_after_pulse);
        
        if length(LOCS)==1
            final_LOCS=[LOCS LOCS]; %this beacuse one topoplot alone produced an error (bad fixing)
            final_PKS=[PKS PKS]; %take also peaks value for marker in the next plot
        elseif isempty(LOCS)
            final_LOCS=[59 59]; %this becasue no peaks produced an error (bad fixing)
            final_PKS=[gfp(1059*5) gfp(1059*5)]; % %take also peaks value at 59 ms for marker in the next plot
        else
            final_LOCS=LOCS;
            final_PKS=PKS;
        end
        
        % plot topoplots at peaks after interpolation
        pop_topoplot(EEG_vis, 1000, final_LOCS , 'After SSP-SIR', [4 length(final_LOCS)], 0,'electrodes','on'); %plots topoplots in one line
        title('\muV');  %visualize topoplot in average reference
        
        hold on
        
        %make window larger
        H=gcf;
        set(H,'units', 'normalized','outerposition',[0 0 1 1]);
        
        %Plot GFP
        subplot(3,1,2) %add a subplot on the bottom
        plot(EEG_vis.times(4950:5300), gfp(4950:5300), 'r', 'linewidth', 3); % plots GFP between -10 +60 ms
        hold on;
        title(['GFP, PC: ' num2str(iii-1)]);
        xlabel('time (ms)');
        ylabel('\muV^2')
        xticks(-20:1:60);
        %set(gca,'box','off') %remove upper xticks and right yticks
        set(gca,'Units', 'centimeters','TickLength',[0.004, 0.004]) %set length x y ticks
        
        hold on
        %insert markers at peak position
        plot(final_LOCS, final_PKS, 'b*', 'LineWidth',2)
        
        %butterfly plot TEPs
        subplot(3,1,3) %add a subplot on the bottom
        for x=1:length(EEG.chanlocs)
            plot(EEG_vis.times(4950:5300), TEPs(x,(4950:5300))) %plot -10 + 60
            hold on
        end
        title(['TEPs -10 +60, PC: ' num2str(iii-1)]);
        xlabel('time (ms)');
        ylabel('\muV')
        xticks(-20:1:60);
        %set(gca,'box','off') %remove upper xticks and right yticks
        set(gca,'Units', 'centimeters','TickLength',[0.004, 0.004]) %set length x y ticks
        
        hold off
        
        % save figure
        saveas(gcf,[current_output_folder '\' current_datasets_savename '_merge_' num2str(iii-1)],'fig')
        saveas(gcf,[current_output_folder '\' current_datasets_savename '_merge_' num2str(iii-1)],'jpg')
        
        
        close
        
    end
    
    
    diary off %stop saving workspace
    
    disp('SSP-SIR ended');
    
end

%USE IN case you want an automatic SSP-SIR PC rejection
% construct the list of PCs to remove, based on the variance over channel.
% Remove component till the varinace drops below (stict) 35. Useful for
% automatic SSP-SIR PCs rejection
% % % suggested_PCs_to_remove=[];
% % % for i=1:length(all_datasets_list)
% % %     suggested_PCs_to_remove(i,1)= find(list_PCs_variance(i,:)<35, 1)-1; %suggested number of PCs to remove
% % % end

%% 7. Final steps of preprocessing after SSP-SIR


for maincounter=1:length(all_datasets_list) %loops through all datasets  %loops through all datasets
    
    %build save name of the dataset
    current_datasets = all_datasets_name{maincounter}; %current dataset name
    curreent_session= extractBetween(current_datasets, 15,16); %current dataset session (T0, T1... etc)
    current_datasets_savename= all_datasets_name{maincounter}; %name with the bids structure
    
    %generate output folder in BIDS structure
    mkdir([main_dir '\derivatives\' extractBefore(current_datasets_savename, '_ses') '\'  'ses-' curreent_session{1} '\eeg\SOUND_SSP-SIR\' current_datasets_savename]); %folder where to save preprocessing intermediates;
    
    %save path for intermediates and figures
    current_output_folder= [main_dir '\derivatives\' extractBefore(current_datasets_savename, '_ses') '\'  'ses-' curreent_session{1} '\eeg\SOUND_SSP-SIR\' current_datasets_savename];
    
    %start saving workspace for SOUND
    diary([current_output_folder '\' current_datasets_savename '_workspace_Final_step']);
    
    
    %If you want to choose manually how many SSP_SIR PCs to
    %remove, you can use something like this
    num_PCs_to_remove=input(['How many components to remove for ' current_datasets_savename ' ? ']);
    EEG = pop_loadset('filename',[current_datasets_savename '_After_SSP-SIR_PC_' num2str(num_PCs_to_remove) '.set'],...
        'filepath', current_output_folder);
    
    
    %If you already made a list you can use something like this instead (automatic, using the suggested PCs)
    %Load the dataset with SSP_SIR PCs removed: (already chosen)
    % %    EEG = pop_loadset('filename',[current_datasets_savename '_After_SSP-SIR_PC_' num2str(suggested_PCs_to_remove{maincounter,2}) '.set'],...
    % %                       'filepath', current_output_folder);
    
    %INTERPOLATION to remove the residual tms-pulse artifact left between
    %-1 and 6ms
    EEG = pop_tesa_removedata( EEG, interp_interval);
    
    % Interpolate missing data around TMS pulse
    % interpolation of missing data aroun TMS pulse using a cubic
    EEG = pop_tesa_interpdata( EEG, 'cubic', [1,1] );
    
    % low pass filter 90 Hz
    EEG = pop_eegfiltnew(EEG, 'hicutoff', low_pass_filt ,'plotfreqz',1);
    
    % Notch at 50 Hz
    EEG = pop_tesa_filtbutter( EEG, notch_filt(1), notch_filt(2), 4, 'bandstop' );
    close  all
    
    % reref to the average
    EEG = pop_reref( EEG, final_ref);
    
    %  DownSample data (5000 Hz to 1000 Hz).
    EEG = pop_resample( EEG, downsample);
    
    % last baseline correction
    EEG = pop_rmbase(EEG, baseline_long);
    
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
    % save final set short
    EEG.setname= [current_datasets_savename '_Final_scaled_avg_ref'];
    pop_saveset(EEG, 'filename',[current_datasets_savename '_Final_scaled_avg_ref'],'filepath', current_output_folder);
   
    
    disp('Preprocessing ended')
    diary off %stop saving workspace
end

disp('All preprocessing ended')
