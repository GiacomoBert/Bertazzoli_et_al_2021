% AN EXAMPLE MATLAB SCRIPT FOR USING DDWIENER AND SOUND
%
% This is an EEGLAB-compatible MATLAB script that demonstrates how DDWiener
% and SOUND can be used to clean a partially noisy EEG dataset.
%
% Required toolbox to run this script: BioSig, TESA, bva-io (plugin to import BrainVision data)
% 18 September 2017 : Tuomas Mutanen and Marta Bortoletto.
% Department of Neuroscience and Biomedical Engineering (NBE), School of Science, Aalto University
% .........................................................................


% Relevant Variables that are created in the script
% EEG_orig %original continous data
% EEG_janela %continous after janelacut
% EEG_artrej %epoched after rejection (148 epochs)
% EEG_clean %after sound
% EEG_sICA %after ICA (performed after sound)
% EEG_clean_muscle %after SSP-SIR

% .........................................................................


%% Define directories and subjects
% STEP 1: Data pre-processing up to artifact rejection

clear;		% prepare a clean workspace with no variables
close all;		% close all open events


% define directories
% all directories must have been created
% addpath C:\Users\Marta\Documents\fieldtrip-20171129
% locpath=('C:\Users\Marta\Documents\eeglab14_1_1b\plugins\dipfit2.3\standard_BESA\standard-10-5-cap385.elp'); %Replace first part with your location of eeglab % locpath = your location of 'standard-10-5-cap385_NoEOG.elp' file (located in eeglab directory)
datalocation =('C:\Users\Neurofisilogia\Documents\MATLAB\Giacomo\GR2016_analyses\all_raw_data'); % datalocation =  write the directory storing your raw data
stepfolder='C:\Users\Neurofisilogia\Documents\MATLAB\Giacomo\GR2016_analyses\Analyses';

figurefolder='C:\Users\Neurofisilogia\Documents\MATLAB\Giacomo\GR2016_analyses\Analyses\Figures';
outputfolder='C:\Users\Neurofisilogia\Documents\MATLAB\Giacomo\GR2016_analyses\Analyses\Final'; % outputfolder = write the directory containing your processed data
SSP_SIRfolder='C:\Users\Neurofisilogia\Documents\MATLAB\Giacomo\GR2016_analyses\Analyses\SSP_SIR';


% find all .vhdr file
ALLsubjects= [dir([datalocation '\*LP.vhdr'])',...
    dir([datalocation '\*RP.vhdr'])',...
    dir([datalocation '\*LF.vhdr'])',...
    dir([datalocation '\*RF.vhdr'])',...
    dir([datalocation '\*SH.vhdr'])'];

%definire elenco soggetti INSERIRE QUI TUTTI DATASET
% INSERIRE IF LOOP PER REREF: QUANDO STIMOLO DX REREF TP9 E VICEVERSA
subjects    =  strtok({ALLsubjects.name}, '.');

%Event selection for last steps from 6
trigger= 132;

 subjects=natsort(subjects); %put datasets in alphabetic order
%% Start loop for processing
% type the subject or a list of subjects to be analyzed (in this case subjects will be sequentially processed)
for isubj=1:numel(subjects)
    subj = subjects{isubj};
    fprintf('Processing: %s \n', [subj '.vhdr']);
    
   %{
    %% 1. Read EEG data from Biosemi, read channel locations and save .set
    %  Use EEGLAB with BioSig toolbox to read in the EEG dataset and the corresponding electrode locations.
    %  Depending on your computer, this may take a moment.
    
    disp('Shell 1 started')
    
    % Start EEGLAB under MATLAB
    eeglab;
    close
    
    % Read in the EEG dataset
    EEG = pop_loadbv(datalocation, [subj '.vhdr']);
    
    % Read in the channel location file and get channel location
    EEG = eeg_checkset(EEG);
    EEG = pop_chanedit(EEG,  'lookup', 'standard-10-5-cap385.elp'); %adds channel locations from standard cap in locpath
    [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    eeglab redraw
    close
    
    % define first 67 electrodes as EEG and others as EOG (VEOG, HEOG, APBsx,
    % APBdx)
    %reference on the nose
    for i = 1:68
        if i < 68
            EEG.chanlocs(i).type = 'EEG';
        else
            EEG.chanlocs(i).type = 'EOG';
            EEG.chanlocs(i).theta = [];
            EEG.chanlocs(i).radius = [];
            EEG.chanlocs(i).X = [];
            EEG.chanlocs(i).Y = [];
            EEG.chanlocs(i).Z = [];
            EEG.chanlocs(i).sph_theta = [];
            EEG.chanlocs(i).sph_phi = [];
            EEG.chanlocs(i).sph_radius = [];
            EEG.chanlocs(i).urchan = [];
        end
    end
    eeglab redraw
    close
    
    % Create lead-field matrix which allows to convert signal between sensors and sources
    [LFM_sphere] = construct_spherical_lead_field(EEG);
    
    %save EEG_orig, dataset and figure
    EEG_orig=EEG;
    EEG = pop_saveset( EEG, 'filename',[subj '_1original.set'],'filepath',stepfolder);
    [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    eeglab redraw
    close
    
    pop_eegplot( EEG, 1, 1, 1);
    saveas(gcf,[figurefolder '\' subj '_original'],'jpg')
    close(gcf)
    
    
    
    % to control the effect of janelacut - step 1
    EEG = eeg_checkset( EEG );
    EEG = pop_epoch( EEG, {  'S132' }, [-0.01        0.03], 'epochinfo', 'yes');
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname','Before_janelacut control','gui','off');
    EEG = eeg_checkset( EEG );
    EEG = pop_rmbase( EEG, [-10  -2]);
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'gui','off');
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset( EEG, 'filename',[subj '_Before_janela.set'],'filepath',stepfolder);
    [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG_Before_janelacontrol=EEG;
    disp('Shell 1 ended')
    
    %% 2. Pre-process the EEG signals
    % reload  original set
    disp('Shell 2 started')
    STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
    EEG = pop_loadset('filename',[subj '_1original.set'],'filepath', stepfolder);
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    eeglab redraw
    close
    
    % % This command defines TMS triggers and eliminates non TMS triggers
    % % (useless, to be removed in future scripts)
    triggers=[];
    triggers_to_erase=[];
    for jj=1:length(EEG.urevent)
        if strcmp(EEG.urevent(jj).type,'S132') %'S150'  'S151'  'S152'  'S160'  'S161'  'S162'
            triggers=cat(2,triggers,EEG.urevent(jj).latency);
        else
            triggers_to_erase=cat(1,triggers_to_erase,jj);
        end
    end
    
    EEG.urevent(triggers_to_erase)=[];%% erase unwanted trigger like boundery
    EEG.event(triggers_to_erase)=[];
    eeglab redraw
    close
    
    % Interpolate time range (Silvia Casarotto method)
    janelacut=[-1, 6]; % define the time window (ms) to be interpolated
    % I have changed this parameter to set it as Silvia---------------
    swindow=[-4, 2];
    %---------------------------------------------------------------
    metodo='moving';
    n=5;
    EEG.data=artefatto_BrainAMP(EEG,triggers,janelacut,swindow,metodo,n);% the file 'artefatto_BrainAMP.m' needs to be in the directory
    clear janelacut swindow metodo n triggerstmp
    eeglab redraw
    close
    EEG_janela=EEG;
    
    
    %to control the effect of janelacut
    EEG = eeg_checkset( EEG );
    EEG = pop_epoch( EEG, {  'S132'  }, [-0.01        0.03], 'epochinfo', 'yes');
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname','After_janelacut control','gui','off');
    EEG = eeg_checkset( EEG );
    EEG = pop_rmbase( EEG, [-10  -2]);
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'gui','off');
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset( EEG, 'filename',[subj '_After_janela.set'],'filepath',stepfolder);
    [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG_After_janelacontrol=EEG;
    disp('Shell 1 ended')
    
    %reload after_janela and before_janela
    STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
    EEG=EEG_After_janelacontrol;
    eeglab redraw
    close
    ALLEEG(2)= EEG_Before_janelacontrol;
    eeglab redraw
    close
    
    % save figure with after-before comparison
    pop_comperp( ALLEEG, 1, 1,2,'addavg','on','addstd','off','subavg','on','diffavg','on','diffstd','off','tplotopt',{'ydir' -1});
    saveas(gcf,[figurefolder '\' strtok(subj,'.') '_janelacontrol'],'fig')
    close(gcf)
    
    %reload data after janela
    ALLEEG = pop_delset( ALLEEG, [2] );
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'retrieve',1,'study',0);
    eeglab redraw
    close
    EEG=EEG_janela;
    eeglab redraw
    close
    
    % High-pass filter (0.1 Hz)
    % We may compare results with different type of filters
    [EEG,~,b] = pop_firws(EEG, 'fcutoff', 0.1, 'ftype', 'highpass', 'wtype', 'hamming', 'forder', 4224, 'minphase', 0);
    
    % Epoch data (-200 to 500 ms)
    EEG = pop_epoch(EEG, {'S132' }, [-0.5 0.5]);
    eeglab redraw
    close
    
    % DownSample data (5000 Hz to 2048 Hz)
    EEG = pop_resample( EEG, 2048);
    
    % Baseline correction
    EEG = pop_rmbase(EEG, [-100 -2]);
    eeglab redraw
    close
    
    %RE-REFERENCE
    %------------------------------------------------------- iNSERIRE if DX E
    %sx
    % Re-reference cortical electrodes to TP10 is timulation is left, reref
    % to TP9 if stimulation is right or sham
    if endsWith(subj, 'LP') == 1 || endsWith(subj, 'LF') == 1
        EEG = pop_reref(EEG, 47); %47 is TP10
    elseif endsWith(subj, 'RP') == 1 || endsWith(subj, 'RF') == 1 || endsWith(subj, 'SH') == 1
        EEG = pop_reref(EEG, 38); %38 is TP9
    end
    
    % Adjust lead-field matrix on the new reference
    LFM_sphere_ref= repmat(LFM_sphere(47,:), size(LFM_sphere, 1)-1, 1);
    LFM_sphere = LFM_sphere([1:46 48:67],:) - LFM_sphere_ref;
    %-----------------------------------------------------------------
    ALLEEG(1) = EEG;
    eeglab redraw
    close
    
    % save dataset
    EEG = pop_saveset( EEG, 'filename',[subj '_2PreProc.set'],'filepath',stepfolder);
    save ([stepfolder '\' subj '_LFM_sphere.txt'], 'LFM_sphere', '-ascii');
    clear eval
    
    % figure
    figure; pop_plottopo(EEG, [1:66] , ' resampled', 0, 'ydir',1); % this figure shows averaged TEP at this stage of processing (interactive plot: click on the channels to zoom)
    saveas(gcf,[figurefolder '\' subj '_TEPs_original'],'fig')
    close(gcf)
    disp('Shell 2 ended')
    
    %% 3. Artifact rejection - This artefact rejection steps have not been used - they should be removed together with the next save dataset
    
    disp('Shell 2 started')
    STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
    EEG = pop_loadset('filename',[subj '_2PreProc.set'],'filepath', stepfolder);
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    
    
    
    % % Artifact rejection to manually reject epochs
    %  % (a window pops up) scroll through epochs and select trials to be rejected by
    %  % clicking on the bad epoch (they turn into yellow). At the end press
    %  % 'reject' button.
    %
    %  pop_eegplot(EEG, 1, 1, 1);
    %
    %  waitfor( findobj('parent', gcf, 'string', 'REJECT'), 'userdata');
    %
    % % changes are saved in dataset1, dataset2 is removed, dataset1 is
    % % retrieved to work on it
    % ALLEEG(1)=EEG;
    % ALLEEG = pop_delset( ALLEEG, [2] );
    % [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'retrieve',1,'study',0);
    
    % eeglab redraw
    
    % save EEG_artrej, dataset and figure
    EEG = pop_saveset( EEG, 'filename',[subj '_3ArtRej.set'],'filepath',stepfolder);
    % EEG_artrej = EEG;
    % EEG_artrej.setname='Artrej';
    % EEG_artrej=EEG;
    %save EEG_orig, dataset and figure
    
    
    
    
    %% 4. Use SOUND to clean the channel-specific noise
    disp('Shell 2 started')
    STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
    EEG = pop_loadset('filename',[subj '_3ArtRej.set'],'filepath', stepfolder);
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    
    LFM_sphere= load([stepfolder '\' subj '_LFM_sphere.txt']);
    
    tic
    disp('Shell 4 started')
    
    % Build the spherical lead field
    %[LFM_sphere] = construct_spherical_lead_field(EEG);
    %
    % % % Re-reference the data and the lead field to the channel with the least noise
    % [~, sigmas] = DDWiener(EEG_evo);
    % [~,bestC] = min(sigmas);
    % [datatmp] = ref_best(EEG_evo, bestC);
    % [LFM_sphere_tmp] = ref_best(LFM_sphere, bestC);
    
    
    % Set the SOUND parameters:
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
    %EEG_clean.setname='after sound';
    
    cat_tri_N = 5;
    
    for i =1:cat_tri_N:size(EEG.data,2)
        
        N_tr= size(EEG.data,3);
        
        % Run the SOUND algorithm:
        i
        tmp_data = squeeze(EEG.data(1:66,i,:));
        cat_tri_real=1;
        for j = 1:(cat_tri_N-1)
            if i+j<= size(EEG.data,2)
                tmp_data = cat(2,tmp_data,squeeze(EEG.data(1:66,i+j,:)));
                cat_tri_real = cat_tri_real+1;
            end
        end
        [corrected_data] = SOUND(tmp_data, LFM_sphere, iter, lambda_value);
        
        for j = 0:(cat_tri_real-1)
            j
            EEG_clean.data(1:66,i+j,:) = reshape(corrected_data(:,(j*N_tr+1):(j+1)*N_tr),[size(corrected_data,1),1, N_tr]);
        end
        
    end
    % (The data contain several channels and trials so this may take several
    % minutes.)
    
    %  figure; pop_plottopo(EEG_clean, [1:20] , ' resampled', 0, 'ydir',1);
    
    ALLEEG(2) = EEG_clean;
    eeglab redraw
    
    % figure;
    %pop_comperp(ALLEEG);
    % % Re-reference the data and the lead field to the channel average:
    % [LFM_sphere_mean] = ref_ave(LFM_sphere);
    % ERP_cleaned = LFM_sphere_mean*x;
    
    %----------------------------------------------------------------
    %save
    EEG.setname= 'After SOUND';
    ALLEEG(2).setname= 'After SOUND';
    EEG = pop_saveset( EEG_clean, 'filename',[subj '_4Sound.set'],'filepath',stepfolder);
    %---------------------------------------------------------------------
    
    disp('Shell 4 ended')
    
    toc
    
    %>____________________________________________________________________________________________________--ARRIVARE
    %FINO A QUI
    
   
    
    
    %% 4b. Artifact rejection
    STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
    EEG = pop_loadset('filename',[subj '_3ArtRej.set'],'filepath', stepfolder);
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    EEG = pop_loadset('filename',[subj '_4Sound.set'],'filepath', stepfolder);
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'retrieve',2,'study',0);
    eeglab redraw
    
    
    % Artifact rejection to manually reject epochs
    % (a window pops up) scroll through epochs and select trials to be rejected by
    % clicking on the bad epoch (they turn into yellow). At the end press
    % 'reject' button.
    
    pop_eegplot(EEG, 1, 1, 1);
    
    waitfor( findobj('parent', gcf, 'string', 'REJECT'), 'userdata');
    
    % changes are saved in dataset2, dataset3 is removed, dataset2 is
    % retrieved to work on it
    ALLEEG(2)=EEG;
    ALLEEG = pop_delset( ALLEEG, [3] );
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'retrieve',2,'study',0);
    eeglab redraw
    
    % save EEG_artrej, dataset
    EEG = pop_saveset( EEG, 'filename',[subj '_4Sound_ar.set'],'filepath',stepfolder);
    EEG_artrej = EEG;
    EEG_artrej.setname='Artrej';
    
    
    
    
    %% 5. ICA
    STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
    EEG = pop_loadset('filename',[subj '_3ArtRej.set'],'filepath', stepfolder);
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    EEG = pop_loadset('filename',[subj '_4Sound_ar.set'],'filepath', stepfolder);
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'retrieve',2,'study',0);
    eeglab redraw
    
    
    %plot before ica butterfly
    % pop_tesa_plot( EEG, 'tepType', 'data', 'xlim', [-100 350], 'ylim', [-30 30], 'CI','off','plotPeak','off' );
    % title(['Before ICA']);
    % saveas(gcf,[figurefolder '\' subj '_BeforeICA'],'jpg')
    % close(gcf)
    
    % ICA on EEG
    EEG = pop_runica(EEG, 'chanind', [1:66], 'interupt','on'); % to exclude EOG that has not been corrected with sound!!!!!!!!!!!!!!
    eeglab redraw
    
    %save ICA all components
    EEG_ALLICA=EEG;
    EEG_ALLICA.setname= 'After soundICA';
    % ALLEEG(2).setname= 'After SoundICA';
    EEG = pop_saveset( EEG_ALLICA, 'filename',[subj '_5allICA.set'],'filepath',stepfolder);
     
    
 
    
    %% 5b. Remove BAD ICA
    
    STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
    EEG = pop_loadset('filename',[subj '_5allICA.set'],'filepath', stepfolder);
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'retrieve',1,'study',0);
    eeglab redraw
    
    %to select ocular IC
    EEG = pop_selectcomps(EEG, [1:size(EEG.icaweights,1)] );%plot interactive topography
    %pop_plotdata(EEG, 0, [1:size(EEG.icasphere,1)] , [1:EEG.trials], 'Averaged ICA', 0, 1, [ -3 3] )% plot avarage
    % ordining on the basis of kurtosis
    
    waitfor( findobj('parent', gcf, 'string', 'OK'), 'userdata'); %to pause the script until ocular IC have been marked
    
    %to remove marked ICA components
    EEG.BadIC=find(EEG.reject.gcompreject==1);
    if isempty(EEG.BadIC)==0
        EEG = pop_subcomp( EEG, EEG.BadIC, 0);
    end
    
    % %plot after ica butterfly
    % pop_tesa_plot( EEG, 'tepType', 'data', 'xlim', [-100 350], 'ylim', [-30 30], 'CI','off','plotPeak','off' );
    % title(['After ICA']);
    % saveas(gcf,[figurefolder '\' subj '_AfterICA'],'jpg')
    % close(gcf)
    
    %save
    EEG.setname= 'After soundICA';
    ALLEEG.setname= 'After SoundICA';
    EEG = pop_saveset( EEG, 'filename',[subj '_5ICA.set'],'filepath',stepfolder);
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'retrieve',1,'study',0);
    eeglab redraw
    
    EEG_sICA = EEG;
    EEG_sICA.setname='after soundICA';
    %}
    
    %% 6. Use SSP-SIR to clean muscle artifacts
    % 150 premotoria
    % 151 parietale
    % 152 motoria
    
    
    %to repeat last steps for each trigger
    for t=1: length(trigger)
        
        % to retrieve data after sound and continue with SSP-SIR
        STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
        % EEG = pop_loadset('filename',[subj '_4Sound_ar.set'],'filepath', stepfolder);
        % [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
        EEG = pop_loadset('filename',[subj '_5ICA.set'],'filepath', stepfolder);
        [ALLEEG, EEG, ~] = eeg_store( ALLEEG, EEG, 0 );
        % EEG = pop_loadset('filename',[subj '_TEPs_SOUNDonly.set'],'filepath', stepfolder);
        % [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
        % [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'retrieve',2,'study',0);
        eeglab redraw
        
        % event selection with variable "trigger", defined at the beginning of the
        % script
        EEG = eeg_checkset( EEG );
        EEG = pop_selectevent( EEG, 'type',{[ 'S' num2str(trigger(t))]},'deleteevents','on','deleteepochs','on','invertepochs','off');
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',num2str(trigger(t)),'gui','off');
        eeglab redraw
        
        
        % Use SSP-SIR to clean muscle artifacts
        % re-load  ICA-sphere
        LFM_sphere= load([stepfolder '\' subj '_LFM_sphere.txt']);
        
        % SSP-SIR
        % when the window pops up, select how many PCA components you want to remove
        % by pressing the spacebar. Each time you press the spacebar, one component
        % is removed.
        % When you have selected all the components you want to remove, click on
        % the side of the window to finish. The window closes.
        
        PC = 1; % Number of artifact dimensions
        
        EEG_clean_muscle = EEG;
        EEG = ALLEEG(2); % the second set is the one with selected trial for one condition
        EEG_clean_ica = EEG;
        
        dimN  =rank(mean(EEG_clean_ica.data(1:66,:,:),3));
        
        %to create TEPs with incremental number of removed components, up to 5
        for iii=1:6
            
            % artifact_topographies has been added to save this variable in the EEG
            % dataset
            % old version in which components to remove are selected manually
            %[~, artifact_topographies, ~, ~, filt_ker,SSP_SIR_operator] = SSP_SIR(mean(EEG_clean_ica.data(1:66,:,:),3), [], LFM_sphere,  dimN - PC, [], 0 ,EEG_clean_ica.srate,EEG.times,[0,50]); %[0, 50] is the time window on which SSP-SIR is applied
            [~, artifact_topographies, ~, ~, filt_ker,SSP_SIR_operator] = SSP_SIR_FISM(mean(EEG_clean_ica.data(1:66,:,:),3), iii, LFM_sphere,  dimN - PC, [], 0 ,EEG_clean_ica.srate,EEG.times,[0,50]); %[0, 50] is the time window on which SSP-SIR is applied
            
            
            
            % to save figures of PC created by SSP-SIR
            saveas(gcf,[SSP_SIRfolder '\' subj '_' num2str(trigger(t)) '_PC' num2str(iii-1)],'fig')
            close(gcf);
            
            
            
            
            for i = 1:size(EEG_clean_ica.data,3)
                i;
                EEG_clean_muscle.data(1:66,:,i) = filt_ker.*(SSP_SIR_operator*EEG_clean_ica.data(1:66,:,i)) +....
                    EEG_clean_ica.data(1:66,:,i) - filt_ker.*EEG_clean_ica.data(1:66,:,i);
            end
            ALLEEG(3) = EEG_clean_muscle;
            eeglab redraw
            
            % to save topographies of removed components. This allows to count how many
            % topographies have been removed
            EEG_clean_muscle.SSPSIR_topographies = artifact_topographies;
            
            
            figure;
            plot(filt_ker(1,:))
            saveas(gcf,[SSP_SIRfolder '\' subj '_filt_ker' num2str(trigger(t)) '_' num2str(iii-1)],'fig')
            close(gcf)
            
            %----------------------------------------------------------------
            save
            EEG.setname= 'After SSP-SIR';
            ALLEEG(5).setname= 'After SSP-SIR';
            EEG = pop_saveset( EEG_clean_muscle, 'filename',[subj '_6SSP-SIR_' num2str(trigger(t)) '_' num2str(iii-1) '.set'],'filepath',stepfolder);
            %---------------------------------------------------------------------
     
          
            %% 7. Last steps of preprocessing after SSP-SIR
            STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
            EEG = pop_loadset('filename',[subj '_4Sound.set'],'filepath', stepfolder);
            [ALLEEG, ~, ~] = eeg_store( ALLEEG, EEG, 0 );
            EEG = pop_loadset('filename',[subj '_5ICA.set'],'filepath', stepfolder);
            [ALLEEG, ~, ~] = eeg_store( ALLEEG, EEG, 0 );
            EEG = pop_loadset('filename',[subj '_6SSP-SIR_' num2str(trigger(t)) '_' num2str(iii-1) '.set'],'filepath', stepfolder);
            [ALLEEG, EEG, ~] = eeg_store( ALLEEG, EEG, 0 );
            [ALLEEG, ~, ~] = pop_newset(ALLEEG, EEG, 1,'retrieve',3,'study',0);
            eeglab redraw
            
            numset= 3; % 2 if you want to use data after sound, 3 if you want to use data after SSP-SIR
            EEG=ALLEEG(numset);
            ALLEEG(numset+1)= ALLEEG(numset); %to create a new dataset, on which run the preprocessing
            [ALLEEG, EEG, ~] = pop_newset(ALLEEG, EEG, 1,'retrieve',numset+1,'study',0);
            eeglab redraw
            
            % low pass filter 70 Hz
            EEG  = pop_basicfilter( EEG,  1:66 , 'Cutoff',  70, 'Design', 'butter', 'Filter', 'lowpass', 'Order',  4 ); % GUI: 02-Mar-2018 11:28:25
            [ALLEEG, EEG, ~] = pop_newset(ALLEEG, EEG, (numset+2),'setname','xx LP 70 Hz','gui','off');
            EEG = eeg_checkset( EEG );
            %EEG_sound_LP70=EEG; %to visualize filters only
            
            % Notch at 50 Hz
            % EEG  = pop_basicfilter( EEG,  1:72 , 'Cutoff',  50, 'Design', 'notch', 'Filter', 'PMnotch', 'Order',  180 ); % GUI: 02-Mar-2018 11:32:49
            % [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, (numset+3),'setname','xx Notch 50Hz','gui','off');
            % EEG = eeg_checkset( EEG );
            EEG_sound_filters=EEG; %to visualize filters only
            EEG_sound_filters.icachansind=[];
            
            %to check effects of filters
            pop_comperp( ALLEEG, 1, 4,5,'addavg','on','addstd','off','subavg','on','diffavg','on','diffstd','off','tplotopt',{'ydir' -1});
            % waitfor(gcf);
            % pop_comperp( ALLEEG, 1, 4,5,'addavg','on','addstd','off','subavg','on','diffavg','on','diffstd','off','tplotopt',{'ydir' -1});
            % waitfor(gcf);
            % EEG = eeg_checkset( EEG );
            % pop_comperp( ALLEEG, 1, 3,5,'addavg','on','addstd','off','subavg','on','diffavg','on','diffstd','off','tplotopt',{'ydir' -1});
            % waitfor(gcf);
            savefig([figurefolder '\' subj '_filters2_' num2str(trigger(t)) '.fig'])
            close(gcf)
            
            % artefact rejection 2
            %  pop_eegplot(EEG, 1, 1, 1);
            %  waitfor( findobj('parent', gcf, 'string', 'REJECT'), 'userdata');
            
            % changes are saved in dataset(numeset+1), dataset(numset+2) is removed, datasetnumset+1) is
            % retrieved to work on it
            ALLEEG(numset+1)=EEG;
            [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'retrieve', numset+1, 'study',0);
            ALLEEG = pop_delset( ALLEEG, numset+2:size(ALLEEG, 2) );
            [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',numset+1,'study',0);
            
            % rereference to linked mastoids
            % EEG=pop_chanedit(EEG, 'append',70,'changefield',{71 'labels' 'TP10'},'load',[],'lookup', locpath);
            % [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
            
            % ----------------------------------------- INSERIRE IF SX E DX
            % rereference to linked mastoids
            % EEG=pop_chanedit(EEG, 'append',70,'changefield',{71 'labels' 'TP10'},'load',[],'lookup', locpath);
            % [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
            
            % EEG=pop_chanedit(EEG, 'setref',{'1:66 68' 'TP10'},'changefield',{68 'type' 'EEG'});
            % [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
            
            % EEG = eeg_checkset( EEG );
            % EEG = pop_reref( EEG, [],'refloc',struct('labels',{'TP10'},'type',{'EEG'},'theta',{108.393},'radius',{0.66489},'X',{-23.3016},'Y',{-70.0758},'Z',{-42.0882},'sph_theta',{-108.393},'sph_phi',{-29.68},'sph_radius',{85},'urchan',{73},'ref',{'TP10'},'datachan',{0}),'exclude',[67] );
            % [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3,'setname','xx LP 70 Hz avgref','gui','off');
            
            % EEG = eeg_checkset( EEG );
            % EEG = pop_reref( EEG, [38 68] ,'exclude',[67] );
            % [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 4,'setname','xx LP 70 Hz mast','gui','off');
            % EEG = eeg_checkset( EEG );
            
            % ALLEEG(numset+1)=EEG;
            % [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve', numset+1, 'study',0);
            % ALLEEG = pop_delset( ALLEEG, [numset+2:size(ALLEEG, 2)]);
            % Re-reference cortical electrodes to TP10 is timulation is left, reref
            % to TP9 if stimulation is right or sham
            
            if endsWith(subj, 'LP') == 1 || endsWith(subj, 'LF') == 1
                EEG = pop_reref( EEG, [],...
                    'refloc',struct(...
                    'labels',{'TP10'},...
                    'type',{'EEG'},...
                    'theta',{108.393},...
                    'radius',{0.66489},...
                    'X',{-23.3016},...
                    'Y',{-70.0758},...
                    'Z',{-42.0882},...
                    'sph_theta',{-108.393},...
                    'sph_phi',{-29.68},...
                    'sph_radius',{85},...
                    'urchan',{47},...
                    'ref',{'TP10'},...
                    'datachan',{0},...
                    'sph_theta_besa',{119.68},...
                    'sph_phi_besa',{-18.393}),...
                    'exclude',67); %compute average reference, exclude VEOG, re-insert TP10
                
                %move TP10 to the correct position in chanlocs structure
                EEG.chanlocs(48:69)= EEG.chanlocs(47:68);
                EEG.chanlocs(47)= EEG.chanlocs(69);
                EEG.chanlocs(69)=[];
                
                %move TP10 to the correct position in data structure
                EEG.data(48:69,:,:)=EEG.data(47:68,:,:);
                EEG.data(47,:,:)=EEG.data(69,:,:);
                EEG.data(69,:,:)=[];
                
            elseif endsWith(subj, 'RP') == 1 || endsWith(subj, 'RF') == 1 || endsWith(subj, 'SH') == 1
                EEG = pop_reref( EEG, [],...
                    'refloc',struct(...
                    'labels',{'TP9'},...
                    'type',{'EEG'},...
                    'theta',{-108.393},...
                    'radius',{0.66489},...
                    'X',{-23.3016},...
                    'Y',{70.0758},...
                    'Z',{-42.0882},...
                    'sph_theta',{108.393},...
                    'sph_phi',{-29.68},...
                    'sph_radius',{85},...
                    'urchan',{38},...
                    'ref',{'TP9'},...
                    'datachan',{0},...
                    'sph_theta_besa',{-119.68},...
                    'sph_phi_besa',{18.393}),...
                    'exclude',67); %compute average reference, exclude VEOG, re-insert TP9
                
                %move TP9 to the correct position in chanlocs struture
                EEG.chanlocs(39:69)= EEG.chanlocs(38:68);
                EEG.chanlocs(38)= EEG.chanlocs(69);
                EEG.chanlocs(69)=[];
                
                %move TP9 to the correct position in data structure
                EEG.data(39:69,:,:)=EEG.data(38:68,:,:);
                EEG.data(38,:,:)=EEG.data(69,:,:);
                EEG.data(69,:,:)=[];
            end
            
            %baseline correction
            EEG = pop_rmbase( EEG, [-100   -2]);
            
            TEPs_SSPSIR= EEG; %to visualize final results
            TEPs_SSPSIR.icachansind=[];
            eeglab redraw
            
            %--------------------------------------------------------------
            %save
            EEG.setname= 'TEPs_SSPSIR';
            ALLEEG(3).setname= 'TEPs_SSPSIR';
            EEG = pop_saveset( EEG, 'filename',[subj '_TEPs_SSPSIR_' num2str(trigger(t)) '_' num2str(iii-1) '.set'],'filepath',SSP_SIRfolder);
            %--------------------------------------------------------------
            
        end
        
        
    end
    
    clear iii  t
    
    
    
end