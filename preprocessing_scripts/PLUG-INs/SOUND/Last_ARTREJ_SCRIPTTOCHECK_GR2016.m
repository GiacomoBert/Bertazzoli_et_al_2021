
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

clear all;		% prepare a clean workspace with no variables
close all;		% close all open events

%eeglab %initialize eeglab

% define directories
% all directories must have been created
%addpath C:\Users\Marta\Documents\fieldtrip-20171129
%locpath=('C:\Users\Marta\Documents\eeglab14_1_1b\plugins\dipfit2.3\standard_BESA\standard-10-5-cap385.elp'); %Replace first part with your location of eeglab % locpath = your location of 'standard-10-5-cap385_NoEOG.elp' file (located in eeglab directory)
%datalocation =('C:\Users\Marta\Dropbox\SoloMe\AAA_LAVORO\Lavoro FbF\1 - FISM\3-Esecuzione progetto FISM\3_Dati\TMS-EEG'); % datalocation =  write the directory storing your raw data
%stepfolder='C:\Users\Marta\Dropbox\AAA_SoloMe\AAA_LAVORO\Lavoro FbF\1 - FISM\3-Esecuzione progetto FISM\4_Analisi\AnalisiSound\Corrected\Steps';
%figurefolder='C:\Users\Marta\Dropbox\AAA_SoloMe\AAA_LAVORO\Lavoro FbF\1 - FISM\3-Esecuzione progetto FISM\4_Analisi\AnalisiSound\Corrected\Figures';
%outputfolder='C:\Users\Marta\Dropbox\AAA_SoloMe\AAA_LAVORO\Lavoro FbF\1 - FISM\3-Esecuzione progetto FISM\4_Analisi\AnalisiSound\Corrected\Final'; % outputfolder = write the directory containing your processed data
%SSP_SIRfolder='C:\Users\Marta\Dropbox\AAA_SoloMe\AAA_LAVORO\Lavoro FbF\1 - FISM\3-Esecuzione progetto FISM\4_Analisi\AnalisiSound\Corrected\Steps\SSP_SIR';

datalocation =('C:\Users\Neurofisilogia\Documents\MATLAB\Giacomo\GR2016_analyses\all_raw_data'); % datalocation =  write the directory storing your raw data
stepfolder='C:\Users\Neurofisilogia\Documents\MATLAB\Giacomo\GR2016_analyses\Analyses';

figurefolder='C:\Users\Neurofisilogia\Documents\MATLAB\Giacomo\GR2016_analyses\Analyses\Figures';
outputfolder='C:\Users\Neurofisilogia\Documents\MATLAB\Giacomo\GR2016_analyses\Analyses\Final'; % outputfolder = write the directory containing your processed data
SSP_SIRfolder='D:\Data\Giacomo\GR2016\Analyses\SSP_SIRD:\Data\Giacomo\GR2016\Analyses\SSP_SIR';



%definire elenco soggetti

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

%assign the number of component to remve from each dataset
%FUTURE IMPLEMENTATION: load all the images in one window.
%ask you to input the numebr of component to remove
%store this number in a structure with the
%association dataset-num of component to remove
%theoretically this loop coul go in the main one

subjects=natsort(subjects);%put datasets in alphabetic order

%IF YOU HAVE ALREADY DONE THE COMPONENT SELECTION:
%load([SSP_SIRfolder '\SSP-SIRcomp_to_remove']);

%IF YOU ALREADY HAVE ALL THE COMPONENT YOU WANT TO REMOVE, use this loop (it does not show the figure):


% for isubj=1:numel(subjects)
%     subj = subjects{isubj};
%     fprintf('Processing: %s \n', subj);
%
%     comp_to_be_rem(isubj).dataset=subj;
%     comp_to_be_rem(isubj).rem_comp=input(['Number of component to remove from ' subj ' (1-6) : ']);
% end

% IF YOU don't HAVE ALL THE COMPONENT YOU WANT TO REMOVE, use this loop
% (it shows all the figures and asks you which one to remove):

for isubj=1:numel(subjects)
    subj = subjects{isubj};
    fprintf('Processing: %s \n', subj);


    comp_to_be_rem(isubj).dataset=subj;
   h1=openfig([subj '_132_PC5.fig'],'invisible'); % open figure
   figure(h1);
   title('6th component');

   h1=openfig([subj '_132_PC4.fig'],'invisible'); % open figure
   figure(h1);
   title('5th component');

   h1=openfig([subj '_132_PC3.fig'],'invisible'); % open figure
   figure(h1);
   title('4th component');

   h1=openfig([subj '_132_PC2.fig'],'invisible'); % open figure
   figure(h1);
   title('3rd component');

   h1=openfig([subj '_132_PC1.fig'],'invisible'); % open figure
   figure(h1);
   title('2nd component');

   h1=openfig([subj '_132_PC0.fig'],'invisible'); % open figure
   figure(h1);
   title('1st component');

    comp_to_be_rem(isubj).rem_comp=input(['Number of component to remove from ' subj ' (1-6) : ']);
    close all

end

save([SSP_SIRfolder '\SSP-SIRcomp_to_remove'],'comp_to_be_rem'); %save list


%% Start loop for processing
% type the subject or a list of subjects to be analyzed (in this case subjects will be sequentially processed)
for isubj=1:numel(subjects)
    subj = subjects{isubj};
    fprintf('Processing: %s \n', subj);
    Ncomp=comp_to_be_rem(isubj).rem_comp; %load number of component to be removed
    
    for t=1: length(trigger)
        
        
        %% to set the number of components removed in SSP-SIR
        %CHANGE!!!! make an if loop that
        % if is a file name es. if '001_GB_HT0_LP'
        % remove components es Ncomp = x
        
        %     % for trigger 132
        %     if sum(subj=='SC01_DX')== 7
        %         Ncomp = 3;
        %     elseif  sum(subj=='SC02_DX')== 7
        %         Ncomp = 4;
        %     elseif  sum(subj=='SC03_DX')== 7
        %         Ncomp = 2;
        %     elseif  sum(subj=='SC04_DX')== 7
        %         Ncomp = 1;
        %     elseif  sum(subj=='SC05_DX')== 7
        %         Ncomp = 2;
        %     elseif  sum(subj=='SC06_DX')== 7
        %         Ncomp = 1;
        %     elseif  sum(subj=='SC07_DX')== 7
        %         Ncomp = 0;
        %     elseif  sum(subj=='SC08_DX')== 7
        %         Ncomp = 0;
        %     elseif  sum(subj=='SC09_DX')== 7
        %         Ncomp = 1;
        %     elseif  sum(subj=='SC10_DX')== 7
        %         Ncomp = 2;
        %     elseif  sum(subj=='SC11_DX')== 7
        %         Ncomp = 1;
        %     elseif  sum(subj=='SC12_DX')== 7
        %         Ncomp = 1;
        %     elseif  sum(subj=='SC13_DX')== 7
        %         Ncomp = 0;
        %     elseif  sum(subj=='SC14_DX')== 7
        %         Ncomp = 2;
        %     elseif  sum(subj=='SC15_DX')== 7
        %         Ncomp = 1;
        %     elseif  sum(subj=='SC21_DX')== 7
        %         Ncomp = 2;
        %     elseif  sum(subj=='SP16_DX')== 7
        %         Ncomp =   2;
        %     elseif  sum(subj=='SP17_DX')== 7
        %         Ncomp =  2;
        %     elseif  sum(subj=='SP18_DX')== 7
        %         Ncomp =   1;
        %     elseif  sum(subj=='SP19_DX')== 7
        %         Ncomp =  1;
        %     end
        %
        
        %to create TEPs with incremental number of removed components, up to 5
        iii= Ncomp + 1;
        
        %% 7. Last steps of preprocessing after SSP-SIR
            STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
            EEG = pop_loadset('filename',[subj '_4Sound.set'],'filepath', stepfolder);
           [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
            EEG = pop_loadset('filename',[subj '_5ICA.set'],'filepath', stepfolder);
           [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
            EEG = pop_loadset('filename',[subj '_6SSP-SIR_' num2str(trigger(t)) '_' num2str(iii-1) '.set'],'filepath', stepfolder);
            [ALLEEG, EEG, ~] = eeg_store( ALLEEG, EEG, 0 );
            [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'retrieve',3,'study',0);
            eeglab redraw
        
        
        numset= 3; % 2 if you want to use data after sound, 3 if you want to use data after SSP-SIR
        EEG=ALLEEG(numset);
        ALLEEG(numset+1)= ALLEEG(numset); %to create a new dataset, on which run the preprocessing
        [ALLEEG, EEG , ~] = pop_newset(ALLEEG, EEG, 1,'retrieve',numset+1,'study',0);
        eeglab redraw
        
        % low pass filter 70 Hz
        EEG  = pop_basicfilter( EEG,  1:66 , 'Cutoff',  70, 'Design', 'butter', 'Filter', 'lowpass', 'Order',  4 ); % GUI: 02-Mar-2018 11:28:25
        [ALLEEG, EEG , ~] = pop_newset(ALLEEG, EEG, (numset+2),'setname','xx LP 70 Hz','gui','off');
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
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHECK IF THIS ARTIFACT REJ IS USEFUL
        
        %     % artefact rejection 2
        %     %to remove muscles during rejection and see EEG signal only
        %
        %     EEGm=EEG;
        %     muscles=EEG.data(68,:,:);
        %     EEG.data(68,:,:)=0;
        %
        %
        %     pop_eegplot(EEG, 1, 1, 0);
        %     waitfor( findobj('parent', gcf, 'string', 'UPDATE MARKS'), 'userdata');
        %
        %     ar2=EEG.reject.rejmanual;
        %     EEG.data(68,:,:)=muscles;
        %     EEG.ar2= find(ar2);
        %
        %     EEG = pop_rejepoch( EEG, [EEG.ar2] ,0);
        %     [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 4,'gui','off');
        %     eeglab redraw
        %
        % %  changes are saved in dataset(numeset+1), dataset(numset+2) is removed, datasetnumset+1) is
        % % retrieved to work on it
        % ALLEEG(numset+1)=EEG;
        % [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'retrieve', numset+1, 'study',0);
        % ALLEEG = pop_delset( ALLEEG, [numset+2:size(ALLEEG, 2)] );
        % [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',numset+1,'study',0);
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %baseline correction
        EEG = pop_rmbase( EEG, [-100   -2]);
        
        TEPs_SSPSIR= EEG; %to visualize final results
        TEPs_SSPSIR.icachansind=[];
        eeglab redraw
        
        %----------------------------------------------------------------
        %save
        EEG.setname= 'TEPs_SSPSIR';
        ALLEEG(3).setname= 'TEPs_SSPSIR';
        EEG = pop_saveset( EEG, 'filename',[subj '_TEPs_SSPSIR_' num2str(trigger(t)) '_' num2str(iii-1) '_ar.set'],'filepath',outputfolder);
        
        %---------------------------------------------------------------------
    end
    
    clear  t
end