function EEG = ARTISTMain_GIACOMO(EEG,cfg)

% ARTISTMain - ARTIST algorithm for fully automated TMS-EEG artifact
%              rejection %find CHANGES for changes to adapt pipeline to
%              the other datasets
% Usage:
% >> EEG = ARTISTMain(EEG, cfg);
%
% Inputs:
%   EEG                  - Raw continuous EEG data in EEGLAB data structure (for one condition)
%   cfg                  - Configuration variable with hyperparameters
% Required fields:
%   cfg.EventCode        - Event code for the TMS pulses in cfg.EventCode (needs to be a string or a numeric value)
%   cfg.TMSEEGrootFolder - Path of the folder to store results
%   Optional fields:
%   cfg.TrialStart       - Start of each epoch in ms (default: -1000)
%   cfg.TrialEnd         - End of each epoch in ms (default: 1500)
%   cfg.PulseLen         - Length of the TMS pulse to be removed and interpolated in
%                          ms (default: 10)
%   cfg.PulseShift       - Time (ms) from zero where you want to start removing the
%                          TMS pulse (make sure it's before the rise of the pulse. Default: -2)
%   cfg.BaseLine         - Baseline time window in ms (default: [-300, -100])
%   cfg.plottimes        - Times (ms) at which topo plots are shown in the butterfly plots (default: [15,25,40,60,75,100,150,200,300])
%   cfg.NameProject      - Name of the project (default: '')
%   cfg.NameCond         - Name of the condition (default: '')
%   cfg.NameSub          - Name of the subject (default: '')
%
% Outputs:
%   EEG                  - clean epoched EEG data
%
% Copyright (C) 2017, Stanford University, W. Wu, K. Corey, and A. Etkin
%
% REFERENCE:
% Wu W, Keller C, Rogasch N, Longwell P, Shpigel E, Rolle C, Etkin A.
% ARTIST: a fully automated artifact rejection algorithm for single-pulse TMS-EEG data.
% Human Brain Mapping, 2018, DOI: 10.1002/hbm.23938.
% -------------
% FOR TESTING, LOAD IN DATASET AND CHANGE ROOT FOLDER
% load('exampleDataset.mat')
% cfg.TMSEEGrootFolder =
% '/Volumes/LabData/GoogleDrive/ARTISTforDistribution'
% cfg.plottimes = [15,25,40,60,75,100,150,200,300]; % ms FOR TOPOPLOT
% -------------
% SETUP - NEED TO HAVE EEGLAB AND WAVELET TOOLBOXES IN THE MATLAB PATH


clc
% SET UP INPUT VARIABLES
if isfield(cfg, 'EventCode')
    if ischar(cfg.EventCode) || isnumeric(cfg.EventCode)
        EventCode = cfg.EventCode;
    else
        error('EventCode NEEDS TO BE A STRING OR NUMBER!');
    end
else
    error('EVENT INFO MISSING!');
end
if isfield(cfg, 'TrialStart')
    TrialStart = cfg.TrialStart;
else
    TrialStart = -1000;% Default start of each epoch: -1000 ms
    cfg.TrialStart = TrialStart;
end
if isfield(cfg, 'TrialEnd')
    TrialEnd = cfg.TrialEnd;
else
    TrialEnd = 1500;% Default end of each epoch: 1500 ms
    cfg.TrialEnd = TrialEnd;
end
T = TrialEnd - TrialStart;
if isfield(cfg, 'PulseLen')
    PulseLen = cfg.PulseLen;
    cfg.TMSlength = cfg.PulseLen; % USED IN CLASSIFYDECAYART
else
    PulseLen = 10;% length of the TMS pulse to be removed and interpolated (ms)
    cfg.PulseLen = PulseLen;
    cfg.TMSlength = cfg.PulseLen; % USED IN CLASSIFYDECAYART
end
if isfield(cfg, 'PulseShift')
    PulseShift = cfg.PulseShift;
else
    PulseShift = -2;% time (ms) from zero where you want to start removing the TMS pulse (make sure it's before the rise of the pulse)
    cfg.PulseShift = PulseShift;
end

if isfield(cfg, 'BaseLine')
    BaseLine = cfg.BaseLine;
else
    BaseLine(1) = -300;
    BaseLine(2) = -100;% Baseline window
    cfg.BaseLine = BaseLine;
end
if ~isfield(cfg, 'plottimes')
    cfg.plottimes = [15,25,40,60,75,100,150,200,300];
end
if ~isfield(cfg, 'TMSEEGrootFolder')
    error('NEED A DIRECTORY TO STORE DATA')
end
if ~isfield(cfg, 'NameProject'); cfg.NameProject = ''; end
if ~isfield(cfg, 'NameCond'); cfg.NameCond = ''; end
if ~isfield(cfg, 'NameSub'); cfg.NameSub = ''; end

disp(['Using ' num2str(cfg.TrialStart) 'ms before TMS pulse'])
disp(['Using ' num2str(cfg.TrialEnd) 'ms after TMS pulse'])
disp(['Using ' num2str(cfg.BaseLine(1)) 'ms to ' num2str(cfg.BaseLine(2)) 'ms for baseline'])
disp(['Will interpolate ' num2str(cfg.PulseLen) 'ms around the TMS pulse'])
disp('If any of these settings are wrong, stop now and rerun');pause(2)

%% SET UP FOLDERS FOR SAVING VISUALIZATION
if ~isdir([cfg.TMSEEGrootFolder filesep cfg.NameProject]);mkdir([cfg.TMSEEGrootFolder filesep cfg.NameProject]);end;cd([cfg.TMSEEGrootFolder filesep cfg.NameProject])
if ~isdir([cfg.TMSEEGrootFolder filesep cfg.NameProject filesep cfg.NameSub]);mkdir([cfg.TMSEEGrootFolder filesep cfg.NameProject filesep cfg.NameSub]);end;cd([cfg.TMSEEGrootFolder filesep cfg.NameProject filesep cfg.NameSub])
if ~isdir([cfg.TMSEEGrootFolder filesep cfg.NameProject filesep cfg.NameSub filesep cfg.NameCond]);mkdir([cfg.TMSEEGrootFolder filesep cfg.NameProject filesep cfg.NameSub filesep cfg.NameCond]);end;cd([cfg.TMSEEGrootFolder filesep cfg.NameProject filesep cfg.NameSub filesep cfg.NameCond])
cfg.folderPath = ([cfg.TMSEEGrootFolder filesep cfg.NameProject filesep cfg.NameSub filesep cfg.NameCond]);
cfg.fullCondName = [cfg.NameProject '_' cfg.NameSub '_' cfg.NameCond];
if ~isdir([cfg.TMSEEGrootFolder filesep cfg.NameProject filesep cfg.NameSub filesep cfg.NameCond filesep 'ICA1']);mkdir([cfg.TMSEEGrootFolder filesep cfg.NameProject filesep cfg.NameSub filesep cfg.NameCond filesep 'ICA1']);end;
if ~isdir([cfg.TMSEEGrootFolder filesep cfg.NameProject filesep cfg.NameSub filesep cfg.NameCond filesep 'ICA2']);mkdir([cfg.TMSEEGrootFolder filesep cfg.NameProject filesep cfg.NameSub filesep cfg.NameCond filesep 'ICA2']);end;
if ~isdir([cfg.TMSEEGrootFolder filesep cfg.NameProject filesep cfg.NameSub filesep cfg.NameCond filesep 'QC']);mkdir([cfg.TMSEEGrootFolder filesep cfg.NameProject filesep cfg.NameSub filesep cfg.NameCond filesep 'QC']);end;

%% REMOVE THE PULSE ARTIFACT AND DOWNSAMPLE TO 1KHZ
disp('REMOVING THE PULSE ARTIFACT')
EEGbeforeRemoval = EEG; % SAVE FOR PLOTTING PURPOSES
cnt = 1;clear PulseStart PulseEnd
for kk = 1:length(EEG.event)
    if strcmp(num2str(EEG.event(kk).type), num2str(EventCode))
        TStart = ceil(EEG.event(kk).latency + TrialStart/1000*EEG.srate);
        if TStart <= 0
            error('Insufficient preTMS data for the first trial! Remove the first event');
        end
        TEnd = ceil(EEG.event(kk).latency + TrialEnd/1000*EEG.srate - 1);
        if TEnd > EEG.pnts
            error('Insufficient postTMS data for the last trial! Remove the last event');
        end
        PulseStart(cnt) = ceil(EEG.event(kk).latency + (PulseShift/1000*EEG.srate));
        PulseEnd(cnt) = ceil(PulseStart(cnt) + floor(PulseLen/1000*EEG.srate));
        x = ceil([TStart : PulseStart(cnt)-1,PulseEnd(cnt) + 1 : TEnd]);
        y = EEG.data(:, x);
        xi = floor(PulseStart(cnt) : PulseEnd(cnt));
        EEG.data(:,xi) = interp1(x, y', xi)';
        if cnt == 1
            firstPulseMS = EEG.event(kk).latency/(EEG.srate/1000);
        end
        cnt = cnt + 1;
    end
end

disp(['FOUND:' num2str(cnt) ' EVENTS TO REMOVE PULSE'])
cfg.srate_orig = EEG.srate;
if EEG.srate >= 1000
    disp('Downsampling to 1KHz')
    EEG = pop_resample(EEG , 1000);
else
    error('Sampling rate should be greater than 1000 Hz!');
end
disp('DOWNSAMPLING TO 1KHz')

%% PLOT THE PULSE ARTIFACT BEFORE AND AFTER PULSE REMOVAL REJECTION
close all; figure
plot(EEGbeforeRemoval.times-firstPulseMS,squeeze(EEGbeforeRemoval.data(10,:,:)),'b','LineWidth',1); hold all; box off;
plot(EEG.times-firstPulseMS,squeeze(EEG.data(10,:,:)),'r','LineWidth',1); hold all; box off;
line([PulseShift/(EEG.srate/1000) PulseShift/(EEG.srate/1000)+PulseLen],[500 500],'Color','k','LineStyle','--','LineWidth',2); hold all;
line([PulseShift/(EEG.srate/1000) PulseShift/(EEG.srate/1000)],[-500 500],'Color','k','LineStyle','--','LineWidth',2); hold all;
text(PulseShift/(EEG.srate/1000), 1500, [num2str(PulseShift/(EEG.srate/1000)) 'ms Pulse Shift, ' num2str(PulseLen) 'ms Pulse Duration']);
xlim([-20 50]); box off; xlabel('Time (ms)'); ylabel('uV');
% REMOVE UNDERLINES SO TITLE DOESN'T HAVE UNDERSCORES
tt = cfg.fullCondName;tm = strfind(tt,'_');for ji = 1:length(tm);tt(tm(ji)) = ' ';end;title(tt);
h = legend({'Before Pulse Removal';'After Pulse Removal'},'Location','SouthEast');legend boxoff
cd(cfg.folderPath);cd('QC');savefig_ARTIST([cfg.fullCondName '_QC_PulseRemoval'],16,16,150,'',4,[10 8]);

%% DETREND THE CONTINUOUS DATA
disp('DETRENDING CONTINUOUS DATA')
stimtrial = [];
for kk = 1:length(EEG.event)
    if strcmp(num2str(EEG.event(kk).type), num2str(EventCode))
        stimtrial = [stimtrial kk];
    end
end
for kk = 1:length(stimtrial)-1
    time1 = round(TrialStart/1000*EEG.srate + EEG.event(stimtrial(kk)).latency);
    time2 = round(TrialStart/1000*EEG.srate + EEG.event(stimtrial(kk + 1)).latency - 1);
    EEG.data(:, time1 : time2) = detrend(EEG.data(:, time1 : time2)','linear')';
end
% LAST EPOCH
time = round([TrialStart/1000*EEG.srate:TrialEnd/1000*EEG.srate] + EEG.event(stimtrial(end)).latency);
EEG.data(:, time) = detrend(EEG.data(:, time)','linear')';


%% EPOCHED DATA
disp('EPOCHING')
EEGepoch = pop_epoch(EEG, {EventCode}, [TrialStart/1000 TrialEnd/1000]);

%% 1ST STAGE: USE ICA ON THE CONTINUOUS DATA TO REMOVE BIG-AMPLITUDE DECAY ARTIFACT
disp('ICA ROUND 1 - REMOVE LARGE AMPLITUDE DECAY ARTIFACT');pause(1)

EEGepoch = pop_reref(EEGepoch,[]);
EEG = pop_reref(EEG,[]);

% INFOMAX ICA
EEG = pop_runica(EEG, 'icatype', 'runica', 'extended', 0, 'pca', EEG.nbchan - 1, 'interupt', 'off');
% COPY THE ICA PARAMETERS TO THE EPOCHED DATA
EEGepoch.icaweights = EEG.icaweights;
EEGepoch.icasphere = EEG.icasphere;
% DETECT BIG-AMPLITUDE DECAY COMPONENTS FROM THE EPOCHED DATA
artcomp = classifydecayart(EEGepoch, cfg);

% PLOT ICA MAPS
ARTIST_plotICA(EEGepoch, cfg, artcomp, 1)

EEG = pop_subcomp(EEG, artcomp);
EEG.artcomp1=artcomp; %CHANGES to save art comp
%% BANDPASS FILTER CONTINUOUS DATA BETWEEN 1 AND 90 Hz %CHANGES: lpass at 90Hz to match the other pipelines
disp('BANDPASS FILTER 1:90Hz');pause(1)
lcut = 1;hcut = 90;
EEG = pop_eegfiltnew(EEG, lcut, 0);
EEG = pop_eegfiltnew(EEG, 0, hcut);

%% NOTCH FILTERING TO REMOVE LINE NOISE %CHANGES: line at 50Hz (Italian line frequency)
disp('NOTCH FILTER (LINE NOISE REMOVAL) at 50 Hz');pause(1)
EEG = pop_eegfiltnew(EEG, 48, 52, 2000*EEG.srate/1000, 1);

%% 2ND STAGE: IDENTIFY AND REMOVE BAD TRIALS AND CHANNELS
disp('REMOVE BAD TRIALS AND CHANNELS');pause(1)
EEG = pop_epoch(EEG, {EventCode}, [TrialStart/1000 TrialEnd/1000]);
% IDENTIFY AND REMOVE BAD TRIALS
[arttrial] = identifyarttrial(EEG, cfg);
disp(['BAD TRIALS: ' num2str(arttrial')]); pause(1)
EEG = pop_select(EEG, 'notrial', arttrial);
EEG.arttrial=arttrial; %CHANGES to save art trials

% IDENTIFY AND REMOVE BAD CHANNELS
% Set cfg.isransac = 1 to remove bad electrode clusters. But Use it
% with caution as ransac tend to remove too many channels 
artchan = identifyartchan(EEG,cfg);
disp(['BAD CHANNELS: ' num2str(artchan)]); pause(1)
EEG = eeg_interp(EEG, artchan);
EEG.artchan=artchan; %CHANGES to save art channels
% %% CHANGES: save indices of bad epochs and bad channels
% save([cfg.TMSEEGrootFolder '\bad_epochs_' cfg.NameSub], 'arttrial');
% save([cfg.TMSEEGrootFolder '\bad_channels_' cfg.NameSub], 'artchan');

%% 3RD STAGE: REMOVE THE REMAINING ARTIFACTS
disp('ICA ROUND 2 - REMOVE REMAINING ARTIFACTS');pause(1)
% COMMON AVERAGE REFERENCE
EEG = pop_reref(EEG, []);
% DETERMINE OPTIMAL COMPONENT NUMBER NUMBER USING PCA
CovM = EEG.data(:, :)*EEG.data(:, :)'/size(EEG.data(:, :), 2);
[~, D] = eig(CovM);
d = sort(diag(D), 'descend');
dd = zeros(1, length(d));
for l = 1:length(d)
    dd(l) = sum(d(1:l));
end
cntNumCompForICA = sum(dd/dd(end) < .999);
% INFOMAX ICA
EEG = pop_runica(EEG, 'icatype', 'runica', 'extended', 0, 'pca', cntNumCompForICA, 'interupt', 'off');
% LOAD TRAINED CLASSIFIER
load classifierweight.mat;
cfg.w = classifier.w;
cfg.b = classifier.b;
artcomp = classifyartcomp(cfg, EEG);

% PLOT ICA MAPS
ARTIST_plotICA(EEG, cfg, artcomp, 2)

EEG = pop_subcomp(EEG, artcomp);
EEG.artcomp2=artcomp; %CHANGES to save ART ICA COMP 2

%% COMMON AVERAGE REFERENCE
disp('COMMON AVERAGE REFERENCING'); pause(1)
EEG = pop_reref(EEG, []);

%% BASELINE CORRECTION w.r.t -1000 ~ -2 ms %%CHANGES 
disp('BASELINE CORRECTION'); pause(1)
EEG = pop_rmbase(EEG, [BaseLine(1), BaseLine(2)]);

%% PLOT TEP AND TOPOGRAPHY
close all;figure;xlimm = [-200 500];colormap(jet);fnameTitle = cfg.fullCondName;
fnameTitletmp = strfind(fnameTitle,'_');for bb = 1:length(fnameTitletmp);fnameTitle(fnameTitletmp(bb)) = ' ';end
tit = sprintf([fnameTitle ' \n Total trials: ' num2str(size(EEG.data,3)) ', Srate: ' num2str(round(cfg.srate_orig)) 'Hz, \n '...
    num2str(length(artchan)) ' Bad Channels, ' num2str(length(arttrial)) ' Bad Trials, '...
    num2str(length(artcomp)) ' Bad Components']);
title(tit,'FontSize',10);
dat = squeeze(mean(EEG.data(:,EEG.times>xlimm(1) & EEG.times<xlimm(2),:),3));
timtopo(dat,EEG.chanlocs,[xlimm(1),xlimm(2),-1.5*max(max(dat)),1.5*max(max(dat))],cfg.plottimes,'',0,0,'shrink','on');box off
cd(cfg.folderPath);cd('QC');savefig_ARTIST([cfg.fullCondName '_TEP'],16,16,150,'',4,[10 8]);
disp('SAVING FILE...')
cd(cfg.folderPath);save(cfg.fullCondName,'EEG');close all
