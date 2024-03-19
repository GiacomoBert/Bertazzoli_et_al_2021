function [arttrial] = identifyarttrial(EEG, cfg)
% FIND BAD TRIALS
time = (EEG.times >= 100)|(EEG.times <= 0);% Only look at baseline and beyond 100 ms 
trialpow = squeeze(mean(EEG.data(:,time,:).^2,2));
meanpow = repmat(median(trialpow, 2), 1, size(EEG.data, 3));
stdpow = repmat(iqr(trialpow, 2), 1, size(EEG.data, 3));
datTmp = squeeze(nanmean(abs(EEG.data),2));
% WARNING IF STD ACROSS EPOCHS IS HIGH
ChanWithLargeSTD = find(nanstd(datTmp')>30);
if ~isempty(ChanWithLargeSTD)
    disp(['WARNING: Channels with high (>30uV) STD: ' num2str(ChanWithLargeSTD)]);pause(2)
end
if isfield(cfg, 'trialthr')
    trialthr = cfg.trialthr;
else
    trialthr = 3;
end
artchannum = repmat(sum((trialpow-meanpow)./stdpow > trialthr, 1), size(EEG.data, 1), 1);
chanthr = round(0.2*EEG.nbchan);
[chan, trial] = find(((trialpow-meanpow)./stdpow > trialthr) & (artchannum < chanthr)); % only interplate trials with less than 10 noisy electrodes. Otherwise simply discard
arttrial = unique(trial);
for kk = 1:length(arttrial)
    artind = (trial == arttrial(kk));
    EEGart = pop_select(EEG, 'trial', arttrial(kk));
    EEGart = eeg_interp(EEGart,chan(artind));
    EEG.data(:, :, arttrial(kk)) = EEGart.data;
end
[~, trial] = find(((trialpow-meanpow)./stdpow > trialthr) & (artchannum >= chanthr)); 
arttrial = unique(trial);

