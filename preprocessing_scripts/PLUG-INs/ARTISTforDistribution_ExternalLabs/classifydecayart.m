function artcomp = classifydecayart(EEG,cfg)
% Reject Big decay artifact
artcomp = [];
W = EEG.icaweights*EEG.icasphere;
Wm = pinv(W);
source = reshape(W*EEG.data(:,:),[size(W,1), size(EEG.data,2), size(EEG.data,3)]);
for ii = 1:size(Wm, 2)
    sourcesig = squeeze(mean(abs(source(ii,(EEG.times<60)&(EEG.times>=cfg.PulseLen),:)), 3));
    if isfield(cfg, 'decaythr')
       I = find(sourcesig>cfg.decaythr); % find components with amplitude greater than a threshold
    else
       I = find(sourcesig>30);
    end
    if ~isempty(I)
        artcomp = [artcomp ii];
    end
end