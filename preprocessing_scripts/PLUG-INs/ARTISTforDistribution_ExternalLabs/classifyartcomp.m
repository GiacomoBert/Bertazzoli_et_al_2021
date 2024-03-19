function [artcomp] = classifyartcomp(cfg, EEG)

%% CALCULATE THE ICA MIXING MATRIX AND SOURCE SIGNALS
W = EEG.icaweights*EEG.icasphere;
Wm = pinv(W);
source = reshape(W*EEG.data(:,:), [size(W,1), size(EEG.data,2), size(EEG.data,3)]);

%% GENERATE THE FEATURE FOR EACH COMPONENT
N = size(cfg.w, 1);
X = zeros(N, size(Wm, 2));
compamp = zeros(1, size(Wm, 2));
for kk = 1:size(Wm, 2)
    Spatial = Wm(:, kk);
    Temporal = squeeze(source(kk, :, :));
    compamp(kk) = max(max(abs(Temporal)));
    temporalfeature = extracttffeatures(Temporal, EEG.times);   
    spatialfeature = extractsfeatures(Spatial, EEG.chanlocs);
    X(:, kk) = [spatialfeature'; temporalfeature'];
end

%% DETECT BAD COMPONENTS
w = cfg.w;
b = cfg.b;

C = sign(w'*X+b*ones(1, size(X, 2)));
artcomp = find(C == -1);
artcomp = setdiff(artcomp, find(compamp <= 3));

