function artchan = identifyartchan(EEG, cfg)
% FIND BAD CHANNELS
%% FIND SINGLE BAD ELECTRODES
artchan = [];
EEGtemp = EEG;
EEGtemp.data = EEG.data - repmat(median(EEG.data, 1),[size(EEG.data, 1), 1, 1]); % REFERENCE TO MEDIAN
thr = 0.02;
if ~isfield(cfg, 'maxiterchanrej')
    cfg.maxiterchanrej = 100;
end
iter = 0;
while true
    maxcorr = maxneighborcorr(EEGtemp);
    artchan1= artchan;
    artchan = unique([artchan find(maxcorr>thr)]);
    if isequal(artchan1, artchan) || isempty(artchan)
        break;
    end
    EEG2 = EEG;
    EEG2 = eeg_interp(EEG2, artchan);
    EEGtemp.data = EEG.data - repmat(mean(EEG2.data, 1),[size(EEG.data, 1), 1, 1]);
    iter = iter + 1;
    if iter >= cfg.maxiterchanrej
        break;
    end
end

%% CAR
EEGtemp = eeg_interp(EEG, artchan);
EEG.data = EEG.data - repmat(mean(EEGtemp.data, 1),[size(EEG.data, 1), 1, 1]);
maxcorr = maxneighborcorr(EEG);
artchan = find(maxcorr>thr);

%% FIND BAD ELECTRODE CLUSTERS USING RANSAC
if ~isfield(cfg, 'isransac')
    cfg.isransac = 0;
end
if cfg.isransac
    times = find((EEG.times < 0)|(EEG.times > 50));
    EEG1 = pop_select(EEG, 'nochannel', artchan);
    remainchan = setdiff(1:size(EEG.data, 1), artchan);
    chanlocs = EEG1.chanlocs;
    locs = [cell2mat({chanlocs.X}); cell2mat({chanlocs.Y}); cell2mat({chanlocs.Z})];
    ransacSubset = floor(0.25*size(EEG1.data, 1));
    ransacSampleSize = 50;
    P = hlp_microcache('cleanchans', @calc_projector, locs, ...
        ransacSampleSize, ransacSubset);
    WRansac = size(EEG1.data, 3);
    ransacCorrelationsT = zeros(length(locs), WRansac);
    
    % Calculate each channel's correlation to its RANSAC reconstruction for
    % each window
    n = size(EEG1.data(:,times,:), 2);
    m = size(EEG1.data, 1);
    p = ransacSampleSize;
    Xwin = EEG1.data(:,times,:);
    parfor k = 1:WRansac
        ransacCorrelationsT(:, k) = ...
            calculateRansacWindow(squeeze(Xwin(:, :, k))', P, n, m, p);
    end
    ransacCorrelationThreshold = 0.75;
    flagged = ransacCorrelationsT < ransacCorrelationThreshold;
    badChannelsFromRansac = remainchan(sum(flagged, 2)/size(flagged, 2) > 0.4);
    
    artchan = [artchan badChannelsFromRansac];
end


%% Helper functions for findNoisyChannels
function P = calc_projector(locs, numberSamples, subsetSize)
% Calculate a bag of reconstruction matrices from random channel subsets

[permutedLocations, subsets] = getRandomSubsets(locs, subsetSize, numberSamples);
randomSamples = cell(1, numberSamples);
parfor k = 1:numberSamples
    tmp = zeros(size(locs, 2));
    slice = subsets(k, :);
    tmp(slice, :) = real(spherical_interpolate(permutedLocations(:, :, k), locs))';
    randomSamples{k} = tmp;
end
P = horzcat(randomSamples{:});

function [permutedLocations, subsets] = getRandomSubsets(locs, subsetSize, numberSamples)
stream = RandStream('mt19937ar', 'Seed', 435656);
numberChannels = size(locs, 2);
permutedLocations = zeros(3, subsetSize, numberSamples);
subsets = zeros(numberSamples, subsetSize);
for k = 1:numberSamples
    subset = randsample(1:numberChannels, subsetSize, stream);
    subsets(k, :) = subset;
    permutedLocations(:, :,  k) = locs(:, subset);
end

function Y = randsample(X, num, stream)
Y = zeros(1, num);
for k = 1:num
    pick = round(1 + (length(X)-1).*rand(stream));
    Y(k) = X(pick);
    X(pick) = [];
end
function rX = calculateRansacWindow(XX, P, n, m, p)
YY = sort(reshape(XX*P, n, m, p),3);
YY = YY(:, :, round(end/2));
rX = sum(XX.*YY)./(sqrt(sum(XX.^2)).*sqrt(sum(YY.^2)));
