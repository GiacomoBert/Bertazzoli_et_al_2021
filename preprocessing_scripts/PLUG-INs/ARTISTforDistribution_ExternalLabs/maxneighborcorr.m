function maxcorr = maxneighborcorr(EEG)
% Find the maximum correlation of each electrode with its neighbors
maxcorr = zeros(1, size(EEG.data, 1));
maxcorrwin = zeros(size(EEG.data, 1), size(EEG.data, 3));
times = find((EEG.times<0)|(EEG.times>50));
corrall = zeros(size(EEG.data, 1), size(EEG.data, 1), size(EEG.data, 3));
for jj =1:size(EEG.data, 3)
    corrall(:, :, jj) = real(corrm(squeeze(EEG.data(:,times,jj)), squeeze(EEG.data(:,times,jj))));
end
for c=1:size(EEG.data, 1)
    for jj = 1:size(EEG.data, 3)
        maxcorrwin(c, jj) = prctile(corrall(c, [1:c-1,c+1:end], jj), 98);
    end
    maxcorr(c) = sum(maxcorrwin(c,:) < 0.4)/size(EEG.data, 3);
end