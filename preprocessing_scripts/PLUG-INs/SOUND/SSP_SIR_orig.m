function [data_correct, artifact_topographies, data_clean, P,filt_ker,SSP_SIR_operator] = SSP_SIR_orig(data, PC, LFM, M, badC, is_ave_ref,fs,timeAxis,TimeROI)

% This functions removes muscle artifacts from TMS-EEG signal using the 
% SSP-SIR approach. The SSP operator is estimated as described in MÄki et
% al. 2011. The SIR and the time-adapted artifact rejection are implemented
% as in Mutanen et al. 2015.

% data: the input data to be cleaned 
% PC: the number of artifact PCs to be removed
% LFM: the lead field used in SIR
% M: the truncation dimension used in the SIR step
% badC: totally rejected channels that are recovered in the SIR step.
% is_ave_ref: if you use average reference, set this to be 1. Otherwise 0.

%To use this function the data should be the standard Nextsim data format.
%The function expects 1450 Hz sampling frequency.

%NOTE: Data and LFM are expeted to have the same reference system!

%Identifying the good channels:
data = double(data);
goodC = setdiff(1:size(data,1),badC);
data = data(goodC,:);

if is_ave_ref
    data = data - repmat(mean(data,1),[size(data,1),1]);
end

%High-pass filtering the data from 100 Hz:

[b,a] = butter(2,100/(fs/2),'high'); 
data_high = (filtfilt(b,a,data'))';

%Estimating the relative-amplitude kernel of the artifact:
%tmp = data_high.^2;


[~, Tind1] = min(abs(timeAxis-TimeROI(1)));
[~, Tind2] = min(abs(timeAxis-TimeROI(2)));

tmp = zeros(size(timeAxis));
tmp(Tind1:Tind2) = 1;
x_scal = round(50/(1000/fs)); %estimating the RMS value in 50-ms sliding windows
x = int8(linspace(-x_scal/2, x_scal/2,x_scal));
for i=1:size(tmp,1)
  %  filt_ker(i,:) = conv(tmp(i,:),ones(1,x_scal),'same')/x_scal;
  filt_ker(i,:) = conv(tmp(i,:),normpdf(1:x_scal,x_scal/2,x_scal/20),'same')/x_scal;
end
filt_ker = sum(filt_ker,1)/size(tmp,1);
filt_ker = filt_ker/max(filt_ker);
%filt_ker = sqrt(filt_ker);
filt_ker = repmat(filt_ker,[size(data_high,1),1]);

%Estimating the artifact-supression matrix P:
[U,S,~] = svd(filt_ker.*data_high,'econ');

if isempty(PC)
tmpfig = figure(666);
 for i =1:(size(U,2));
%        [~, ~,~,~] = spectrogram(U(:,i)'*data,2^5,ceil(0.75*2^5),[],fs,'yaxis');

subplot(1,2,1);
tmpsig = U(:,i)'*data;
%        newtimef( tmpsig, length(tmpsig), [timeAxis(1) timeAxis(end)], fs, 0 );
 [tf,freqs,times] = timefreq(tmpsig,fs,'freqs',[0 500]);
 imagesc(timeAxis,freqs,abs(tf));
 size(tf)
% ylim([0,200])
%[~,~,~]= plot_time_freq(U(:,i)'*data,0,[-8,0]);
% [~, ~,~,~] = spectrogram(U(:,i)'*data,2^5,ceil(0.75*2^5),[],fs,'yaxis');
subplot(1,2,2);
 plot(diag(S))
 hold on;
 plot(i,S(i,i),'r.','MarkerSize',10)
 hold off;
%[~,t1] = min(abs(ti)); [~,t2] = min(abs(ti-20));
%[~,f1] = min(abs(fr-30)); [~,f2] = min(abs(fr-100));
%pause();

w = waitforbuttonpress;
if w == 0
    disp('Button click')
    break;
else
    disp('Key press')
end
%ASR(i) = sqrt(mean(mean(pf(t1:t2,f2:end),1),2))/sqrt(mean(mean(pf(t2:end,1:f1),1),2))
 end
% subtightplot(2,2,3:4);

 close(tmpfig);
PC = i -1
M = rank(data) -PC;
end


P = eye(size(data,1)) - U(:,1:PC)*(U(:,1:PC))';
artifact_topographies = U(:,1:PC);
%Suppressing the artifacts:
data_clean = P*data;

%Performing SIR for the suppressed data:
L = LFM(goodC,:);
if is_ave_ref
    L = L - repmat(mean(L,1),[size(L,1),1]);
end
PL = P*L;

tau_proj = PL*PL';
[U,S,V] = svd(tau_proj);
S_inv = zeros(size(S));
S_inv(1:M,1:M) = diag(1./diag(S(1:M,1:M)));
tau_inv = V*S_inv*U';
suppr_data_SIR = LFM*(PL)'*tau_inv*data_clean;

SSP_SIR_operator =  LFM*((PL)'*(tau_inv*P));

% %Performing SIR for the original data:
tau_proj = L*L';
[U,S,V] = svd(tau_proj);
S_inv = zeros(size(S));
S_inv(1:M,1:M) = diag(1./diag(S(1:M,1:M)));
tau_inv = V*S_inv*U';
orig_data_SIR = LFM*(L)'*tau_inv*data;

filt_ker = repmat(filt_ker(1,:),[size(orig_data_SIR,1),1]);


data_correct = filt_ker.*suppr_data_SIR + orig_data_SIR - filt_ker.*orig_data_SIR;

end

