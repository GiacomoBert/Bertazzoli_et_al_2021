function feature = extracttffeatures(data, times)

warning off;

data = double(data);
baseline = (times <= -100) & (times >= -300);
database = data(baseline,:)'; % 
databasevec = database(:);

% %% FOR RECHARGE ARTIFACT
% rechargeT = (times <= 40) & (times >= 1);
% datarecharge = data(rechargeT,:)';

%% POST-TMS
T = find(times>=1);
datapost = data(T,:);
datavec = datapost(:)';

%% DOWNSAMPLE TO 200HZ
fs = 1000/(times(2)-times(1));
factor = max(floor(fs/200),1);
datavec = datavec(1:factor:end);
fs2 = round(fs/factor);


%% CALCULATE FEATURES
% Power within 0~30 ms postTMS
T = (times >= 1) & (times <= 30);
stimabs = log(mean(abs(mean(data(T,:), 2)), 1)/std(databasevec));

% Proc Spectrum for Channel
[pxx, freq] = pwelch(datavec, ones(1, fs2), [], fs2, fs2);
pxx = 10*log10(pxx * fs2/2);

% The average log band power between 8 and 13 Hz
p = 0;
for i = 8:13
    p = p + pxx(find(freq == i,1));
end
Hz8_13 = p / (13-8+1);

% The average log band power between 14 and 20 Hz
p = 0;
for i = 14:20
    p = p + pxx(find(freq == i,1));
end
Hz14_20 = p / (20-14+1);

% The average log band power between 14 and 30 Hz
p = 0;
for i = 14:30
    p = p + pxx(find(freq == i,1));
end
Hz14_30 = p / (30-14+1);

% The average log band power between 14 and 30 Hz
p = 0;
for i = 31:50
    p = p + pxx(find(freq == i,1));
end
Hz31_50 = p / (50-31+1);

% The average log band power between 5 and 7 Hz
p = 0;
for i = 4:7
    p = p + pxx(find(freq == i,1));
end
Hz4_7 = p / (7-4+1);

% The average log band power between 65 and 80 Hz
p = 0;
for i = 65:80
    p = p + pxx(find(freq == i,1));
end
Hz65_80 = p / (80-65+1);

% The average log band power at 60 Hz
p = 0;
for i = 60:60
    p = p + pxx(find(freq == i,1));
end
Hz60 = p / (60-60+1);

% lambda and FitError: deviation of a component's spectrum from
% a protoptypical 1/frequency curve
p1.x = 2; %first point: value at 2 Hz
p1.y = pxx(find(freq == p1.x,1));

p2.x = 3; %second point: value at 3 Hz
p2.y = pxx(find(freq == p2.x,1));

%third point: local minimum in the band 5-13 Hz
p3.y = min(pxx(find(freq == 5,1):find(freq == 13,1)));
p3.x = freq(find(pxx == p3.y,1));

%fourth point: min - 1 in band 5-13 Hz
p4.x = p3.x - 1;
p4.y = pxx(find(freq == p4.x,1));

%fifth point: local minimum in the band 33-39 Hz
p5.y = min(pxx(find(freq == 33,1):find(freq == 39,1)));
p5.x = freq(find(pxx == p5.y,1));

%sixth point: min + 1 in band 33-39 Hz
p6.x = p5.x + 1;
p6.y = pxx(find(freq == p6.x,1));

%fifth point: local minimum in the band 65-80 Hz
p7.y = min(pxx(find(freq == 65,1):find(freq == 80,1)));
p7.x = freq(find(pxx == p7.y,1));

%sixth point: min + 1 in band 65-80 Hz
p8.x = p7.x + 1;
p8.y = pxx(find(freq == p8.x,1));

pX = [p1.x; p2.x; p3.x; p4.x; p5.x; p6.x; p7.x; p8.x];
pY = [p1.y; p2.y; p3.y; p4.y; p5.y; p6.y; p7.y; p8.y];

myfun = @(x,xdata)(exp(x(1))./ xdata.^exp(x(2))) - x(3);
xstart = [4, -2, 54];
fittedmodel = lsqcurvefit(myfun,xstart,double(pX),double(pY), [], [], optimset('Display', 'off'));

%FitError: mean squared error of the fit to the real spectrum in the band 2-40 Hz.
f1 = 8; f2 = 15;
ts_8to15 = freq(find(freq == f1) : find(freq == f2));
fs_8to15 = pxx(find(freq == f1) : find(freq == f2));
fiterror = log(norm(myfun(fittedmodel, ts_8to15)-fs_8to15)^2);

%FitError: mean squared error of the fit to the real spectrum in the band 2-40 Hz.
f1 = 65; f2 = 80;
ts_65to80 = freq(find(freq == f1) : find(freq == f2));
fs_65to80 = pxx(find(freq == f1) : find(freq == f2));
fiterror_65to80 = log(norm(myfun(fittedmodel, ts_65to80)-fs_65to80)^2);

%lambda: parameter of the fit
lambda = fittedmodel(2);

% Averaged local skewness 15s
T = (times >= 1);
abs_local_scewness = abs(skewness(data(T,:)));
mean_abs_local_scewness_15 = log(mean(abs_local_scewness));

% TEP 180 - 200 amp & 100-120 amp
T = (times >= 141) & (times <= 220);
TEP200 = log(mean(abs(mean(data(T,:), 2)), 1)/std(databasevec));
T = (times >= 61) & (times <= 120);
TEP100 = log(mean(abs(mean(data(T,:), 2)), 1)/std(databasevec));
T = (times >= 31) & (times <= 60);
TEP45 = log(mean(abs(mean(data(T,:), 2)), 1)/std(databasevec));

T = (times >= 61) & (times <= 100);
TEP80 = log(mean(mean(abs(data(T,:)), 1), 2)/std(databasevec));
T = (times >= 161) & (times <= 200);
TEP180 = log(mean(mean(abs(data(T,:)), 1), 2)/std(databasevec));
T = (times >= 141) & (times <= 180);
TEP160 = log(mean(mean(abs(data(T,:)), 1), 2)/std(databasevec));
T = (times >= 41) & (times <= 80);
TEP60 = log(mean(mean(abs(data(T,:)), 1), 2)/std(databasevec));
T1 = (times >= 81) & (times <= 120);
T2 = (times >= 181) & (times <= 220);
TEPdiff100200 = log(abs(mean(mean(data(T1,:), 1), 2)-mean(mean(data(T2,:), 1), 2))/std(databasevec));
TEP200max = max([TEP200, TEP180, TEP160]);
TEP100max = max([TEP100, TEP80]);

% kurtosis
kurtosis = trimmean(kurt(data), 1);

% maximum amplitude
maxamp = log(max(abs(datavec)));

% wavlet detection of EKG artifact
wt = modwt(data(:)/std(datavec),5);
wtrec = zeros(size(wt));
wtrec(6,:) = wt(6,:);
y = imodwt(wtrec, 'sym4');
y = abs(y).^2;
[pks1, ~]= findpeaks(y, 'MINPEAKHEIGHT', 2, 'MINPEAKDISTANCE', round(600/1000*fs));
[pks2, ~]= findpeaks(y, 'MINPEAKHEIGHT', 2, 'MINPEAKDISTANCE', round(50/1000*fs));
if (length(pks1) > 0.8*size(data,2)*size(data,1)/fs) && (length(pks2) - length(pks1) < 1.1*size(data,2)*size(data,1)/fs)
    EKG = 1;
else
    EKG = 0;
end


% time entropy average
timeentropy = zeros(1, size(data, 2));
for ii = 1:size(data, 2)
    timeentropy(ii) = entropy(data(:,ii));
end
timeentropyavg = mean(timeentropy);

% % recharge magnitude
% rechargemag = log(mean(mean(abs(datarecharge)))/std(databasevec));

% Append Features
feature = [stimabs, maxamp, TEP45, TEP200, TEP100, mean_abs_local_scewness_15, Hz31_50, Hz14_30, Hz8_13, Hz4_7, lambda, fiterror, EKG];%rechargemag

