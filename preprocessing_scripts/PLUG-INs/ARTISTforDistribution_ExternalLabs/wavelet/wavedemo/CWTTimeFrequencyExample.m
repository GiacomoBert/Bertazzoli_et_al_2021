%% Time-Frequency Analysis with the Continuous Wavelet Transform 
%
% This example shows how to use the continuous wavelet transform (CWT) to
% analyze signals jointly in time and frequency. The example explains the
% relationship between scale and frequency. The example discusses the
% localization of transients where the CWT outperforms the short-time
% Fourier transform (STFT). The example also shows how to synthesize
% time-frequency localized signal approximations using the inverse CWT.
%
% The CWT is compared with the STFT in a number of the examples. You must
% have the Signal Processing Toolbox(TM) in order to run the |spectrogram|
% examples.

%% Converting Scale to Frequency
% The key relationship for converting scale to frequency is $f=f_0/(s
% \Delta t)$ where $f_0$ is the center frequency of the wavelet in
% cycles/unit time, $s$ is the scale, and $\Delta t$ is the sampling
% interval. Both <matlab:doc('cwtft') |cwtft|> and <matlab:doc('cwt')
% |cwt|> convert scales to frequencies for you and output the result. For
% wavelets supported by |cwt|, you can find the center frequency using
% |centfrq| and convert scale to frequency using |scal2frq|. For example,
% find the frequency corresponding to a scale of 0.006 for the complex
% Morlet wavelet with a center frequency of 1.5 Hz and bandwidth
% parameter of 1. Note the scale in units of time is $s \Delta t$. Assume a
% a sampling period (interval) of 1 millisecond.

f0 = centfrq('cmor1-1.5');
s = 6;
dt = 0.001;
Freq = scal2frq(s,'cmor1-1.5',dt);
f = f0/(s*dt);

%%
% Compare the value of |Freq| and |f|. They are both equal to 250 Hz. Note
% that frequency is not simply the reciprocal of the scale.
%
% For wavelets supported by |cwtft|, the function |cwtftinfo| contains the
% information about scale-to-frequency conversion. You can also find the
% information in the |cwtft| reference. 

%% Selecting a Wavelet and Scales
% Using the CWT for time-frequency analysis works best when the magnitude
% of the wavelet's Fourier transform is symmetric with respect to its
% center frequency. Additionally, use a complex-valued wavelet whose
% Fourier transform is nonzero only on the positive frequencies. Based on
% these criteria, the best wavelets for time-frequency analysis, in
% general, are analytic wavelets such as the bump and analytic Morlet
% wavelets supported by |cwtft|.
%
% After you select the wavelet, you must select the scales over which to
% compute the CWT.
%
% Both in time and frequency, the dilation (stretching and shrinking)
% operation on the wavelet is multiplicative. Therefore, the natural
% spacing of scales for the CWT is logarithmic where there is a constant
% ratio between successive elements. This differs from linearly-spaced
% scales where there is a constant difference between adjacent elements.
% For the discrete wavelet and wavelet packet transforms, the ratio is 2,
% so that the scales are powers of two. For the CWT, choose any ratio
% greater than one. Common choices are $2^{1/12}$, $2^{1/16}$, $2^{1/32}$
% where 12, 16, and 32 are referred to as the number of voices per octave.
% To construct a logarithmically spaced scale vector, choose the number of
% voices per octave and total number of octaves. You want the range of
% normalized scales (assuming the sampling period to be one) to be greater
% than one and less than the length of your input signal. The following
% assumes a unit sampling period and constructs a vector of scales with 16
% voices per octave to cover a range of approximately 8 octaves.

numOctave = 8;
numVoices = 16;
a0 = 2^(1/numVoices);
scales = a0.^(1:numOctave*numVoices);
fprintf('Number of octaves is %1.3f\n',max(log2(scales))-min(log2(scales)));
%%
% The largest scale is 256, which is appropriate as long as the input
% signal has 256 samples.

%%
% Alternatively, you can choose a smallest scale, $s_0$, and construct the
% scale vector as follows. A common choice for $s_0$ is $2 \Delta t$.

s0 = 2;
a0 = 2^(1/numVoices);
scales = s0*a0.^(0:numOctave*numVoices);
fprintf('Number of octaves is %1.3f\n',max(log2(scales))-min(log2(scales)));
%%
% The minimum scale is 2 and the maximum scale is 512.

%%
% The helper function |helperCWTTimeFreqVector| accepts a minimum desired
% postive frequency and maximum frequency and returns a logarithmically
% spaced scale vector for use with the CWT. 
%
% The following examples illustrate applications of the CWT for
% time-frequency analysis.

%% Time-Frequency Analysis of Modulated Signals
% Load a quadratic chirp signal and show a plot of its spectrogram. The
% signal's frequency begins at approximately 500 Hz at t = 0, decreases to
% 100 Hz at t=2, and increases back to 500 Hz at t=4. The sampling
% frequency is 1 kHz.

load quadchirp;
[S,F,T] = spectrogram(quadchirp,100,98,128,1000);
helperCWTTimeFreqPlot(S,T,F,'surf','STFT of Quadratic Chirp','Seconds','Hz')
%%
% Obtain a time-frequency plot of this signal using the CWT. Use the bump
% wavelet to obtain the CWT. Construct a scale vector with a minimum
% desired frequency of 20 Hz and maximum desired frequency of 500 Hz. Use
% 32 voices per octave. Plot the result as a function of time and
% frequency.

f0 = 5/(2*pi);
scales = helperCWTTimeFreqVector(20,500,f0,0.001,32);
cwtquadchirp = cwtft({quadchirp,0.001},'wavelet','bump','scales',scales);
helperCWTTimeFreqPlot(cwtquadchirp.cfs,tquad,cwtquadchirp.frequencies,...
    'surf','CWT of Quadratic Chirp -- Bump Wavelet','Seconds','Hz')

%% 
% The CWT with the bump wavelet produces a time-frequency analysis very
% similar to the STFT.
%
% Frequency and amplitude modulation occur frequently in natural signals.
% Use the CWT to obtain a time-frequency analysis of an echolocation pulse
% emitted by a big brown bat (Eptesicus Fuscus). The sampling interval is 7
% microseconds. Use the bump wavelet with logarithmically spaced scales and
% 32 voices per octave. Thanks to Curtis Condon, Ken White, and Al Feng of
% the Beckman Center at the University of Illinois for the bat data and
% permission to use it in this example.

load batsignal
t = 0:DT:(numel(batsignal)*DT)-DT;
f0 = 5/(2*pi);
scales = helperCWTTimeFreqVector(1e3,71428.57,f0,DT,32);
cwtbat = cwtft({batsignal,DT},'wavelet','bump','scales',scales);
helperCWTTimeFreqPlot(cwtbat.cfs,t.*1e6,cwtbat.frequencies./1000,...
    'surf','Bat Echolocation (CWT)','Microseconds','kHz')
%%
% Obtain and plot the STFT of the bat data.
[S,F,T] = spectrogram(batsignal,50,48,128,1/DT);
helperCWTTimeFreqPlot(S,T.*1e6,F./1e3,'surf','Bat Echolocation (STFT)',...
    'Microseconds','kHz')
%%
% For both the simulated and natural modulated signals, the CWT provides
% results similar to the STFT.

%% Detection of Transients in Oscillations Using the CWT
% There are certain situations in time-frequency analysis where the CWT can
% provide a more informative time-frequency transform than the short-time
% Fourier transform. One such situation occurs when the signal is corrupted
% by transients. The appearance and disappearance of these transients often
% has physical significance. Therefore, it is important to be able to
% localize these transients in addition to characterizing oscillatory
% components in the signal. To simulate this, create a signal consisting of
% two sine waves with frequencies of 150 and 200 Hz. The sampling frequency
% is 1 kHz. The sine waves have disjoint time supports. The 150-Hz sine
% wave occurs between 100 and 300 milliseconds. The 200-Hz sine wave occurs
% from 700 milliseconds to 1 second. Additionally, there are two transients
% at 222 and 800 milliseconds. The signal is corrupted by noise.

rng default;
dt = 0.001;
t = 0:dt:1-dt;
addNoise = 0.025*randn(size(t));
x = cos(2*pi*150*t).*(t>=0.1 & t<0.3)+sin(2*pi*200*t).*(t>0.7);
x = x+addNoise;
x([222 800]) = x([222 800 ])+[-2 2];
figure;
plot(t.*1000,x);
xlabel('Milliseconds'); ylabel('Amplitude');
%%
% Zoom in on the two transients to see that they represent disturbances in
% the oscillations at 150 and 200 Hz.

subplot(2,1,1)
plot(t(184:264).*1000,x(184:264));
ylabel('Amplitude')
title('Transients')
axis tight;
subplot(2,1,2)
plot(t(760:840).*1000,x(760:840));
ylabel('Amplitude')
axis tight;
xlabel('Milliseconds')

%%
% Obtain and plot the CWT using the analytic Morlet wavelet.

f0 = 6/(2*pi);
scales = helperCWTTimeFreqVector(1,500,f0,0.001,32);
cwty = cwtft({x,dt},'wavelet','morl','scales',scales,'padmode','symw');
figure
helperCWTTimeFreqPlot(cwty.cfs,t.*1e3,cwty.frequencies,...
    'contour','Analytic Morlet Wavelet','Milliseconds','Hz')


%%
% The analytic Morlet wavelet exhibits poorer frequency localization than
% the bump wavelet, but superior time localization. This makes the Morlet
% wavelet a better choice for transient localization. Plot the
% magnitude-squared fine scale coefficients to demonstrate the localization
% of the transients.


plot(t.*1000,abs(cwty.cfs(1,:)).^2);
xlabel('Milliseconds'); ylabel('Magnitude');
grid on; 
title('Analytic Morlet Wavelet -- Fine Scale Coefficients');
hold off
%%
% The wavelet shrinks to enable time localization of the transients with a
% high degree of accuracy while stretching to permit frequency localization
% of the oscillations at 150 and 200 Hz.

%%
% The STFT can only localize the transients to the width of the window. You
% have to plot the STFT in decibels (dB) to be able to visualize the
% transients.

[S,F,T] = spectrogram(x,50,48,128,1000);
surf(T.*1000,F,20*log10(abs(S)),'edgecolor','none'); view(0,90);
axis tight;
shading interp;
colorbar
xlabel('Time'), ylabel('Hz');
title('STFT')
%%
% The transients appear in the STFT only as a broadband increase
% in power. Compare short-time power estimates obtained from the STFT
% before (centered at 183 msec) and after (centered at 223 msec) the
% appearance of the first transient.

figure;
plot(F,20*log10(abs(S(:,80))));
hold on;
plot(F,20*log10(abs(S(:,100))),'r');
legend('T = 183 msec','T = 223 msec')
xlabel('Hz');
ylabel('dB');
hold off;
%%
% The STFT does not provide the ability to localize transients to the same
% degree as the CWT.
%% Removing A Time-Localized Frequency Component Using the Inverse CWT
% Create a signal consisting of exponentially weighted sine waves. There
% are two 25-Hz components -- one centered at 0.2 seconds and one centered
% at 0.5 seconds. There are two 70-Hz components -- one centered at 0.2 and
% one centered at 0.8 seconds. Clearly, the first 25-Hz and 70-Hz
% components co-occur in time.

t = 0:1/2000:1-1/2000;
dt = 1/2000;
x1 = sin(50*pi*t).*exp(-50*pi*(t-0.2).^2);
x2 = sin(50*pi*t).*exp(-100*pi*(t-0.5).^2);
x3 = 2*cos(140*pi*t).*exp(-50*pi*(t-0.2).^2);
x4 = 2*sin(140*pi*t).*exp(-80*pi*(t-0.8).^2);
x = x1+x2+x3+x4;
figure;
plot(t,x)
title('Superimposed Signal')

%%
% Obtain the CWT using the bump wavelet and display the result as a
% function of time and frequency.

s0 = 2;
a0 = 2^(1/32);
scales = (s0*a0.^(32:7*32)).*dt;
cwtx = cwtft({x,dt},'Scales',scales,'Wavelet',{'bump',[4 0.9]});
figure;
contour(t,cwtx.frequencies,abs(cwtx.cfs))
xlabel('Seconds'), ylabel('Hz');
title('Analytic CWT using Bump Wavelet')
hcol = colorbar;
hcol.Label.String = 'Magnitude';

%%
% Remove the 25 Hz component which occurs from approximately 0.05 to 0.35
% seconds by zeroing out the CWT coefficients. Use the inverse CWT
% (<matlab:doc('icwtft') |icwtft|>) to reconstruct an approximation to the
% signal.

freqidx = (cwtx.frequencies>19 & cwtx.frequencies<31);
timeidx = 100:700;
reconcwt = cwtx;
reconcwt.cfs(freqidx,timeidx) = 0;
xrec = icwtft(reconcwt);
%%
% Display the CWT of the reconstructed signal. The initial 25-Hz component
% is removed.

cwtxrec = cwtft({xrec,dt},'Scales',scales,'Wavelet',{'bump',[4 0.9]});
figure;
contour(t,cwtxrec.frequencies,abs(cwtxrec.cfs))
xlabel('Seconds'), ylabel('Hz');
title('1st 25-Hz Component Removed')
colorbar;
%%
% Plot the original signal and the reconstruction.

subplot(2,1,1);
plot(t,x);
title('Original Signal');
subplot(2,1,2);
plot(t,xrec)
title('Signal with first 25-Hz component removed');
%%
% Finally, compare the reconstructed signal with the original signal
% without the 25-Hz component centered at 0.2 seconds.
y = x2+x3+x4;
figure;
plot(t,xrec)
hold on
plot(t,y,'r--')
legend('Inverse CWT approximation','Original Signal Without 25-Hz');
hold off


%% Determining Exact Frequency Through the Analytic CWT
% When you obtain the wavelet transform of a sine wave using an analytic
% wavelet, the analytic CWT coefficients actually encode the frequency.
%
% To illustrate this, consider an otoacoustic emission obtained from a
% human ear. Otoacoustic emissions (OAEs) are emitted by the cochlea (inner
% ear) and their presence are indicative of normal hearing. Load and plot
% the OAE data. The data are sampled at 20 kHz.

load dpoae
plot(t.*1000,dpoaets)
xlabel('Milliseconds')
ylabel('Amplitude')
%%
% The emission was evoked by a stimulus beginning at 25 milliseconds and
% ending at 175 milliseconds. Based on the experimental parameters, the
% emission frequency should be 1230 Hz. Obtain and plot the CWT.

f0 = 5/(2*pi);
scales = helperCWTTimeFreqVector(1,2000,f0,1/20000,16);
dpoaeCWT = cwtft({dpoaets,1/20000},'wavelet','bump','scales',scales);
helperCWTTimeFreqPlot(dpoaeCWT.cfs,t.*1000,dpoaeCWT.frequencies,...
    'surf','CWT of OAE','Milliseconds','Hz')
%%
% You can investigate the time evolution of the OAE by finding the CWT
% coefficients closest in frequency to 1230 Hz and examining their
% magnitudes as a function of time. Plot the magnitudes along with time
% markers designating the beginning and end of the evoking stimulus.

[~,idx1230] = min(abs(dpoaeCWT.frequencies-1230));
cfsOAE = dpoaeCWT.cfs(idx1230,:);
plot(t.*1000,abs(cfsOAE));
hold on
AX = gca;
plot([25 25],[AX.YLim(1) AX.YLim(2)],'r')
plot([175 175],[AX.YLim(1) AX.YLim(2)],'r')
xlabel('msec'), ylabel('Magnitude');
title('CWT Coefficient Magnitudes')
%%
% There is some delay between the onset of the evoking stimulus and the
% OAE. Once the evoking stimulus is terminated, the OAE immediately begins
% to decay in magnitude.
%%
% Another way to isolate the emission is to use the inverse CWT to
% reconstruct a frequency-localized approximation in the time domain.
%
% Construct a frequency-localized emission approximation by extracting the
% CWT coefficients corresponding to frequencies between 1150 and 1350 Hz.
% Use these coefficients and invert the CWT. Plot the original data along
% with the reconstruction and markers indicating the beginning and end of
% the evoking stimulus.

Scidx = find(dpoaeCWT.frequencies>= 1150 & dpoaeCWT.frequencies <=1350);
emissionCWT = dpoaeCWT;
tmpcfs = zeros(size(emissionCWT.cfs));
tmpcfs(Scidx,:) = emissionCWT.cfs(Scidx,:);
emissionCWT.cfs = tmpcfs;
xrec = icwtft(emissionCWT);
figure
plot(t.*1000,dpoaets)
hold on;
plot(t.*1000,xrec,'r')
AX = gca;
ylim = AX.YLim;
plot([25 25],ylim,'k')
plot([175 175],ylim,'k')
xlabel('Milliseconds')
ylabel('Amplitude')
title('Frequency-Localized Reconstruction of Emission')
%%
% In the time-domain data, you clearly see how the emission ramps on and
% off at the application and termination of the evoking stimulus.
%
% It is important to note that even though a range of frequencies were
% selected for the reconstruction, the analytic wavelet transform actually
% encodes the exact frequency of the emission. To demonstrate this, take
% the Fourier transform of the emission approximation reconstructed from
% the analytic CWT.

xdft = fft(xrec);
freq = 0:2e4/numel(xrec):1e4;
xdft = xdft(1:numel(xrec)/2+1);
figure
plot(freq,abs(xdft));
xlabel('Hz'); ylabel('Magnitude')
title('Fourier Transform of CWT-Based Signal Approximation');
[~,maxidx] = max(abs(xdft));
fprintf('The frequency is %4.2f Hz\n',freq(maxidx));
%% Conclusions and Further Reading
% In this example you learned how to use the CWT to obtain a time-frequency
% analysis of a 1-D signal using an analytic wavelet with |cwtft|. You
% learned how to create logarithmically spaced scale vectors appropriate
% for the discretized version of the CWT. You saw examples of signals where
% the CWT provides similar results to the STFT and an example where the CWT
% can provide more interpretable results than the STFT. Finally, you
% learned how to reconstruct time-scale (frequency) localized
% approximations to a signal using |icwtft|.
%% 
% For more information on the CWT see the Wavelet Toolbox documentation.
%
% Reference: Mallat, S. "A Wavelet Tour of Signal Processing: The Sparse
% Way", Academic Press, 2009.
%% Appendix
% The following helper functions are used in this example.
%
% * <matlab:edit('helperCWTTimeFreqPlot.m') helperCWTTimeFreqPlot>
% * <matlab:edit('helperCWTTimeFreqVector.m') helperCWTTimeFreqVector>

displayEndOfDemoMessage(mfilename)













