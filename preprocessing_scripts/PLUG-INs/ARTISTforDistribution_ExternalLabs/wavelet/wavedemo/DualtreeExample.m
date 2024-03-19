%% Dual-Tree Wavelet Transforms
% This example shows how the dual-tree discrete wavelet transform (DWT)
% provides advantages over the critically sampled DWT for signal and image
% processing. The dual-tree DWT is implemented as two separate two-channel
% filter banks. To gain the advantages described in this example, you
% cannot arbitrarily choose the scaling and wavelet filters used in the two
% trees. The lowpass (scaling) and highpass (wavelet) filters of one tree,
% $\{h_0,h_1\}$, must generate a scaling function and wavelet that are
% approximate Hilbert transforms of the scaling function and wavelet
% generated by the lowpass and highpass filters of the other tree,
% $\{g_0,g_1\}$. Therefore, the complex-valued scaling function and wavelet
% formed from the two trees are approximately analytic.
%
% As a result, the dual-tree (complex) DWT exhibits less shift variance and
% more directional selectivity than the critically sampled DWT with only a
% $2^d$ redundancy factor for $d$-dimensional data. The redundancy in the
% dual-tree DWT is significantly less than the redundancy in the
% undecimated (stationary) DWT.
%
% This example illustrates the approximate shift invariance of the
% dual-tree DWT, the selective orientation of the dual-tree analyzing
% wavelets in 2-D, and the use of the dual-tree DWT in image denoising.

%% Near Shift Invariance of the Dual-Tree DWT
% The DWT suffers from shift variance, meaning that small shifts in the
% input signal or image can cause significant changes in the distribution
% of signal/image energy across scales in the DWT coefficients. The complex
% dual-tree DWT is approximately shift invariant. 

%%
% To demonstrate this on a test signal, construct two shifted discrete-time
% impulses 128 samples in length. One signal has the unit impulse at sample
% 60, while the other signal has the unit impulse at sample 64. Both
% signals clearly have unit energy ($\ell^2$ norm).

kronDelta1 = zeros(128,1);
kronDelta1(60) = 1;
kronDelta2 = zeros(128,1);
kronDelta2(64) = 1;

%%
% Obtain the DWT and dual-tree DWT of the two signals down to level 3 with
% wavelet and scaling filters of length 14. Extract the level-3 detail
% coefficients for comparison.

J = 3;   
dwt1 = dddtree('dwt',kronDelta1,J,'sym7');
dwt2 = dddtree('dwt',kronDelta2,J,'sym7');
dwt1Cfs = dwt1.cfs{J};
dwt2Cfs = dwt2.cfs{J};
      
dt1 = dddtree('cplxdt',kronDelta1,J,'dtf3');
dt2 = dddtree('cplxdt',kronDelta2,J,'dtf3');      
dt1Cfs = dt1.cfs{J}(:,:,1)+1i*dt1.cfs{J}(:,:,2);
dt2Cfs = dt2.cfs{J}(:,:,1)+1i*dt2.cfs{J}(:,:,2);

%%
% Plot the absolute values of the DWT and dual-tree DWT coefficients for
% the two signals at level 3 and compute the energy (squared $\ell^2$
% norms) of the coefficients.

figure;
subplot(1,2,1)
stem(abs(dwt1Cfs),'markerfacecolor',[0 0 1]);
title(['DWT Squared 2-norm = ' num2str(norm(dwt1Cfs,2)^2)]);
subplot(1,2,2)
stem(abs(dwt2Cfs),'markerfacecolor',[0 0 1])
title(['DWT Squared 2-norm = ' num2str(norm(dwt2Cfs,2)^2)]);

figure;
subplot(1,2,1)
stem(abs(dt1Cfs),'markerfacecolor',[0 0 1]);
title(['Dual-tree DWT Squared 2-norm = ' num2str(norm(dt1Cfs,2)^2)]);       
subplot(1,2,2)
stem(abs(dwt2Cfs),'markerfacecolor',[0 0 1])
title(['Dual-tree DWT Squared 2-norm = ' num2str(norm(dt2Cfs,2)^2)]);

%%
% Note the four sample shift in the signal has caused an almost 6.5% change
% in the energy of the level-3 DWT wavelet coefficients. However, the
% dual-tree wavelet coefficients show only a 0.3% change in energy.

%%
% To demonstrate the utility of approximate shift invariance in real data,
% we analyze an electrocardiogram (ECG) signal. The sampling interval for
% the ECG signal is 1/180 seconds. The data are taken from Percival &
% Walden (2000), p.125 (data originally provided by William Constantine and
% Per Reinhall, University of Washington). For convenience, we take the
% data to start at t=0.

load wecg
dt = 1/180;
t = 0:dt:(length(wecg)*dt)-dt;
figure;
plot(t,wecg)
xlabel('Seconds'); ylabel('Millivolts');

%%
% The large positive peaks approximately 0.7 seconds apart are the R waves
% of the cardiac rhythm. First, decompose the signal using the critically
% sampled DWT and plot the original signal along with the level-2 and
% level-3 wavelet coefficients. The level-2 and level-3 coefficients were
% chosen because the R waves are isolated most prominently in those scales.

J = 6; 
dtDWT1 = dddtree('dwt',wecg,J,'farras');
details = zeros(2048,3);
details(2:4:end,2) = dtDWT1.cfs{2};
details(4:8:end,3) = dtDWT1.cfs{3};
subplot(311)
stem(t,details(:,2),'Marker','none','ShowBaseline','off')
title('Level 2'); ylabel('mV');
subplot(312)
stem(t,details(:,3),'Marker','none','ShowBaseline','off')
title('Level 3'); ylabel('mV');
subplot(313)
plot(t,wecg); title('Original Signal');
xlabel('Seconds'); ylabel('mV');

%%
% Repeat the above analysis for the dual-tree transform. In this case, just
% plot the real part of the dual-tree coefficients at levels 2 and 3.

dtcplx1 = dddtree('cplxdt',wecg,J,'dtf3');
details = zeros(2048,3);
details(2:4:end,2) = dtcplx1.cfs{2}(:,1,1)+1i*dtcplx1.cfs{2}(:,1,2);
details(4:8:end,3) = dtcplx1.cfs{3}(:,1,1)+1i*dtcplx1.cfs{3}(:,1,2);
subplot(311)
stem(t,real(details(:,2)),'Marker','none','ShowBaseline','off')
title('Level 2'); ylabel('mV');
subplot(312)
stem(t,real(details(:,3)),'Marker','none','ShowBaseline','off')
title('Level 3'); ylabel('mV');
subplot(313)
plot(t,wecg); title('Original Signal');
xlabel('Seconds'); ylabel('mV');
%%
% Both the critically sampled and dual-tree wavelet transforms localize an
% important feature of the ECG waveform to similar scales.

%%
% An important application of wavelets in 1-D signals is to obtain an
% analysis of variance by scale. It stands to reason that this analysis of
% variance should not be sensitive to circular shifts in the input signal.
% Unfortunately, this is not the case with the critically sampled DWT. To
% demonstrate this, we circularly shift the ECG signal by 4 samples,
% analyze the unshifted and shifted signals with the critically sampled
% DWT, and calculate the distribution of energy across scales.

wecgShift = circshift(wecg,4);
dtDWT2 = dddtree('dwt',wecgShift,J,'farras');
sigenrgy = norm(wecg,2)^2;
enr1 = cell2mat(cellfun(@(x)(norm(x,2)^2/sigenrgy)*100,dtDWT1.cfs,'uni',0));
enr2 = cell2mat(cellfun(@(x)(norm(x,2)^2/sigenrgy)*100,dtDWT2.cfs,'uni',0));
levels = {'D1';'D2';'D3';'D4';'D5';'D6';'A6'};
enr1 = enr1(:);
enr2 = enr2(:);
table(levels,enr1,enr2,'VariableNames',{'Level','enr1','enr2'})


%%
% Note that the wavelet coefficients at levels 3 and 4 show approximately
% 3% changes in energy between the original and shifted signal. Next, we
% repeat this analysis using the complex dual-tree wavelet transform.

dtcplx2 = dddtree('cplxdt',wecgShift,J,'dtf3');
cfs1 = cellfun(@squeeze,dtcplx1.cfs,'uni',0);
cfs2 = cellfun(@squeeze,dtcplx2.cfs,'uni',0);
cfs1 = cellfun(@(x) complex(x(:,1),x(:,2)),cfs1,'uni',0);
cfs2 = cellfun(@(x) complex(x(:,1),x(:,2)),cfs2,'uni',0);
dtenr1 = cell2mat(cellfun(@(x)(norm(x,2)^2/sigenrgy)*100,cfs1,'uni',0));
dtenr2 = cell2mat(cellfun(@(x)(norm(x,2)^2/sigenrgy)*100,cfs2,'uni',0));
dtenr1 = dtenr1(:);
dtenr2 = dtenr2(:);
table(levels,dtenr1,dtenr2, 'VariableNames',{'Level','dtenr1','dtenr2'})

%%
% The dual-tree transform produces a consistent analysis of variance by
% scale for the original signal and its circularly shifted version.

%% Directional Selectivity in Image Processing
% The standard implementation of the DWT uses separable filtering of the
% columns and rows of the image. The LH, HL, and HH wavelets for
% Daubechies' least-asymmetric phase wavelet with 4 vanishing moments
% (sym4) are shown in the following figure.

figure;
J = 5;                      
L = 3*2^(J+1);
N = L/2^J;
Y = zeros(L,3*L);
dt = dddtree2('dwt',Y,J,'sym4');
dt.cfs{J}(N/3,N/2,1) = 1;
dt.cfs{J}(N/2,N/2+N,2) = 1;
dt.cfs{J}(N/2,N/2+2*N,3) = 1;
dwtImage = idddtree2(dt);
imagesc(dwtImage); axis xy; axis off;
title({'Critically Sampled DWT';'2-D separable wavelets (sym4) -- LH, HL, HH'});

%%
% Note that the LH and HL wavelets have clear horizontal and vertical
% orientations respectively. However, the HH wavelet on the far right mixes
% both the +45 and -45 degree directions, producing a checkerboard
% artifact. This mixing of orientations is due to the use of real-valued
% separable filters. The HH real-valued separable filter has passbands in
% all four high frequency corners of the 2-D frequency plane.

%%
% The dual-tree DWT achieves directional selectivity by using wavelets that
% are approximately analytic, meaning that they have support on only one
% half of the frequency axis. In the dual-tree DWT, there are six subbands
% for both the real and imaginary parts. The six real parts are formed by
% adding the outputs of column filtering followed by row filtering of the
% input image in the two trees. The six imaginary parts are formed by
% subtracting the outputs of column filtering followed by row filtering.

%%
% The filters applied to the columns and rows may be from the same filter
% pair, $\{h_0, h_1\}$ or $\{g_0, g_1\}$, or from different filter pairs,
% $\{h_0, g_1\}, \{g_0, h_1\}$. The following code shows the orientation of
% the 12 wavelets corresponding to the real and imaginary parts of the
% complex oriented dual-tree DWT.

J = 4;
L = 3*2^(J+1);
N = L/2^J;
Y = zeros(2*L,6*L);
wt = dddtree2('cplxdt',Y,J,'dtf3');
wt.cfs{J}(N/2,N/2+0*N,2,2,1) = 1;
wt.cfs{J}(N/2,N/2+1*N,3,1,1) = 1;
wt.cfs{J}(N/2,N/2+2*N,1,2,1) = 1;
wt.cfs{J}(N/2,N/2+3*N,1,1,1) = 1;
wt.cfs{J}(N/2,N/2+4*N,3,2,1) = 1;
wt.cfs{J}(N/2,N/2+5*N,2,1,1) = 1;
wt.cfs{J}(N/2+N,N/2+0*N,2,2,2) = 1;
wt.cfs{J}(N/2+N,N/2+1*N,3,1,2) = 1;
wt.cfs{J}(N/2+N,N/2+2*N,1,2,2) = 1;
wt.cfs{J}(N/2+N,N/2+3*N,1,1,2) = 1;
wt.cfs{J}(N/2+N,N/2+4*N,3,2,2) = 1;
wt.cfs{J}(N/2+N,N/2+5*N,2,1,2) = 1;
waveIm = idddtree2(wt);
imagesc(waveIm); axis off;
title('Complex Oriented Dual-Tree 2-D Wavelets');

%%
% The top row of the preceding figure shows the six directional wavelets of
% the real oriented dual-tree wavelet transform. The second row shows the
% imaginary parts. Together the real and imaginary parts form the complex
% oriented dual-tree wavelet transform. The real and imaginary parts are
% oriented in the same direction. You can use |dddtree2| with the
% |'realdt'| option to obtain the real oriented dual-tree DWT, which uses
% only the real parts. Using the real oriented dual-tree wavelet transform,
% you can achieve directional selectivity, but you do not gain the full
% benefit of using analytic wavelets such as approximate shift invariance.

%% Edge Representation in Two Dimensions
% The approximate analyticity and selective orientation of the complex
% dual-tree wavelets provide superior performance over the standard 2-D DWT
% in the representation of edges in images. To illustrate this, we analyze
% test images with edges consisting of line and curve singularities in
% multiple directions using the critically sampled 2-D DWT and the 2-D
% complex oriented dual-tree transform. First, analyze an image of an
% octagon, which consists of line singularities.

load woctagon;
figure;
imagesc(woctagon); colormap gray; 
title('Original Image'); axis equal; axis off;

%%
% Decompose the image down to level 4 and reconstruct an image
% approximation based on the level-4 detail coefficients.

dtcplx = dddtree2('cplxdt',woctagon,4,'dtf3');
dtDWT = dddtree2('dwt',woctagon,4,'farras');

dtcplx.cfs{1} = zeros(size(dtcplx.cfs{1}));
dtcplx.cfs{2} = zeros(size(dtcplx.cfs{2}));
dtcplx.cfs{3} = zeros(size(dtcplx.cfs{3}));
dtcplx.cfs{5} = zeros(size(dtcplx.cfs{5}));

dtDWT.cfs{1} = zeros(size(dtDWT.cfs{1}));
dtDWT.cfs{2} = zeros(size(dtDWT.cfs{2}));
dtDWT.cfs{3} = zeros(size(dtDWT.cfs{3}));
dtDWT.cfs{5} = zeros(size(dtDWT.cfs{5}));

dtImage = idddtree2(dtcplx);
dwtImage = idddtree2(dtDWT);
subplot(121)
imagesc(dtImage); axis equal; axis off; colormap gray;
title('Complex Oriented Dual-Tree');
subplot(122)
imagesc(dwtImage); axis equal; axis off; colormap gray;
title('DWT')
%%
% Next, analyze an octagon with hyperbolic edges. The edges in the
% hyperbolic octagon are curve singularities.

load woctagonHyperbolic;
figure;
imagesc(woctagonHyperbolic); colormap gray; 
title('Octagon with Hyperbolic Edges'); axis equal; axis off;

%%
% Again, decompose the image down to level 4 and reconstruct an image
% approximation based on the level-4 detail coefficients for both the
% standard 2-D DWT and the complex oriented dual-tree DWT.

dtcplx = dddtree2('cplxdt',woctagonHyperbolic,4,'dtf3');
dtDWT = dddtree2('dwt',woctagonHyperbolic,4,'farras');

dtcplx.cfs{1} = zeros(size(dtcplx.cfs{1}));
dtcplx.cfs{2} = zeros(size(dtcplx.cfs{2}));
dtcplx.cfs{3} = zeros(size(dtcplx.cfs{3}));
dtcplx.cfs{5} = zeros(size(dtcplx.cfs{5}));

dtDWT.cfs{1} = zeros(size(dtDWT.cfs{1}));
dtDWT.cfs{2} = zeros(size(dtDWT.cfs{2}));
dtDWT.cfs{3} = zeros(size(dtDWT.cfs{3}));
dtDWT.cfs{5} = zeros(size(dtDWT.cfs{5}));

dtImage = idddtree2(dtcplx);
dwtImage = idddtree2(dtDWT);
subplot(121)
imagesc(dtImage); axis equal; axis off; colormap gray;
title('Complex Oriented Dual-Tree');
subplot(122)
imagesc(dwtImage); axis equal; axis off; colormap gray;
title('DWT')

%%
% Note that the ringing artifacts evident in the 2-D critically sampled DWT
% are absent in the 2-D complex dual-tree transforms of both images. The
% complex oriented dual-tree DWT more faithfully reproduces line and curve
% singularities.

%% Image Denoising
% Because of the ability to isolate distinct orientations in separate
% subbands, the dual-tree DWT is often able to outperform the standard
% separable DWT in applications like image denoising. To demonstrate this,
% use the helper function |helperCompare2DDenoising|. The helper function
% loads an image and adds zero-mean white Gaussian noise with $\sigma =
% 25$. For a user-supplied range of thresholds, the function compares
% denoising using soft thresholding for the critically sampled DWT, the
% real oriented dual-tree DWT, and the complex oriented dual-tree DWT. For
% each threshold value, the root-mean-square (RMS) error and peak
% signal-to-noise ratio (PSNR) are displayed.

numex = 3;
helperCompare2DDenoising(numex,0:2:100,'PlotMetrics');

%%
% Both the real oriented and the complex oriented dual-tree DWTs outperform
% the standard DWT in RMS error and PSNR.

%%
% Next, obtain the denoised images for a threshold value
% of 25, which is equal to the standard deviation of the additive noise.

numex = 3;
helperCompare2DDenoising(numex,25,'PlotImage');

%%
% With a threshold value equal to the standard deviation of the additive 
% noise, the complex oriented dual-tree transform provides a PSNR almost
% 4 dB higher than the standard 2-D DWT.

%% Summary
% We have shown that the dual-tree DWT possesses the desirable properties
% of near shift invariance and directional selectivity not achievable with
% the critically sampled DWT. We have demonstrated how these properties can
% result in improved performance in signal analysis, the representation of
% singularities in images, and image denoising. In addition to the real
% oriented and complex oriented dual-tree DWT, |dddtree| and |dddtree2|
% also support the double-density wavelet transform and dual-tree
% double-density wavelet transforms, which are additional examples of
% overcomplete wavelet filter banks (frames) with advantages over the
% standard DWT.
%% Further Reading
%  Kingsbury, N.G. "Complex Wavelets for Shift Invariant Analysis and
%  Filtering of Signals". Journal of Applied and Computational Harmonic
%  Analysis. Vol 10, Number 3, May 2001, pp. 234-253.
%
%  Percival, D.B. and A.T. Walden. "Wavelet Methods for Time Series
%  Analysis", Cambridge University Press, 2000.
%
%  Selesnick, I., Baraniuk, R.G., and N.G. Kingsbury. "The Dual-Tree
%  Complex Wavelet Transform." IEEE Signal Processing Magazine. Vol. 22,
%  Number 6, November, 2005, pp. 123-151.
%
%  Selesnick, I. "The Double Density DWT". Wavelets in Signal and Image
%  Analysis: From Theory to Practice (A.A Petrosian, F.G. Meyer, eds.),
%  Norwell, MA: Kluwer Academic Publishers, 2001, pp. 39-66. 
%
%  Selesnick, I. "The Double-Density Dual-Tree Wavelet Transform". IEEE
%  Transactions on Signal Processing. Vol. 52, Number 5, May 2004, pp.
%  1304-1314.

%% Appendix
% The following helper functions are used in this example.
%
% * <matlab:edit('helperCompare2DDenoising.m') helperCompare2DDenoising.m>
% * <matlab:edit('helperDDDT2Deno.m') helperDDDT2Deno.m>

displayEndOfDemoMessage(mfilename)
      


