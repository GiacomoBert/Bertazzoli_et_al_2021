function [yse,yre,yce] = helperCompare2DDenoising(numEX,t,PlotFlag) 
% This function helperCompare2DDenoising is only in support of
% DualtreeExample. It may change in a future release.
%   
%   
%   M. Misiti, Y. Misiti, G. Oppenheim, J.M. Poggi  26-Nov-2012.
%   Last Revision: 25-Jul-2013.
%   Copyright 1995-2013 The MathWorks, Inc.


if nargin<1 , numEX = 1; end

% Load image
%-----------
switch numEX
    case 1 , X = imread('glasses.jpg');
    case 2 , X = imread('wpeppers.jpg');
    case 3 , X = imread('ironthree.jpg');
    case 4 , load mask;
    case 5 , load laure;
    case 6 , load finger;
    case 7, load crtcol
    case 8, load porche
    case 9, load cathe_1.mat
    otherwise , load facets
end
if ismember(numEX,1:3)
    X = double(round(0.299*X(:,:,1) + 0.587*X(:,:,2) + 0.114*X(:,:,3)));
end

pause(0.1)
hWait = waitbar(0,'Please wait...');

rng default;
XN = X + 25*randn(size(X));

[se,se_psnr,se_mse,se_maxerr,se_L2RAT] = local_denoising('dwt',X,XN,t); %#ok<*NASGU>
[re,re_psnr,re_mse,re_maxerr,re_L2RAT] = local_denoising('realdt',X,XN,t);
waitbar(0.25,hWait);
[ce,ce_psnr,ce_mse,ce_maxerr,ce_L2RAT] = local_denoising('cplxdt',X,XN,t);
waitbar(0.5,hWait);
LW = 3;
FS = 9;
[mse,idx_se] = min(se); yse = helperDDDT2Deno('dwt',XN,t(idx_se)); %#ok<*ASGLU>
[mre,idx_re] = min(re); yre = helperDDDT2Deno('realdt',XN,t(idx_re));
waitbar(0.75,hWait);
[mce,idx_ce] = min(ce); yce = helperDDDT2Deno('cplxdt',XN,t(idx_ce));
%--------------------------------------------------------------------------
if (strcmp(PlotFlag,'PlotMetrics') || strcmp(PlotFlag,'plotmetrics'))
figure('Units','normalized','Position',[0.01  0.4  0.3  0.3]); 
hold on
plot(t,se,'m-','LineWidth',LW)
plot(t,re,'b--','LineWidth',LW)
plot(t,ce,'r:','LineWidth',LW)
hold off
grid
title('RMS Error vs. Threshold Value')
xlabel('Threshold Value');
ylabel('RMS Error');
legend(...
    'Standard 2-D', ...
    'Real Oriented 2-D', ...
    'Complex Oriented 2-D');
box on
%--------------------------------------------------------------------------
figure('Units','normalized','Position',[0.54  0.4  0.3  0.3]);
hold on
plot(t,se_psnr,'m-','LineWidth',LW)
plot(t,re_psnr,'b--','LineWidth',LW)
plot(t,ce_psnr,'r:','LineWidth',LW)
hold off
grid
title('PSNR vs. Threshold Value')
xlabel('Threshold Value');
ylabel('PSNR (dB)');
legend(...
    'Standard 2-D', ...
    'Real Oriented 2-D', ...
    'Complex Oriented 2-D', ...
    'Location','Best');
box on
elseif (strcmp(PlotFlag,'PlotImage') || strcmp(PlotFlag,'plotimage'))
%--------------------------------------------------------------------------


fig = figure('Units','Normalized','Position',[0.10  0.10  0.42  0.80]); 
set(fig,'DefaultAxesXTick',[],'DefaultAxesYTick',[],'DefaultAxesBox','On');
if ~exist('map','var')
    NBC = min([max(abs([X(:);XN(:)])),255]);
    map = gray(NBC);    
end
colormap(map);


ax = zeros(1,5);
ax(1) = subplot(3,2,1); imagesc(X);  title('Original Image');
ax(2) = subplot(3,2,2); imagesc(XN); title('Noisy Image'); 
ax(3) = subplot(3,2,3); imagesc(yse); 
title({'Denoised Image','Standard 2-D'}); 
xlabel({sprintf('THR: %7.2f',t(idx_se)) ' - ' sprintf('RMSE: %7.2f',mse), ...
    sprintf('PSNR: %7.2f',se_psnr(idx_se))},'FontSize',FS)

ax(4) = subplot(3,2,4); imagesc(yre);
title({'Denoised Image','Real Oriented 2-D Dual-Tree'});
xlabel({sprintf('THR: %7.2f',t(idx_re)) ' - ' sprintf('RMSE: %7.2f',mre), ...
    sprintf('PSNR: %7.2f',re_psnr(idx_re))},'FontSize',FS)
ax(5) = subplot(3,2,5); imagesc(yce);
title({'Denoised Image','Complex Oriented 2-D Dual-Tree'});
xlabel({sprintf('THR: %7.2f',t(idx_ce)) ' - ' sprintf('RMSE: %7.2f',mce), ...
    sprintf('PSNR: %7.2f',ce_psnr(idx_ce))},'FontSize',FS)
set(ax,'XTick',[],'YTick',[],'Box','On');
end

delete(hWait)

%--------------------------------------------------------------------------
function [err,psnr,mse,maxerr,L2RAT] = local_denoising(meth,X,XN,t)

N = length(t);
err = zeros(1,N);
psnr = zeros(1,N);
mse = zeros(1,N);
maxerr = zeros(1,N);
L2RAT = zeros(1,N);
for k = 1:N
    y = helperDDDT2Deno(meth,XN,t(k));
    err(k) = sqrt(mean(mean((y-X).^2)));
    [psnr(k),mse(k),maxerr(k),L2RAT(k)] = psnr_mse_maxerr(X,y);    
end
%--------------------------------------------------------------------------
