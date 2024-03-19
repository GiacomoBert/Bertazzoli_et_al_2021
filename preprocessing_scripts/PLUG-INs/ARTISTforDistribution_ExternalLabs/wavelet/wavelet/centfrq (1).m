function [freq,xval,recfreq] = centfrq(wname,iter,~)
%CENTFRQ Wavelet center frequency.
%   FREQ = CENTFRQ('wname') returns the center frequency in hertz
%   of the wavelet function 'wname' (see WAVEFUN).
%
%   For FREQ = CENTFRQ('wname',ITER), ITER is the number
%   of iterations used by the WAVEFUN function to compute
%   the wavelet.
%
%   [FREQ,XVAL,RECFREQ] = CENTFRQ('wname',ITER, 'plot')   
%   returns in addition the associated center frequency based 
%   approximation RECFREQ on the 2^ITER points grid XVAL 
%   and plots the wavelet function and RECFREQ.
%
%   See also SCAL2FRQ, WAVEFUN, WFILTERS. 

%   M. Misiti, Y. Misiti, G. Oppenheim, J.M. Poggi 04-Mar-98.
%   Last Revision: 20-Oct-2014.
%   Copyright 1995-2014 The MathWorks, Inc.

% Check arguments.
if nargin==1, iter = 8; end

% Retrieve wavelet.
wname = deblankl(wname);
wtype = wavemngr('type',wname);
switch wtype
  case 1 , [~,psi,xval] = wavefun(wname,iter);
  case 2 , [~,psi,~,~,xval] = wavefun(wname,iter);
  case 3 , [~,psi,xval] = wavefun(wname,iter);
  case 4 , [psi,xval] = wavefun(wname,iter);
  case 5 , [psi,xval] = wavefun(wname,iter);
end

T = max(xval)-min(xval);         % T is the size of the domain of psi.
n = length(psi);
psi = psi-mean(psi);             % psi is numerically centered.
psiFT = fft(psi);                % computation of the modulus
sp = (abs(psiFT));               % of the FT.

% Compute arg max of the modulus of the FT (center frequency).
Is_BIOR31 = isequal(wname,'bior3.1');
if ~Is_BIOR31
    [vmax,indmax] = max(sp);
    if indmax > n/2
        indmax = n-indmax+2;         % indmax is always >= 2.
    end
else
    [~,I,~] = localmax(sp(:)',1,false);
    indmax = I(1);                   % first local max 
    vmax = sp(indmax);    
end
per = T/(indmax-1);              % period corresponding to the maximum.		         
freq = 1/per;                    % associated frequency.

if nargin > 2
    % plots and  computation of the associated reconstructed signal.
    psiFT(sp<vmax) = 0;
    if Is_BIOR31 , psiFT(sp>vmax) = 0; end
    recfreq = ifft(psiFT);
    if wtype <= 4
        plot(xval,psi,'-b',xval, ...
            0.75*max(abs(psi))*real(recfreq)/max(abs(recfreq)),'-r'),
        title(['Wavelet ',wname,' (blue) and Center frequency',...
            ' based approximation'])
        xlabel(['Period: ',num2str(per), '; Cent. Freq: ', num2str(freq)])
    else
        subplot(211)
        plot(xval,real(psi),'-b',xval, ...
            0.75*max(abs(psi))*real(recfreq)/max(abs(recfreq)),'-r')
        title(['Wavelet ',wname,' (blue) and Center frequency',...
            ' based approximation'])
        ylabel('Real parts')
        xlabel(['Period: ',num2str(per), '; Cent. Freq: ', num2str(freq)])
        
        subplot(212),plot(xval,imag(psi),'-b',xval, ...
            0.75*max(abs(psi))*imag(recfreq)/max(abs(recfreq)),'-r')
        ylabel('Imaginary parts')
        xlabel(['Period: ',num2str(per), '; Cent. Freq: ', num2str(freq)])
    end
else
    recfreq = [];
end
