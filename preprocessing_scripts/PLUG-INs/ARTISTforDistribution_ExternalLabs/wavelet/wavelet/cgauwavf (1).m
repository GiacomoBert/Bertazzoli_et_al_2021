function [psi,X] = cgauwavf(LB,UB,N,NumWAW)
%CGAUWAVF Complex Gaussian wavelet.
%   [PSI,X] = CGAUWAVF(LB,UB,N,P) returns values of the Pth 
%   derivative of the complex Gaussian function 
%   F = Cp*exp(-i*x)*exp(-x^2) on an N point regular grid for the 
%   interval [LB,UB]. Cp is such that the 2-norm of the Pth 
%   derivative of F is equal to 1. P can be integer values 
%   from 1 to 8.
%
%   Output arguments are the wavelet function PSI
%   computed on the grid X.
%   [PSI,X] = CGAUWAVF(LB,UB,N) is equivalent to
%   [PSI,X] = CGAUWAVF(LB,UB,N,1).
%
%   These wavelets have an effective support of [-5 5].
%
%   ----------------------------------------------------
%   If you have access to the Symbolic Math Toolbox(TM),
%   you may specify a value of P > 8. 
%   ----------------------------------------------------
%
%   See also GAUSWAVF, WAVEINFO.

%   M. Misiti, Y. Misiti, G. Oppenheim, J.M. Poggi 01-Jun-99.
%   Last Revision 08-May-2012.
%   Copyright 1995-2012 The MathWorks, Inc.

% Check arguments.
%-----------------
nbIn = nargin;
switch nbIn
    case 0 
        if ~exist('maple') , nmax = 8; else nmax = 45; end %#ok<EXIST>
        psi = nmax;
        % psi contains the number max for Gaussian wavelets.
        % This number depends of Symbolic Toolbox

    case 2
        error(message('Wavelet:FunctionInput:NotEnough_ArgNum'));

    case 3 , NumWAW = 1;

    case 4 
        if ischar(NumWAW)
            [~,NumWAW] = wavemngr('fam_num',NumWAW);
            NumWAW = wstr2num(NumWAW);
        end

    otherwise
        error(message('Wavelet:FunctionInput:TooMany_ArgNum'));
end
if errargt(mfilename,NumWAW,'int')
    error(message('Wavelet:FunctionArgVal:Invalid_ArgVal'));
end

% Compute values of the Complex Gauss wavelet.
X = linspace(LB,UB,N);  % wavelet support.
if find(NumWAW==1:8)
  X2 = X.^2;
  F0 = exp(-X2);
  F1 = exp(-1i*X);
  F2 = (F1.*F0)/(exp(-1/2)*2^(1/2)*pi^(1/2))^(1/2);
end

switch NumWAW
  case 1
    psi = F2.*(-1i-2*X)*2^(1/2);

  case 2
    psi = 1/3*F2.*(-3+4*1i*X+4*X2)*6^(1/2);

  case 3
    psi = 1/15*F2.*(7*1i+18*X-12*1i*X.^2-8*X.^3)*30^(1/2);
                
  case 4
    psi = 1/105*F2.*(25-56*1i*X-72*X.^2+32*1i*X.^3+16*X.^4)*210^(1/2);

  case 5
    psi = 1/315*F2.*(-81*1i-250*X+280*1i*X.^2+240*X.^3-80*1i*X.^4-32*X.^5)*210^(1/2);

  case 6
    psi = 1/3465*F2.*(-331+972*1i*X+1500*X.^2-1120*1i*X.^3-720*X.^4+192*1i*X.^5+64*X.^6)*2310^(1/2);

  case 7
    psi = 1/45045*F2.*(1303*1i+4634*X-6804*1i*X.^2-7000*X.^3+3920*1i*X.^4+2016*X.^5-448*1i*X.^6-128*X.^7)*30030^(1/2);

  case 8
    psi = 1/45045*F2.*(5937-20848*1i*X-37072*X.^2+36288*1i*X.^3+28000*X.^4-12544*1i*X.^5-5376*X.^6+1024*1i*X.^7+256*X.^8)*2002^(1/2);

  otherwise
    if ~exist('maple') %#ok<EXIST>
        msg = getWavMSG('Wavelet:ToolNeeded:SymbolicTBX');
        errargt(mfilename,msg,'msg')
        error(message('Wavelet:ToolNeeded:SymbolicTBX'));
    end
    y = sym('t');
    f = exp(-1i*y).*exp(-(y.*y));
    d = diff(f,NumWAW);
    for j = 1:length(X)
        t = X(j); %#ok<NASGU>
        psi(j) = eval(d);
    end
end
intL2 = sum(psi.*conj(psi));
norL2 = intL2(end)*(X(2)-X(1));
psi   = psi/sqrt(norL2);
