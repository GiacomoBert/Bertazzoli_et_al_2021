function [coefs,varargout] = cwt(SIG,scales,WAV,varargin)
%CWT Continuous 1-D wavelet transform.
%   COEFS = CWT(S,SCALES,WNAME) returns the continuous wavelet transform of
%   the real-valued signal S. The wavelet transform is computed for the
%   specified scales SCALES using the analyzing wavelet WNAME. SCALES is a
%   positive scalar or 1-D vector with positive elements. The string WNAME
%   denotes a wavelet recognized by WAVEMNGR. COEFS is a matrix with the
%   number of rows equal to the length of SCALES and number of columns
%   equal to the length of the input signal. The k-th row of COEFS
%   corresponds to the CWT coefficients for the k-th element in the scales
%   vector.
%
%   COEFS = CWT(S,SCALES,WNAME,PLOTMODE) plots the continuous wavelet 
%   transform coefficients.
%   Coefficients are colored using PLOTMODE.
%   PLOTMODE = 'lvl' (By scale) or 
%   PLOTMODE = 'glb' (All scales) 
%   PLOTMODE = 'abslvl' or 'lvlabs' (Absolute value by scale) or
%   PLOTMODE = 'absglb' or 'glbabs' (Absolute value over all scales)
%
%   CWT(...,'plot') is equivalent to CWT(...,'absglb')
%
%   You get 3-D plots (surfaces) using the same PLOTMODE options listed
%   above preceded by '3D'. For example: COEFS = CWT(...,'3Dplot') or COEFS
%   = CWT(...,'3Dlvl'). You must have at least two elements in SCALES in
%   order to use the '3Dplot' option.
%
%   COEFS = CWT(S,SCALES,WNAME,PLOTMODE,XLIM) use the coefficients within
%   the x-axis limits XLIM to determine the color scaling of the plot.
%   XLIM is a two-element vector [x1 x2] with 1 <= x1 < x2 <=
%   length(S).
%
%   COEFS = CWT(S,SCALES,WNAME,'scal') plots the scalogram for the
%   continuous wavelet transform. The scalogram is displayed as a image
%   plot. Use 'scalCNT' to obtain a contour plot of the scalogram. You must
%   have at least two elements in SCALES in order to use the 'scalCNT'
%   option.
%
%   [COEFS,SC] = CWT(..., 'scal') returns the scalogram. You must specify
%   either 'scal' or 'scalCNT' to output the scalogram. You must have at
%   least two elements in SCALES in order to use the 'scalCNT' option.
%
%   [COEFS,FREQUENCIES] = CWT(S,SCALES,WNAME,SAMPLINGPERIOD) returns the
%   frequencies in cycles per unit time corresponding to the scales and the
%   analyzing wavelet WNAME. SAMPLINGPERIOD is a positive real-valued
%   scalar. If the units of SAMPLINGPERIOD are seconds, the frequencies are
%   in hertz.
%
%   [COEFS,SC,FREQUENCIES] = CWT(S,SCALES,WNAME,SAMPLINGPERIOD,'scal')
%   returns the scalogram and the frequencies corresponding to the
%   scales and the analyzing wavelet. You can also use the flag 'scalCNT'
%   to output the scalogram if you have at least two elements in SCALES.
%   SAMPLINGPERIOD is only used in the conversion of scales to frequencies,
%   specifying SAMPLINGPERIOD does not affect the appearance of plots
%   generated by CWT.
%
%   %Example 1: 
%   %   Create a signal consisting of two disjoint sine waves with 
%   %   frequencies of 150 and 200 Hz with additive noise. The signal is 
%   %   sampled at 1 kHz. The signal has two transients located at 0.2 and 
%   %   0.8 seconds. Obtain the CWT with a complex Morlet wavelet and plot 
%   %   the result as a function of time and frequency.
%   %
%   dt = 0.001;
%   t = 0:dt:1-dt;
%   x = ...
%   cos(2*pi*150*t).*(t>=0.1 & t<0.3)+sin(2*pi*200*t).*(t>0.7);
%   y = x+0.05*randn(size(t));
%   y([200 800]) = y([200 800])+2;
%   a0 = 2^(1/32);
%   scales = 2*a0.^(0:6*32);
%   [cfs,frequencies] = cwt(y,scales,'cmor1-1.5',dt);
%   contour(t,frequencies,abs(cfs)); grid on;
%   xlabel('Seconds'); ylabel('Hz');
%
%   %Example 2
%   %   Obtain and plot the scalogram of the Kobe earthquake data recorded
%   %   in Australia on 16 January 1995. The sampling period is one second. 
%   %   Output the frequencies in Hz and use them to construct a 
%   %   time-frequency plot of the scalogram.
%   %   
%   load kobe
%   a0 = 2^(1/32);
%   scales = a0.^(4*32:8*32);
%   [cfs,sc,frequencies] = cwt(kobe,scales,'cmor1-1.5',1,'scal');
%   figure;
%   contour(1:numel(kobe),frequencies,sc)
%   ylabel('Hz'); xlabel('Sample');
%   grid on
%   title('Kobe Earthquake Data -- Time-Frequency Plot of Scalogram')
%
%   See also CWTFT, MODWT, SWT, WAVEDEC, WAVEFUN, WAVEINFO, WCODEMAT, WPDEC.

%   M. Misiti, Y. Misiti, G. Oppenheim, J.M. Poggi 12-Mar-96.
%   Copyright 1995-2015 The MathWorks, Inc.

nargoutchk(0,3);
% Check scales.
%--------------
err = 0;
if isempty(scales) ,         err = 1;
elseif min(size(scales))>1 , err = 1;
elseif min(scales)<eps,      err = 1;
end
if err
    errargt(mfilename, ...
        getWavMSG('Wavelet:FunctionArgVal:Invalid_ScaVal'),'msg');
    error(message('Wavelet:FunctionArgVal:Invalid_ScaVal'))
end

% Check signal.
%--------------
if isnumeric(SIG)
    ySIG    = SIG;
    lenSIG  = length(ySIG);
    xSIG    = (1:lenSIG);
    stepSIG = 1;
    
elseif isstruct(SIG)
    try
        ySIG = SIG.y;
    catch ME  %#ok<NASGU>
        err = 1;
    end
    if err~=1
        lenSIG = length(ySIG);
        try
            xSIG = SIG.x; stepSIG = xSIG(2)-xSIG(1);
        catch ME  %#ok<NASGU>
            try
                stepSIG = SIG.step;
                xSIG = (0:stepSIG:(lenSIG-1)*stepSIG);
            catch ME  %#ok<NASGU>
                try
                    xlim = SIG.xlim;
                    xSIG = linspace(xlim(1),xlim(2),lenSIG);
                    stepSIG = xSIG(2)-xSIG(1);
                catch ME  %#ok<NASGU>
                    xSIG = (1:lenSIG); stepSIG = 1;
                end
            end
        end
    end
    
elseif iscell(SIG)
    ySIG = SIG{1};
    xATTRB  = SIG{2};
    lenSIG  = length(ySIG);
    len_xATTRB = length(xATTRB);
    if len_xATTRB==lenSIG
        xSIG = xATTRB; 
        stepSIG = xSIG(2)-xSIG(1);

    elseif len_xATTRB==2
        xlim = xATTRB;
        xSIG = linspace(xlim(1),xlim(2),lenSIG);
        stepSIG = xSIG(2)-xSIG(1);

    elseif len_xATTRB==1
        stepSIG = xATTRB;
        xSIG = (0:stepSIG:(lenSIG-1)*stepSIG);
    else
        xSIG = (1:lenSIG); stepSIG = 1;
    end
else
    err = 1;
end
if err
    errargt(mfilename, ...
        getWavMSG('Wavelet:FunctionArgVal:Invalid_SigVal'),'msg');
    error(message('Wavelet:FunctionArgVal:Invalid_SigVal'))
end

% Parse variable input arguments
params = parseinputs(lenSIG,stepSIG,varargin{:});
DT = params.samplingPeriod;
validateattributes(DT,{'double'},{'positive'});
plotmode = params.plotmode; 


% Check wavelet.
%---------------
getINTEG = 1;
getWTYPE = 1;
if ischar(WAV)
    precis = 10; % precis = 15;
    [val_WAV,xWAV] = intwave(WAV,precis);
    stepWAV = xWAV(2)-xWAV(1);
    wtype = wavemngr('type',WAV);
    if wtype==5 , val_WAV = conj(val_WAV); end
    getINTEG = 0;
    getWTYPE = 0;

elseif isnumeric(WAV)
    val_WAV = WAV;
    lenWAV  = length(val_WAV);
    xWAV = linspace(0,1,lenWAV);
    stepWAV = 1/(lenWAV-1);
    
elseif isstruct(WAV)
    try
        val_WAV = WAV.y; 
    catch ME  %#ok<NASGU>
        err = 1; 
    end
    if err~=1
        lenWAV = length(val_WAV);
        try
            xWAV = WAV.x; stepWAV = xWAV(2)-xWAV(1);
        catch ME  %#ok<NASGU>
            try
                stepWAV = WAV.step;
                xWAV = (0:stepWAV:(lenWAV-1)*stepWAV);
            catch ME  %#ok<NASGU>
                try
                    xlim = WAV.xlim;
                    xWAV = linspace(xlim(1),xlim(2),lenWAV);
                    stepWAV = xWAV(2)-xWAV(1);
                catch ME  %#ok<NASGU>
                    xWAV = (1:lenWAV); stepWAV = 1;
                end
            end
        end
    end
    
elseif iscell(WAV)
    if isnumeric(WAV{1})
        val_WAV = WAV{1};
    elseif ischar(WAV{1})
        precis  = 10;
        val_WAV = intwave(WAV{1},precis);
        wtype = wavemngr('type',WAV{1});        
        getINTEG = 0;
        getWTYPE = 0;
    end
    xATTRB  = WAV{2};
    lenWAV  = length(val_WAV);
    len_xATTRB = length(xATTRB);
    if len_xATTRB==lenWAV
        xWAV = xATTRB; stepWAV = xWAV(2)-xWAV(1);

    elseif len_xATTRB==2
        xlim = xATTRB;
        xWAV = linspace(xlim(1),xlim(2),lenWAV);
        stepWAV = xWAV(2)-xWAV(1);

    elseif len_xATTRB==1
        stepWAV = xATTRB;
        xWAV = (0:stepWAV:(lenWAV-1)*stepWAV);
    else
        xWAV = linspace(0,1,lenWAV);
        stepWAV = 1/(lenWAV-1);
    end
end
if err
    errargt(mfilename, ...
        getWavMSG('Wavelet:FunctionArgVal:Invalid_WavVal'),'msg');
    error(message('Wavelet:FunctionArgVal:Invalid_WavVal'))
end
xWAV = xWAV-xWAV(1);
xMaxWAV = xWAV(end);
if getWTYPE ,  wtype = 4; end
if getINTEG ,  val_WAV = stepWAV*cumsum(val_WAV); end

ySIG   = ySIG(:)';
nb_SCALES = length(scales);
coefs     = zeros(nb_SCALES,lenSIG);
ind  = 1;

    
    
for k = 1:nb_SCALES
    a = scales(k);
    a_SIG = a/stepSIG;
    j = 1+floor((0:a_SIG*xMaxWAV)/(a_SIG*stepWAV));     
    if length(j)==1 , j = [1 1]; end
    f            = fliplr(val_WAV(j));
    coefs(ind,:) = -sqrt(a)*wkeep1(diff(wconv1(ySIG,f)),lenSIG);
    ind          = ind+1;
end

    frequencies = [];
    if ~isempty(DT) && ischar(WAV)
        frequencies = scal2frq(scales,WAV,DT);
        
    end
    
 switch nargout
    case 2
        if strncmpi(plotmode,'scal',1)
            SC = wscalogram(' ',coefs,scales,ySIG,xSIG);
            varargout{1} = SC;
        elseif (~strncmpi(plotmode,'scal',1) && ~isempty(frequencies))
            varargout{1} = frequencies;
        else
            error(message('Wavelet:FunctionOutput:UnsupportedOutCWT'));
        end
    case 3
        if (strncmpi(plotmode,'scal',1) && ~isempty(DT))
            SC = wscalogram(' ',coefs,scales,ySIG,xSIG);
            varargout{1} = SC;
            varargout{2} = frequencies;
        else
            error(message('Wavelet:FunctionOutput:UnsupportedOutCWT'));
        end
    
end   

% Test for plots.
%----------------
    if isempty(plotmode) , return; end

% Display Continuous Analysis.
%-----------------------------
dummyCoefs = coefs;
NBC = 240;
if strncmpi('3D',plotmode,2)
    dim_plot = '3D';
elseif strncmpi('scal',plotmode,4)
    dim_plot = 'SC';    
else
    dim_plot = '2D';
end

if (isequal(wtype,5) && ~strncmpi(plotmode,'scal',1))
   if ~isempty(strfind(plotmode,'lvl')) 
       plotmode = 'lvl';
   else
       plotmode = 'glb';   
   end
end

switch plotmode
    case {'lvl','3Dlvl'}
        lev_mode  = 'row';   abs_mode  = 0;   msg_Ident = 'By_scale';
        
    case {'glb','3Dglb'}
        lev_mode  = 'mat';   abs_mode  = 0;   msg_Ident = '';
        
    case {'abslvl','lvlabs','3Dabslvl','3Dlvlabs'}
        lev_mode  = 'row';   abs_mode  = 1;   msg_Ident = 'Abs_BS';
        
    case {'absglb','glbabs','plot','2D','3Dabsglb','3Dglbabs','3Dplot','3D'}
        lev_mode  = 'mat';   abs_mode  = 1;   msg_Ident = 'Abs';
        
    case {'scal','scalCNT'}
        lev_mode  = 'mat';   abs_mode  = 1;   msg_Ident = 'Abs';
        
    otherwise
        plotmode  = 'absglb';
        lev_mode  = 'mat';   abs_mode  = 1;   msg_Ident = 'Abs';
        dim_plot  = '2D';
end
if ~isempty(msg_Ident) , msg_Ident = [msg_Ident '_']; end
msg_Ident = [msg_Ident 'Values_of'];

if abs_mode , dummyCoefs = abs(dummyCoefs); end
if (~isempty(params.xlim) && ~isequal(plotmode,'scal') && ~isequal(plotmode,'scalCNT'))
    xlim = params.xlim;
    if xlim(2)<xlim(1) , xlim = xlim([2 1]); end    
    if xlim(1)<1      , xlim(1) = 1;   end
    if xlim(2)>lenSIG , xlim(2) = lenSIG; end
    indices = xlim(1):xlim(2);
    switch plotmode
      case {'glb','absglb'}
        cmin = min(min(dummyCoefs(:,indices)));
        cmax = max(max(dummyCoefs(:,indices)));
        dummyCoefs(dummyCoefs<cmin) = cmin;
        dummyCoefs(dummyCoefs>cmax) = cmax;

      case {'lvl','abslvl'}
        cmin = min(dummyCoefs(:,indices),[],2);
        cmax = max(dummyCoefs(:,indices),[],2);
        for k=1:nb_SCALES
            ind = dummyCoefs(k,:)<cmin(k);
            dummyCoefs(k,ind) = cmin(k);
            ind = dummyCoefs(k,:)>cmax(k);
            dummyCoefs(k,ind) = cmax(k);
        end
    end
elseif strcmpi(plotmode,'scalCNT')
    if ~isempty(params.nbcl) , nbcl =  params.nbcl; end
end

nb    = min(3,nb_SCALES);
level = '';
for k=1:nb , level = [level ' '  num2str(scales(k))]; end %#ok<AGROW>
if nb<nb_SCALES , level = [level '...']; end
nb     = ceil(nb_SCALES/10);
ytics  = 1:nb:nb_SCALES;
tmp    = scales(1:nb:nb*length(ytics));
ylabs  = num2str(tmp(:),'%0.2f');
plotPARAMS = {NBC,lev_mode,abs_mode,ytics,ylabs,'',xSIG};
  

switch dim_plot
  case 'SC'
      if ~exist('nbcl','var') , nbcl = 10; end
      switch plotmode
          case 'scal',     typePLOT = 'image';
          case 'scalCNT' , typePLOT = 'contour';
      end
      wscalogram(typePLOT,coefs,scales,ySIG,xSIG,'nbcl',nbcl);     
      
      
  case '2D'
    if wtype<5
        titleSTR = getWavMSG(['Wavelet:divCMDLRF:' msg_Ident],level);
        plotPARAMS{6} = titleSTR;
        axeAct = gca;
        plotCOEFS(axeAct,dummyCoefs,plotPARAMS);
    else
        axeAct = subplot(2,2,1);
        titleSTR = getWavMSG('Wavelet:divCMDLRF:Real_Part');
        plotPARAMS{6} = titleSTR;
        plotCOEFS(axeAct,real(dummyCoefs),plotPARAMS);
        axeAct = subplot(2,2,2);
        titleSTR = getWavMSG('Wavelet:divCMDLRF:Imag_Part');
        plotPARAMS{6} = titleSTR;
        plotCOEFS(axeAct,imag(dummyCoefs),plotPARAMS);
        axeAct = subplot(2,2,3);
        titleSTR = getWavMSG('Wavelet:divCMDLRF:Modulus_of');
        plotPARAMS{6} = titleSTR;
        plotCOEFS(axeAct,abs(dummyCoefs),plotPARAMS);
        axeAct = subplot(2,2,4);
        titleSTR = getWavMSG('Wavelet:divCMDLRF:Angle_of');
        plotPARAMS{6} = titleSTR;
        plotCOEFS(axeAct,angle(dummyCoefs),plotPARAMS);
    end
    colormap(pink(NBC));

  case '3D'
    if wtype<5
        titleSTR = getWavMSG(['Wavelet:divCMDLRF:' msg_Ident],level);
        plotPARAMS{6} = titleSTR;
        axeAct = gca;
        surfCOEFS(axeAct,dummyCoefs,plotPARAMS);
    else
        axeAct = subplot(2,2,1);
        titleSTR = 'Real part';
        plotPARAMS{6} = titleSTR;
        surfCOEFS(axeAct,real(dummyCoefs),plotPARAMS);
        axeAct = subplot(2,2,2);
        titleSTR = 'Imaginary part';
        plotPARAMS{6} = titleSTR;
        surfCOEFS(axeAct,imag(dummyCoefs),plotPARAMS);
        axeAct = subplot(2,2,3);
        titleSTR = 'Modulus';
        plotPARAMS{6} = titleSTR;
        surfCOEFS(axeAct,abs(dummyCoefs),plotPARAMS);
        axeAct = subplot(2,2,4);
        titleSTR = 'Angle';
        plotPARAMS{6} = titleSTR;
        surfCOEFS(axeAct,angle(dummyCoefs),plotPARAMS);
    end
    
end        
            
%----------------------------------------------------------------------
function plotCOEFS(axeAct,coefs,plotPARAMS)

[NBC,lev_mode,abs_mode,ytics,ylabs,titleSTR] = deal(plotPARAMS{1:6});

coefs = wcodemat(coefs,NBC,lev_mode,abs_mode);
image(coefs);
set(axeAct, ...
        'YTick',ytics, ...
        'YTickLabel',ylabs, ...
        'YDir','normal', ...
        'Box','On' ...
        );
title(titleSTR,'Parent',axeAct,'fontsize',9);
xlabel(getWavMSG('Wavelet:divCMDLRF:TimeORSpace'),'Parent',axeAct);
ylabel(getWavMSG('Wavelet:divCMDLRF:Scales_a'),'Parent',axeAct);
%----------------------------------------------------------------------
function surfCOEFS(axeAct,coefs,plotPARAMS)

[NBC,~,~,ytics,ylabs,titleSTR] = deal(plotPARAMS{1:6});

surf(coefs);
set(axeAct, ...
        'YTick',ytics, ...
        'YTickLabel',ylabs, ...
        'YDir','normal', ...
        'Box','On' ...
        );
title(titleSTR,'Parent',axeAct);
xlabel(getWavMSG('Wavelet:divCMDLRF:TimeORSpace'),'Parent',axeAct);
ylabel(getWavMSG('Wavelet:divCMDLRF:Scales_a'),'Parent',axeAct);
zlabel('COEFS','Parent',axeAct);

xl = [1 size(coefs,2)];
yl = [1 size(coefs,1)];
zl = [min(min(coefs)) max(max(coefs))];
set(axeAct,'XLim',xl,'YLim',yl,'ZLim',zl,'view',[-30 40]);

colormap(pink(NBC));
shading('interp')
%----------------------------------------------------------------------
function params = parseinputs(lenSIG,stepSIG,varargin)
% Set defaults and fill based on varargin
params.samplingPeriod = [];
params.samplingInterval = [];
params.plotmode = [];
params.nbcl = [];
params.xlim = [];

% Determine if plotmode is specified
tfplotmode = cell2mat(cellfun(@ischar,varargin,'uni',0));
    if any(tfplotmode)
        params.plotmode = varargin{tfplotmode>0};
        idxplotmode = find(tfplotmode);
    end

% Record whether the plotmode is for the scalogram
scalplot = strncmpi(params.plotmode,'scal',1);

%Find how many scalars are specified
tfscalar = cell2mat(cellfun(@isscalar,varargin,'uni',0));
idxscalar = find(tfscalar);
numscalar = nnz(tfscalar);

% Handle case of 1 or 2 scalars and whether or not the scalar occurs
% before or after a string for the scalogram

    if ((numscalar>1 && ~scalplot) || numscalar>2)
        error(message('Wavelet:FunctionInput:Invalid_ScalNum'));
    end
    
    if (numscalar == 1 && ~scalplot)
        params.samplingPeriod = varargin{idxscalar};
    elseif (numscalar == 2 && scalplot)
        params.samplingPeriod = varargin{idxscalar(1)};
        params.nbcl = varargin{idxscalar(end)};
    elseif (numscalar == 1 && scalplot)
        if idxscalar>idxplotmode
            params.nbcl = varargin{idxscalar};
        else
            params.samplingPeriod = varargin{idxscalar};
        end
    end

% Remove any scalar inputs to make parsing of XLIM easier    
varargin(idxscalar) = [];
       
% Determine whether a sampling period is specified elsewhere        
    if (numel(stepSIG) == 1) && (stepSIG ~=1)
        params.samplingInterval = stepSIG;
    end

% Cannot have two different sampling intervals specified    
tfSampInterval = ~isempty(params.samplingInterval);    
    if (tfSampInterval && ~isempty(params.samplingPeriod))
        if (abs(params.samplingInterval-params.samplingPeriod)>sqrt(eps))
            error(message('Wavelet:FunctionInput:SampPeriodsIncompatible'));
        end
    end
    
    if tfSampInterval
        params.samplingPeriod = params.samplingInterval;
    end
% Determine if x-limits for the plot have been specified    
tfxlim = cell2mat(cellfun(@isnumeric, varargin,'uni',0));
    if any(tfxlim)
        params.xlim = varargin{tfxlim>0};
        if (params.xlim(1) < 1 || params.xlim(2) > lenSIG || ...
                ~all(params.xlim == round(params.xlim)))
            error(message('Wavelet:FunctionInput:InvalidPlotlimits'));
        end
    end

    
 

