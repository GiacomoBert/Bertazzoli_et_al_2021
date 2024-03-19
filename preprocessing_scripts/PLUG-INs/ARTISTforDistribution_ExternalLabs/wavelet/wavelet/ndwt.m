function varargout = ndwt(x,level,varargin)
%NDWT Nondecimated 1-D wavelet transform.
%   NDWT will be removed in a future release of MATLAB. Use one of the
%   following functions instead:
%       <a href="matlab:help modwt">modwt</a>
%       <a href="matlab:help swt">swt</a>

% Error in R2015a

error(message('Wavelet:warnobsolete:ErrorReplaceNDWT'));
nbIn = length(varargin);
if nbIn < 1
    error(message('MATLAB:narginchk:notEnoughInputs'));
elseif nbIn > 3
    error(message('MATLAB:narginchk:tooManyInputs'));
end
nextArg = 2;
if ischar(varargin{1})     % Wavelet name
    [LoD,HiD,LoR,HiR] = wfilters(varargin{1});
    
elseif iscell(varargin{1}) % Wavelet Filters
    [LoD,HiD,LoR,HiR] = deal(varargin{1}{:});
else                       % Wavelet Filters
    if isfield(varargin{1},'LoD') && isfield(varargin{1},'HiD') && ...
            isfield(varargin{1},'LoR') && isfield(varargin{1},'HiR')
        LoD = varargin{1}.LoD; HiD = varargin{1}.HiD;
        LoR = varargin{1}.LoR; HiR = varargin{1}.HiR;
    else
        error(message('Wavelet:FunctionArgVal:Invalid_ArgVal'));
    end
end
lf = length(LoD);

dwtEXTM = 'sym';
while nbIn>=nextArg
    argName = varargin{nextArg};
    argVal  = varargin{nextArg+1};
    nextArg = nextArg + 2;
    switch argName
        case 'mode' , dwtEXTM = argVal;
    end
end

% Initialization.
rowvect =  size(x,1)<=1;
lx = length(x);
longs = zeros(1,level+2);
longs(end) = lx;
dec = cell(1,level+1);
% idx = 1;
idx = level+1;
for k=1:level
    lx = length(x);
    x = wextend('1d',dwtEXTM,x,lf-1,'b');
    lkeep = lx+lf-1;
    dec{idx} = wkeep(conv(HiD,x),lkeep);
    x = wkeep(conv(LoD,x),lkeep);
    idx = idx-1;
end
dec{idx} = x;
for k = 1:level+1 , longs(k) = length(dec{k}); end

% Non Decimated Wavelet Transform.
wt.rowvect = rowvect;
wt.level = level;
wt.mode  = dwtEXTM;
wt.filters.LoD = LoD;
wt.filters.HiD = HiD;
wt.filters.LoR = LoR;
wt.filters.HiR = HiR;
wt.dec = dec;
wt.longs = longs;

switch nargout
    case 1
        varargout{1} = wt;
    case 2        
        if ~rowvect , catDIR = 1; longs = longs'; else catDIR = 2; end
        varargout{1} = cat(catDIR,dec{:});
        varargout{2} = longs;
    case 3        
        NbCd1 = length(wt.dec{1});
        A = wt.dec{end}(end-NbCd1+1:end);
        D = zeros(level,NbCd1);
        for k = 1:level
            D(k,:) = wt.dec{k}(end-NbCd1+1:end);
        end
        varargout = {A,D,wt};        
end
