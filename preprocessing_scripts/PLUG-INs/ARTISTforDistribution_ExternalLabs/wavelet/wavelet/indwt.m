function x = indwt(varargin)
%INDWT Inverse nondecimated 1-D wavelet transform.
%   INDWT will be removed in a future release of MATLAB. Use one of the
%   following functions instead:
%       <a href="matlab:help imodwt">imodwt</a>
%       <a href="matlab:help iswt">iswt</a>

% Error in R2015a
error(message('Wavelet:warnobsolete:ErrorReplaceINDWT'));
nbIn = nargin;
narginchk(1,6);
if ~isstruct(varargin{1})
    error(message('Wavelet:FunctionArgVal:Invalid_ArgVal'));
end

wt = varargin{1};
ndCell = wt.dec;
LoR = wt.filters.LoR;
HiR = wt.filters.HiR;
LX  = wt.longs(end);
level = wt.level;
nextArg = 2;

nam_Rec = 's';
lev_Rec = 0;
while nbIn>=nextArg
    argName = lower(varargin{nextArg});
    argVal  = varargin{nextArg+1};
    nextArg = nextArg + 2;
    switch argName
        case {'a','d','ca','cd'}
            nam_Rec = argName;
            lev_Rec = argVal;
    end
end

fistIDX  = 1;
lastSTEP = level;
switch nam_Rec
    case 's'
        
    case 'cd'
        num2keep = level+1-(lev_Rec-1);
        x = ndCell{num2keep};
        return;
        
    case 'ca'
        if lev_Rec==level
             x = ndCell{1};
             return
        end
        for k = level+1:-1:2+(level-lev_Rec)
            ndCell{k} = zeros(size(ndCell{k}));
        end
        lastSTEP = level-lev_Rec;
        
    case 'd'
        set2zero = (1:level+1);
        num2keep = level+2-lev_Rec;
        set2zero = setdiff(set2zero,num2keep);
        for k = set2zero
            ndCell{k} = zeros(size(ndCell{k}));
        end        
        
    case 'a'
        for k = level+1:-1:2+(level-lev_Rec)
            ndCell{k} = zeros(size(ndCell{k}));
        end
end

idx = fistIDX;
for k=1:lastSTEP
    a = conv(ndCell{idx},LoR);
    d = conv(ndCell{idx+1},HiR);
    ndCell{idx+1} = (a+d)/2;
    if idx<level
        ndCell{idx+1} = wkeep(ndCell{idx+1},length(ndCell{idx+2}),'c');
    end
    idx = idx+1;
end
if isequal(nam_Rec,'ca') && ~isequal(lev_Rec,0)
    x = ndCell{idx}; return; 
end
x = wkeep(ndCell{end},LX,'c');

