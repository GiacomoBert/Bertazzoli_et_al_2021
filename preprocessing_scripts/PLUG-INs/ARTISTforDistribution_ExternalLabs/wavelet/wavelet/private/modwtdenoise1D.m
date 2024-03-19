function [xd,cxd] = modwtdenoise1D(x,wav,lev,softHard,scal)
%   [xd,cxd] = modwtdenoise1D(x,wav,lev,softHard,scal) returns denoised
%   1D signal xd and MODWT coefficients cxd.   

 % Only support 'mln' at present
    if ~strcmpi(scal,'mln')
        error(message('Wavelet:modwt:ScalingError'));
    end

    if (isrow(x) || iscolumn(x))    
    % Obtain the MODWT
        wt = modwt(x,wav,lev);
    else
        wt = x;
    end
    
validateattributes(wt,{'double'},{'2d','real'});

% Check to see that the level for denoising does not exceed the level
% of the transform
    if (lev >= size(wt,1))
        error(message('Wavelet:modwt:InvalidDenoiseLevel'));
    end    

% Determine the level dependent thresholds
    for kk = 1:lev
        madest(kk) = sqrt(2)*median(abs(wt(kk,:)))/0.6745;
    end
% Thresholds    
thr = sqrt(2*madest.^2./2.^(1:lev)*log(length(x)));

% Threshold MODWT coefficients
    for kk = 1:lev
        wt(kk,:) = wthresh(wt(kk,:),softHard,thr(kk));
    end
% Invert the MODWT    
xd = imodwt(wt,wav);
cxd = wt;