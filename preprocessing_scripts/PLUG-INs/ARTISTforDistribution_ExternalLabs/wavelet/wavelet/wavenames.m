function W = wavenames(T)
%WAVENAMES Wavelet names information.
%   W = WAVENAMES(T) returns a cell array which contains
%   the name of all wavelets of type T. The valid values
%   for T are:
%       - 'all'  : all wavelets.
%       - 'lazy' : the "lazy" wavelet.
%       - 'orth' : orthogonal wavelets.
%       - 'bior' : biorthogonal wavelets.
%
%    W = WAVENAMES is equivalent to W = WAVENAMES('all').

%   M. Misiti, Y. Misiti, G. Oppenheim, J.M. Poggi 21-Jun-2003.
%   Last Revision: 20-Dec-2010.
%   Copyright 1995-2010 The MathWorks, Inc.

wnameCell_LAZY = {'lazy'};

wnameCell_ORTH = {...
    'haar', ...
    'db1','db2','db3','db4','db5','db6','db7','db8','db9','db10' , ...
    'sym2','sym3','sym4','sym5','sym6','sym7','sym8',...
    'coif1','coif2','coif3','coif4','coif5',  ...
    };

wnameCell_BIOR = {...
    'bior1.1', 'bior1.3', 'bior1.5', ...
    'bior2.2', 'bior2.4', 'bior2.6', 'bior2.8', ...
    'bior3.1', 'bior3.3', 'bior3.5', 'bior3.7', 'bior3.9',  ...
    'bior4.4', 'bior5.5', 'bior6.8' , ...
    'rbio1.1', 'rbio1.3', 'rbio1.5', ...
    'rbio2.2', 'rbio2.4', 'rbio2.6', 'rbio2.8', ...
    'rbio3.1', 'rbio3.3', 'rbio3.5', 'rbio3.7', 'rbio3.9', ...
    'rbio4.4', 'rbio5.5', 'rbio6.8', ...
    'cdf1.1','cdf1.3','cdf1.5', ...
    'cdf2.2','cdf2.4','cdf2.6', ...
    'cdf3.1','cdf3.3','cdf3.5', ...
    'cdf4.2','cdf4.4','cdf4.6', ...
    'cdf5.1','cdf5.3','cdf5.5', ...
    'cdf6.2','cdf6.4','cdf6.6', ...
    'bs3','9.7','rbs3','r9.7' , ...
    };

if nargin<1 , T = 'all'; end
switch lower(T)
    case 'lazy' , W = wnameCell_LAZY';
    case 'orth' , W = wnameCell_ORTH';
    case 'bior' , W = wnameCell_BIOR';
    case 'all'  , 
        W = [wnameCell_LAZY,wnameCell_ORTH,wnameCell_BIOR]';
    otherwise 
        error(message('Wavelet:FunctionArgVal:Unknown_TypWav'));
end
