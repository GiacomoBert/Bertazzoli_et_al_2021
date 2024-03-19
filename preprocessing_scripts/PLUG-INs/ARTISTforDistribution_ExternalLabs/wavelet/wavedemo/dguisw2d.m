function varargout = dguisw2d(varargin)
%DGUISW2D Shows 2-D SWT de-noising wavelet GUI tools in the Wavelet Toolbox.
%
% This is a slideshow file for use with wshowdrv.m
% To see it run, type 'wshowdrv dguisw2d',
%
%   See also SWT2, ISWT2.

%   M. Misiti, Y. Misiti, G. Oppenheim, J.M. Poggi 05-May-99.
%   Last Revision: 12-Apr-2012.
%   Copyright 1995-2012 The MathWorks, Inc.

% Initialization and Local functions if necessary.
if nargin>0
    action = varargin{1};
    switch action
        case 'auto'    , wshowdrv('#autoMode',mfilename,'close');
        case 'gr_auto' , wshowdrv('#gr_autoMode',mfilename,'close');

        case 'getFigParam'
            figName  = 'Wavelet GUI Example: SW2D';
            showType = 'command';
            varargout = {figName,showType};

        case 'slidePROC_Init'
            figHandle = varargin{2};
            localPARAM = wtbxappdata('get',figHandle,'localPARAM');
            if ~isempty(localPARAM)
                active_fig = localPARAM{1};
                delete(active_fig);
                wtbxappdata('del',figHandle,'localPARAM');
            end

        case 'slidePROC'
            [figHandle,idxSlide]  = deal(varargin{2:end});
            localPARAM = wtbxappdata('get',figHandle,'localPARAM');
            if isempty(localPARAM)
                active_fig = sw2dtool;
                wenamngr('Inactive',active_fig);
                localPARAM = {active_fig};
                wtbxappdata('set',figHandle,'localPARAM',localPARAM);
                wshowdrv('#modify_cbClose',figHandle,active_fig,'sw2dtool');
            else
                active_fig = deal(localPARAM{:});
            end
            numDEM = idxSlide-1;
            demoSET = {...
                'noiswom'      , 'haar', 3 , {'penallo',46.12} , 'BW' ; ...
                'noiswom'      , 'haar', 5 , {'penallo',48.62} , 'BW' ; ...
                'noiswom'      , 'db3' , 4 , {'penallo',NaN}   , 'BW'; ...
                'nbarb1'       , 'db1' , 4 , {} ,                'BW'; ...
                'noissi2d'     , 'db1' , 2 , {} ,                'BW'; ...
                'jellyfish256' , 'db1' , 3 , {'penalhi',38} ,    'COL' ...
                };
            nbDEM = size(demoSET,1);
            if ismember(numDEM,1:nbDEM)
                paramDEM = demoSET(numDEM,:);
                sw2dtool('demo',active_fig,paramDEM{:});
                wenamngr('Inactive',active_fig);
            end
    end
    return
end

if nargout<1,
    wshowdrv(mfilename)
else
    idx = 0;	slide(1).code = {}; slide(1).text = {};

    %========== Slide 1 ==========
    idx = idx+1;
    slide(idx).code = {
        'figHandle = gcf;'
        [mfilename ,'(''slidePROC_Init'',figHandle);']
        '' };

    %========== Slide 2 to 7 ==========
    for idx = 2:7
        slide(idx).code = {
            [mfilename ,'(''slidePROC'',figHandle,', int2str(idx), ');']
            };
    end

    varargout{1} = slide;

end

