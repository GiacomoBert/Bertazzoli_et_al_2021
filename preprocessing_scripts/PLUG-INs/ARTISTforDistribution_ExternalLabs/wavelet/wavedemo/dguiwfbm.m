function varargout = dguiwfbm(varargin)
%DGUIWFBM Shows Fractional Brownian motion generation tool in the Wavelet Toolbox.
%
% This is a slideshow file for use with wshowdrv.m
% To see it run, type 'wshowdrv dguiwfbm', 

%   M. Misiti, Y. Misiti, G. Oppenheim, J.M. Poggi 10-Sep-2003.
%   Last Revision: 23-May-2012.
%   Copyright 1995-2012 The MathWorks, Inc.

% Initialization and Local functions if necessary.
if nargin>0
	action = varargin{1};
	switch action
	  case 'auto'    , wshowdrv('#autoMode',mfilename,'close');
	  case 'gr_auto' , wshowdrv('#gr_autoMode',mfilename,'close');

	  case 'getFigParam'
		  figName  = getWavMSG('Wavelet:wavedemoMSGRF:DGUI_WFBM');
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
			  active_fig = wfbmtool;
			  wenamngr('Inactive',active_fig);
			  localPARAM = {active_fig};
			  wtbxappdata('set',figHandle,'localPARAM',localPARAM);
			  wshowdrv('#modify_cbClose_NEW',figHandle,active_fig,'wfbmtool');
		  else
			  active_fig = deal(localPARAM{:});
		  end
		  numDEM = idxSlide-1;
          demoSET = {...
                  0.2 , 'generate'    ; ...
                  0.2 , 'statistics'  ; ...
                  0.5 , 'generate'    ; ...
                  0.5 , 'statistics'  ; ...
                  0.9 , 'generate'    ; ...
                  0.9 , 'statistics'    ...                  
              };
          nbDEM = size(demoSET,1);
          if ~ismember(numDEM,(1:nbDEM)) , return; end
          paramDEM = demoSET(numDEM,:);
          wfbmtool('demoPROC',active_fig,[],[],paramDEM);
          wfigmngr('storeValue',active_fig,'File_Save_Flag',1);
          wenamngr('Inactive',active_fig);
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
      if any(find(idx==[4 6]))
          slide(idx).idxPrev = idx-2;
      end
  end
  
  varargout{1} = slide;
  
end

