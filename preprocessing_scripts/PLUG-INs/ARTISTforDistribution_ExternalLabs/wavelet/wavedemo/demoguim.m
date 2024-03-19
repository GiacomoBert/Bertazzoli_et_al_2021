function out1 = demoguim(option,in2,in3) %#ok<INUSD>
%DEMOGUIM Main GUI-mode examples menu in the Wavelet Toolbox.
%   DEMOGUIM creates the window for GUI mode examples.

%   DEMOGUIM('auto') shows all the GUI mode examples
%   in automatic mode.
%
%   DEMOGUIM('gr_auto') shows all the GUI mode examples
%   in automatic mode: first in the increasing slide order
%   and then in the decreasing slide order.
%
%   DEMOGUIM('loop') shows all the GUI mode examples
%   in loop mode.

%   M. Misiti, Y. Misiti, G. Oppenheim, J.M. Poggi 12-Mar-96.
%   Last Revision: 26-Aug-2014.
%   Copyright 1995-2014 The MathWorks, Inc.
%   $Revision: 1.17.4.8.2.4 $ $Date: 2011/04/03 13:38:05 $

if nargin==0 , option = 'create'; end

tag_guim_tool = 'Guim_Tool';
tag_sub_close = 'Guim_Close';  % tag for sub_window

switch option
    case 'create'
        win = demoguimwin;
        set(win,'CloseRequestFcn',@(o,e)cmdwavedemo ('closeguidem',win))
        if nargout>0 , out1 = win; end

    case {'dw1d','wp1d','cw1d','cwim','dw2d','wp2d','wvdi','wpdi',...
          'sw1d','de1d','re1d','cf1d','sw2d','cf2d','sext','iext',...
          'nwav','wfbm','wfus','mdw1d','wmul','wmspca'}
        if nargin==1
            win = wfindobj('figure','Tag',tag_guim_tool);
        else
            win = in2;
        end
        switch option
            case {'futureDEMOS'} % Template
                uiwait(msgbox('This example is not available.',...
                    'Caution','warn','modal'));
                return
        end
        if ishandle(win)
            set(wfindobj(win,'Style','pushbutton'),'Enable','Off');
        end
        demoname = ['dgui' option];
		addclose = [mfilename '(''endshow'');'];
		fig = wshowdrv(demoname,addclose,tag_sub_close);
		if nargout>0 , out1 = fig; end

    case {'auto','gr_auto'}
        %**************************************%
        %** OPTION = 'auto' and 'gr_auto'    **%
		%** All examples in automatic modes. **%
        %**************************************%
		lstDEMOS = {...
				'dw1d','wp1d','cw1d','cwim','dw2d','wp2d','wvdi','wpdi', ...
				'sw1d','de1d','re1d','cf1d','sw2d','cf2d','sext','iext',...
                'nwav','wfbm','wfus','mdw1d','wmul','wmspca'};
		win = demoguim('create');
		set(wfindobj(win,'Style','pushbutton'),'Enable','Off');
		stop = 0;
        while stop==0
            for k=1:length(lstDEMOS)
                feval(['dgui',lstDEMOS{k}],option);
                if ~ishandle(win) , stop = 1; break; end
            end
            if ~isequal(stop,1) && nargin==2 && isequal(in2,'loop')
                stop = 0 ;
            else
                stop = 1;
            end
        end
        cmdwavedemo('closeguidem',win)

    case 'loop'
        %**************************************************************%
        %** OPTION = 'loop' - loop for all examples (automatic mode) **%
        %**************************************************************%
        demoguim('auto','loop');

    case 'endshow'
        %**************************************************%
        %** OPTION = 'endshow' - used to finish wshowdrv **%
        %**************************************************%
        win = wfindobj('figure','Tag',tag_guim_tool);
        set(findobj(win,'Style','pushbutton'),'Enable','on');
		mousefrm(win,'arrow');
        drawnow;

    case 'close'
        %***********************************************%
        %** OPTION = 'close' - close demoguim window  **%
        %***********************************************%
        win = wfindobj('figure','Tag',tag_guim_tool);
        cmdwavedemo('closeguidem',win)

    otherwise
        errargt(mfilename,getWavMSG('Wavelet:moreMSGRF:Unknown_Opt'),'msg');
        error(message('Wavelet:FunctionArgVal:Invalid_ArgVal'));
end

