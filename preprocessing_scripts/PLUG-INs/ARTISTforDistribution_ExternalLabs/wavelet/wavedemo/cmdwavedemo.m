function cmdwavedemo(option,varargin)
%CMDWAVEDEMO Command options for WAVEDEMO function.

%   M. Misiti, Y. Misiti, G. Oppenheim, J.M. Poggi 25-Aug-2012.
%   Last Revision: 26-Aug-2014.
%   Copyright 1995-2014 The MathWorks, Inc.
%   $Revision:  $  $Date:  $

tag_dem_tool  = 'Demo_Tool';
tag_btn_close = 'Demo_Close';       % tag for sub_window
tag_guim_tool = 'Guim_Tool';
tag_sce_tool  = 'dscedw1d';
tag_btn_closeBIS = 'Pus_Close_Win'; % tag for sub_window
tag_sub_close = 'Guim_Close';       % tag for sub_window

switch option
    case 'close'
        %***********************************************%
        %** OPTION = 'close' - close wavedemo window  **%
        %***********************************************%
        mousefrm(0,'watch')
        win = wfindobj('figure','Tag',tag_dem_tool);

        % Closing all opened main analysis windows.
        %------------------------------------------
        pus_handles = wfindobj(0,'Style','pushbutton');
        hdls        = findobj(pus_handles,'flat','Tag',tag_btn_close);
        for i=1:length(hdls)
            try eval(get(hdls(i),'Callback')); end
        end
        try %#ok<*TRYNC>
            hdls = findobj(pus_handles,'flat','Tag',tag_btn_closeBIS);
            for k=1:length(hdls)
                try 
                    FigChild = get(hdls(k),'Parent');
                    eval(get(FigChild,'CloseRequestFcn'));
                end
            end
        end
        try
            WfigPROP = wtbxappdata('get',win,'WfigPROP');
            FigChild = WfigPROP.FigChild;
            for k=1:length(FigChild)
                try eval(get(FigChild,'CloseRequestFcn')); end
            end
        end
        win_gui = wfindobj('figure','Tag',tag_guim_tool);
        try cmdwavedemo('closeguidem',win_gui); end 
        win_sce = wfindobj('figure','Tag',tag_sce_tool);
        try delete(win_sce); end
       
        % Closing the wavedemo window.
        %-----------------------------
        try delete(win); end
        mextglob('clear');
        wtbxmngr('clear');
        mousefrm(0,'arrow');
        
    case 'closeguidem'
        %***********************************************%
        %** OPTION = 'close' - close demoguim window  **%
        %***********************************************%
        mousefrm(0,'watch')

        % Closing all opened main analysis windows.
        %------------------------------------------
        pus_handles = wfindobj(0,'Style','pushbutton');
        hdls        = findobj(pus_handles,'flat','Tag',tag_sub_close);
        for i=1:length(hdls)
            hdl = hdls(i);
            try par = get(hdl,'Parent'); end             %#ok<*TRYNC>
            try eval(get(hdl,'Callback')); end
            try delete(par); end
        end

        % Closing the demoguim window.
        %-----------------------------
        try delete(varargin{1}); end

        win_ini = wfindobj('figure','Tag',tag_dem_tool);
        if ~isempty(win_ini)
            pus_handles = findobj(win_ini,'Style','pushbutton');
            set(pus_handles,'Enable','on');
        else
            mextglob('clear')
            wtbxmngr('clear')
        end
        mousefrm(0,'arrow');
        
end