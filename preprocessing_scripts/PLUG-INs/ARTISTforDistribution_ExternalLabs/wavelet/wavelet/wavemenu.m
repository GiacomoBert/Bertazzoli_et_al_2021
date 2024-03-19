function varargout = wavemenu(varargin)
% WAVEMENU Start the Wavelet Toolbox graphical user interface tools.
%    WAVEMENU launches a menu for accessing the various 
%    graphical tools provided in the Wavelet Toolbox.
%
%    In addition, WAVEMENU(COLOR) let you choose the color
%    preferences. Available values for COLOR are:
%        'k', 'w' , 'y' , 'r' , 'g', 'b' , 'std' (or 's')
%        and 'default' (or 'd').
%
%    WAVEMENU is equivalent to WAVEMENU('default')

% Last Modified by GUIDE v2.5 20-Nov-2010 15:06:19
%
%   M. Misiti, Y. Misiti, G. Oppenheim, J.M. Poggi 12-Mar-96.
%   Last Revision: 07-Aug-2013.
%   Copyright 1995-2013 The MathWorks, Inc.
%   $Revision: 1.22.4.21 $ $Date: 2013/08/23 23:45:16 $


%*************************************************************************%
%                BEGIN initialization code - DO NOT EDIT                  %
%                ----------------------------------------                 %
%*************************************************************************%
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @wavemenu_OpeningFcn, ...
                   'gui_OutputFcn',  @wavemenu_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end
%----------------- Begin solution for G575274 --------------------%
A = wtbxmngr('is_on');
F = wfindobj('Type','figure');
if A
    OK = isempty(F);
    if ~OK
        OK = true;
        for k = 1:length(F)
            if isequal('wavemenu_Win',get(F(k),'Tag'));
                OK = false;
                break
            end
        end
    end
    if OK
        if isappdata(0,'Def_WGlob_Struct'),rmappdata(0,'Def_WGlob_Struct'); end
        if isappdata(0,'Wavelets_Info'),rmappdata(0,'Wavelets_Info'); end 
        if isappdata(0,'WTBX_Glob_Info'),rmappdata(0,'WTBX_Glob_Info'); end 
        if isappdata(0,'DWT_Attribute'),rmappdata(0,'DWT_Attribute'); end 
    end
end
%----------------- End of solution for G575274 -------------------%
if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
%*************************************************************************%
%                END initialization code - DO NOT EDIT                    %
%*************************************************************************%


%*************************************************************************%
%                BEGIN Opening Function                                   %
%                ----------------------                                   %
% --- Executes just before wavemenu is made visible.                      %
%*************************************************************************%
function wavemenu_OpeningFcn(hObject,eventdata,handles,varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to wavemenu (see VARARGIN)

% Choose default command line output for wavemenu
handles.output = hObject;
set(hObject,'WindowStyle','normal');

% Update handles structure
guidata(hObject, handles);

%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%
% TOOL INITIALISATION Introduced manually in the automatic generated code %
%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%
Init_Tool(hObject,eventdata,handles,varargin{:});

% Force translations because of g290327,318823.
 wtranslate(mfilename,hObject)
%*************************************************************************%
%                END Opening Function                                     %
%*************************************************************************%


%*************************************************************************%
%                BEGIN Output Function                                    %
%                ---------------------                                    %
% --- Outputs from this function are returned to the command line.        %
%*************************************************************************%
function varargout = wavemenu_OutputFcn(hObject,eventdata,handles) %#ok<INUSL>
% varargout  cell array for returning output args (see VARARGOUT);

% Get default command line output from handles structure
varargout{1} = handles.output;
%*************************************************************************%
%                END Output Function                                      %
%*************************************************************************%


%=========================================================================%
%                BEGIN Callback Functions                                 %
%                ------------------------                                 %
%=========================================================================%
%-------------------------------------------------------------------
function Pus_Btn_CreateFcn(hObject,eventdata,handles) %#ok<INUSD,DEFNU>
% if isunix
%     bkColor = mextglob('get','Def_UICBkColor');
%     set(hObject,bkColor);
% end
%-------------------------------------------------------------------
function Pus_TOOL_Callback(hObject,eventdata,handles,ToolName) %#ok<INUSL,DEFNU>

mousefrm(0,'watch');
switch ToolName
    case 'dw1dtool'    , dw1dtool;
    case 'wp1dtool'    , wp1dtool;
    case 'cw1dtool'    , cw1dtool;
    case 'cwimtool'    , cwimtool;
    case 'dw2dtool'    , dw2dtool;
    case 'wp2dtool'    , wp2dtool;
    case 'wvdtool'     , wvdtool;
    case 'wpdtool'     , wpdtool;
    case 'sw1dtool'    , sw1dtool;
    case 'de1dtool'    , de1dtool;
    case 're1dtool'    , re1dtool;
    case 'cf1dtool'    , cf1dtool;
    case 'sw2dtool'    , sw2dtool;
    case 'cf2dtool'    , cf2dtool;
    case 'sigxtool'    , sigxtool;
    case 'imgxtool'    , imgxtool;
    case 'wfbmtool'    , wfbmtool;
    case 'wfustool'    , wfustool;
    case 'nwavtool'    , nwavtool;
    case 'wlifttool'   , wlifttool;
    case 'wmspcatool'  , wmspcatool;
    case 'wmuldentool' , wmuldentool;        
    case 'mdw1dtool'   , mdw1dtool;
    case 'wc2dtool'    , wc2dtool;
    case 'dw3dtool'    , dw3dtool;
    case 'cwtfttool'   , cwtfttool; 
    case 'wmp1dtool'   , wmp1dtool;
    case 'cwtfttool2'  , cwtfttool2;
end
mousefrm(0,'arrow');
pause(0.1);
%-------------------------------------------------------------------
function Pus_Close_Win_Callback(hObject,eventdata,handles) %#ok<INUSD>

% Closing all opened main analysis windows.
%------------------------------------------
fig = gcbf;
wfigmngr('close',fig);

% Closing the wavemenu window.
%-----------------------------
try
    delete_Callback;
    delete(fig);
catch ME    %#ok<NASGU>
end
%-------------------------------------------------------------------
function delete_Callback(hObject,eventdata,handles) %#ok<INUSD>

mextglob('clear');
wtbxmngr('clear');
mousefrm(0,'arrow');
%=========================================================================%
%                END Callback Functions                                   %
%=========================================================================%


%=========================================================================%
%                BEGIN Tool Initialization                                %
%                -------------------------                                %
%=========================================================================%
function Init_Tool(hObject,eventdata,handles,varargin) %#ok<INUSL>

% Check for first call.
%----------------------
LstMenusInFig  = findall(get(hObject,'Children'),'flat','Type','uimenu');
lstTagsInFig = get(LstMenusInFig,'tag');
idxMenuFile = find(strcmp('figMenuFile',lstTagsInFig),1);
extendFLAG = isempty(idxMenuFile);

nbIN = length(varargin);
switch nbIN
    case 0     , varargin{1} = [];
    case {1,2,3} 
    otherwise
        error(message('Wavelet:FunctionInput:TooMany_ArgNum'))
end

if ~wtbxmngr('is_on') , wtbxmngr('ini'); end
first = ~mextglob('is_on');
if first 
    mextglob('ini',varargin{:});
elseif ~isempty(varargin{1})
    mextglob('pref',varargin{:});    
else
    return
end

if extendFLAG
    wfigmngr('extfig',hObject,'ExtMainFig_WTBX');
end

% Set CLOSE functions.
%---------------------
set(hObject,'CloseRequestFcn',@Pus_Close_Win_Callback)
MenusInFig  = findall(hObject,'Type','uimenu');
% LabelsInFig = get(MenusInFig,'label');
% TMP = strfind(LabelsInFig,'Close');
% idxMenuClose = cat(2,TMP{:});
TagsInFig = get(MenusInFig,'tag');
idxMenuClose =  find(strcmp(TagsInFig,'figMenuClose'));
if ~isempty(idxMenuClose)
    hMenu_Close = MenusInFig(idxMenuClose);
    set(hMenu_Close,'Callback',@Pus_Close_Win_Callback)
end

% Set colors and fontes for the figure.
%---------------------------------------
wfigmngr('set_FigATTRB',hObject,'wavemenu');

if extendFLAG
    redimfigATTRB = wtbxappdata('get',hObject,'redimfigATTRB');
    if isempty(redimfigATTRB)
        redimfig('On',hObject,[0.8 1.2],'left');
        wtbxappdata('set',hObject,'redimfigATTRB',true);
    end
end

set(hObject,'DeleteFcn',@delete_Callback);
%=========================================================================%
%                END Tool Initialization                                  %
%=========================================================================%
