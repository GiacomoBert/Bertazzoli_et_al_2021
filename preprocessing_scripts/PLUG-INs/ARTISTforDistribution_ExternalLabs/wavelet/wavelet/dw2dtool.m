function varargout = dw2dtool(option,varargin)
%DW2DTOOL Discrete wavelet 2-D tool.
%   VARARGOUT = DW2DTOOL(OPTION,VARARGIN)

%   M. Misiti, Y. Misiti, G. Oppenheim, J.M. Poggi 12-Mar-96.
%   Last Revision 08-May-2012.
%   Copyright 1995-2012 The MathWorks, Inc.

% Test inputs.
%-------------
if nargin==0 , option = 'create'; end
[option,winAttrb] = utguidiv('ini',option,varargin{:});

% Default values.
%----------------
max_lev_anal = 8;
default_nbcolors = 128;

% Tag property of objects.
%-------------------------
tag_m_exp_wrks  = 'm_exp_wrks';
tag_pus_anal  = 'Pus_Anal';
tag_pus_deno  = 'Pus_Deno';
tag_pus_comp  = 'Pus_Comp';
tag_pus_hist  = 'Pus_Hist';
tag_pus_stat  = 'Pus_Stat';
tag_pop_declev= 'Pop_DecLev';
tag_pus_vis   = 'Pus_Visu';
tag_pus_big   = 'Pus_Big';
tag_pus_rec   = 'Pus_Rec';
tag_pop_viewm = 'Pop_ViewM';
tag_txt_full  = 'Txt_Full';
tag_pus_full  = ['Pus_Full.1';'Pus_Full.2';'Pus_Full.3';'Pus_Full.4'];

tag_axefigutil = 'Axe_FigUtil';
tag_linetree   = 'Tree_lines';
tag_txttree    = 'Tree_txt';
tag_axeimgbig  = 'Axe_ImgBig';
tag_axeimgini  = 'Axe_ImgIni';
tag_axeimgvis  = 'Axe_ImgVis';
tag_axeimgsel  = 'Axe_ImgSel';
tag_axeimgdec  = 'Axe_ImgDec';
tag_axeimgsyn  = 'Axe_ImgSyn';
% tag_axeimghdls = 'Img_Handles';

% MemBloc0 of stored values.
%---------------------------
n_InfoInit   = 'DW2D_InfoInit';
% ind_filename = 1;
% ind_pathname = 2;
nb0_stored   = 2;

% MemBloc1 of stored values.
%---------------------------
n_param_anal   = 'DWAn2d_Par_Anal';
% ind_img_name   = 1;
% ind_wav_name   = 2;
% ind_lev_anal   = 3;
% ind_img_t_name = 4;
% ind_img_size   = 5;
% ind_nbcolors   = 6;
ind_act_option = 7;
% ind_simg_type  = 8;
% ind_thr_val    = 9;
nb1_stored     = 9;

% MemBloc3 of stored values.
%---------------------------
n_miscella      = 'DWAn2d_Miscella';
ind_graph_area  =  1;
% ind_pos_axebig  =  2;
% ind_pos_axeini  =  3;
% ind_pos_axevis  =  4;
% ind_pos_axedec  =  5;
% ind_pos_axesyn  =  6;
% ind_pos_axesel  =  7;
ind_view_status =  8;
ind_save_status =  9;
ind_sel_funct   = 10;
nb3_stored      = 10;

% Miscellaneous values.
%----------------------
% square_viewm  = 1;
% tree_viewm    = 2;

switch option
    case 'create'

        % Get Globals.
        %--------------
        [Def_Txt_Height,Def_Btn_Height,Def_Btn_Width,Pop_Min_Width,  ...
         X_Spacing,Y_Spacing,Def_FraBkColor] = ...
            mextglob('get',...
              'Def_Txt_Height','Def_Btn_Height','Def_Btn_Width','Pop_Min_Width',  ...
              'X_Spacing','Y_Spacing','Def_FraBkColor' ...
              );

        % Variables initialization.
        %--------------------------
        dw2d_PREFS = wtbutils('dw2d_PREFS');
        Col_BoxTitleSel = dw2d_PREFS.Col_BoxTitleSel;
        linDW2D_Color   = Col_BoxTitleSel;
        select_funct = 'dw2dimgs(''get_img'');';

        % Wavelet 2-D window initialization.
        %-----------------------------------
        win_title = getWavMSG('Wavelet:dw2dRF:NamWinDW2D');
        [win_dw2dtool,pos_win,win_units,str_numwin,...
             pos_frame0,Pos_Graphic_Area] = ...
                wfigmngr('create',win_title,winAttrb,'ExtFig_Tool',mfilename,1,1,0);
        set(win_dw2dtool,'Tag','DW2DTOOL');    
        if nargout>0 , varargout{1} = win_dw2dtool; end
		
		% Add Coloration Mode Submenu.
		%-----------------------------
		wfigmngr('add_CCM_Menu',win_dw2dtool);

		% Add Help for Tool.
		%------------------
		wfighelp('addHelpTool',win_dw2dtool,getWavMSG('Wavelet:dw2dRF:HLP_TwoDim'),'DW2D_GUI');

		% Add Help Item.
		%----------------
		wfighelp('addHelpItem',win_dw2dtool,getWavMSG('Wavelet:dw2dRF:HLP_WkImages'),'DW2D_WORKING');
		wfighelp('addHelpItem',win_dw2dtool,getWavMSG('Wavelet:dw2dRF:HLP_LoadSave'),'DW2D_LOADSAVE');

        % Menu construction for current figure.
        %--------------------------------------
		[m_files,m_load,m_save] = ...
			wfigmngr('getmenus',win_dw2dtool,'file','load','save');
        set(m_save,'Enable','Off');

        m_loadtst = uimenu(m_files,...
            'Label',getWavMSG('Wavelet:commongui:Lab_Example'), ...
            'Position',3,             ...
            'Separator','Off'         ...
            );
        m_imp_wrks = uimenu(m_files,...
            'Label',getWavMSG('Wavelet:commongui:Lab_Import'), ...
            'Position',4,'Separator','On'    ...
            );
        m_exp_wrks = uimenu(m_files,...
            'Label',getWavMSG('Wavelet:commongui:Lab_Export'),  ...
            'Position',5,'Enable','Off','Separator','Off', ...
            'Tag',tag_m_exp_wrks  ...
            );       
        
        uimenu(m_load,...
            'Label',getWavMSG('Wavelet:commongui:Str_Image'),    ...
            'Position',1,             ...
            'Callback',['dw2dmngr(''load_img'',' str_numwin ');']  ...
            );

        uimenu(m_load,...
            'Label',getWavMSG('Wavelet:commongui:Str_Coefficients'), ...
            'Position',2,             ...
            'Callback',['dw2dmngr(''load_cfs'',' str_numwin ');']  ...
                                );
        uimenu(m_load,...
            'Label',getWavMSG('Wavelet:dw2dRF:Lab_Decomposition'), ...
            'Position',3,              ...
            'Callback',['dw2dmngr(''load_dec'','  str_numwin ');']   ...
            );
                            
		uimenu(m_save,...
            'Label',getWavMSG('Wavelet:dw2dRF:Syn_Img'),...
            'Position',1,                 ...
            'Callback',['dw2dmngr(''save_synt'',' str_numwin ');'] ...
            );
        uimenu(m_save,...
            'Label',getWavMSG('Wavelet:commongui:Str_Coefficients'), ...
            'Position',2,             ...
            'Callback',['dw2dmngr(''save_cfs'',' str_numwin ');']  ...
            );
        uimenu(m_save,...
            'Label',getWavMSG('Wavelet:commongui:Str_Decomp'), ...
            'Position',3,              ...
            'Callback',['dw2dmngr(''save_dec'',' str_numwin ');']   ...
            );
        Men_Save_APP = uimenu(m_save,...
            'Label',getWavMSG('Wavelet:dw2dRF:Str_app_sig'), ...
            'Position',4,             ...
            'Separator','On',...
            'Tag','Men_Save_APP' ...
            );
        uimenu(Men_Save_APP,...
            'Label',getWavMSG('Wavelet:dw2dRF:AppLevX',1), ...
            'Position',1, ...
            'Callback',['dw2dmngr(''save_app'',' str_numwin ');']  ...
            );
         Men_Save_APP_CFS = uimenu(m_save,...
            'Label',getWavMSG('Wavelet:dw2dRF:Lab_AppCfs'), ...
            'Position',5,             ...
            'Separator','On',...
            'Tag','Men_Save_APP_CFS' ...
            );
        uimenu(Men_Save_APP_CFS,...
            'Label',getWavMSG('Wavelet:dw2dRF:Lab_AppCfsOf',1), ...
            'Position',1,             ...
            'Callback',['dw2dmngr(''save_app_cfs'',' str_numwin ');']  ...
            );

        % Submenu of test signals.
        %-------------------------
        m_demoIDX = uimenu(m_loadtst, ...
           'Label',getWavMSG('Wavelet:dw2dRF:Lab_IndImg'), ...
           'Tag','Lab_IndImg','Position',1);        
        m_demoCOL = uimenu(m_loadtst, ...
           'Label',getWavMSG('Wavelet:dw2dRF:Lab_ColImg'), ...
           'Tag','Lab_ColImg','Position',2);  
        
        beg_call_str = ['dw2dmngr(''demo'',' str_numwin];
        tab  = char(9);
        lab_str = getWavMSG('Wavelet:dw2dRF:Lab_Example',3,'sym4', ...
            tab,getWavMSG('Wavelet:moreMSGRF:EX2D_Name_detail'));
        cba     = [beg_call_str ',''detail'',''sym4'',3,''BW'');'];
        uimenu(m_demoIDX,'Label',lab_str,'Callback',cba);

        lab_str = getWavMSG('Wavelet:dw2dRF:Lab_Example',2,'haar', ...
            tab,getWavMSG('Wavelet:moreMSGRF:EX2D_Name_woman'));        
        cba     = [beg_call_str ',''woman'',''haar'',2,''BW'');'];
        uimenu(m_demoIDX,'Label',lab_str,'Callback',cba);

        lab_str = getWavMSG('Wavelet:dw2dRF:Lab_Example',3,'haar', ...
            tab,getWavMSG('Wavelet:moreMSGRF:EX2D_Name_detfingr'));        
        cba     = [beg_call_str ',''detfingr'',''haar'',3,''BW'');'];
        uimenu(m_demoIDX,'Label',lab_str,'Callback',cba);

        lab_str = getWavMSG('Wavelet:dw2dRF:Lab_Example',3,'haar', ...
            tab,getWavMSG('Wavelet:moreMSGRF:EX2D_Name_tire'));               
        cba     = [beg_call_str ',''tire'',''haar'',3,''BW'');'];
        uimenu(m_demoIDX,'Label',lab_str,'Callback',cba);
        
        lab_str = getWavMSG('Wavelet:dw2dRF:Lab_Example',3,'sym3', ...
            tab,getWavMSG('Wavelet:moreMSGRF:EX2D_Name_chess'));                       
        cba     = [beg_call_str ',''chess'',''sym3'',3);'];
        uimenu(m_demoIDX,'Label',lab_str,'Callback',cba);

        lab_str = getWavMSG('Wavelet:dw2dRF:Lab_Example',2,'bior3.7', ...
            tab,getWavMSG('Wavelet:moreMSGRF:EX2D_Name_barb'));                       
        cba     = [beg_call_str ',''wbarb'',''bior3.7'',2,''BW'');'];
        uimenu(m_demoIDX,'Label',lab_str,'Callback',cba);

        lab_str = getWavMSG('Wavelet:dw2dRF:Lab_Example',3,'bior5.5', ...
            tab,getWavMSG('Wavelet:moreMSGRF:EX2D_Name_facets'));                               
        cba     = [beg_call_str ',''facets'',''bior5.5'',3);'];
        uimenu(m_demoIDX,'Separator','on','Label',lab_str,'Callback',cba);

        lab_str = getWavMSG('Wavelet:dw2dRF:Lab_Example',3,'bior4.4', ...
            tab,getWavMSG('Wavelet:moreMSGRF:EX2D_Name_geometry'));                               
        cba     = [beg_call_str ',''geometry'',''bior4.4'',3);'];
        uimenu(m_demoIDX,'Label',lab_str,'Callback',cba);

        lab_str = getWavMSG('Wavelet:dw2dRF:Lab_Example',4,'db3', ...
            tab,getWavMSG('Wavelet:moreMSGRF:EX2D_Name_sinsin'));                               
        cba     = [beg_call_str ',''sinsin'',''db3'',4);'];
        uimenu(m_demoIDX,'Label',lab_str,'Callback',cba);

        lab_str = getWavMSG('Wavelet:dw2dRF:Lab_Example',3,'coif2', ...
            tab,getWavMSG('Wavelet:moreMSGRF:EX2D_Name_Tartan'));                               
        cba     = [beg_call_str ',''tartan'',''coif2'',3);'];
        uimenu(m_demoIDX,'Label',lab_str,'Callback',cba);

        lab_str = getWavMSG('Wavelet:dw2dRF:Lab_Example',3,'db1', ...
            tab,getWavMSG('Wavelet:moreMSGRF:EX2D_Name_Mandel'));                               
        cba     = [beg_call_str ',''mandel'',''db1'',3);'];
        uimenu(m_demoIDX,'Separator','on','Label',lab_str,'Callback',cba);

        lab_str = getWavMSG('Wavelet:dw2dRF:Lab_Example',3,'db1', ...
            tab,getWavMSG('Wavelet:moreMSGRF:EX2D_Name_Julia'));                                       
        cba     = [beg_call_str ',''julia'',''db1'',3);'];
        uimenu(m_demoIDX,'Label',lab_str,'Callback',cba);

        lab_str = getWavMSG('Wavelet:dw2dRF:Lab_Example',5,'coif1', ...
            tab,getWavMSG('Wavelet:moreMSGRF:EX2D_Name_Ifs'));                               
        cba     = [beg_call_str ',''wifs'',''coif1'',5);'];
        uimenu(m_demoIDX,'Label',lab_str,'Callback',cba);
        
        lab_str = getWavMSG('Wavelet:dw2dRF:Lab_Example',3,'db2', ...
            tab,getWavMSG('Wavelet:moreMSGRF:EX2D_Name_Nwoman'));                               
        cba     = [beg_call_str ',''noiswom'',''db2'',3,''BW'');'];
        uimenu(m_demoIDX,'Separator','on','Label',lab_str,'Callback',cba);

        lab_str = getWavMSG('Wavelet:dw2dRF:Lab_Example',3,'bior6.8', ...
            tab,getWavMSG('Wavelet:moreMSGRF:EX2D_Name_NSinsin'));                                       
        cba     = [beg_call_str ',''noissi2d'',''bior6.8'',3,''BW'');'];
        uimenu(m_demoIDX,'Label',lab_str,'Callback',cba);

        lab_str = getWavMSG('Wavelet:dw2dRF:Lab_Example',2,'haar', ...
            tab,getWavMSG('Wavelet:moreMSGRF:EX2D_Name_Mandrill'));                               
        cba     = [beg_call_str ',''wmandril'',''haar'',2);'];
        uimenu(m_demoIDX,'Separator','on','Label',lab_str,'Callback',cba);

        lab_str = getWavMSG('Wavelet:dw2dRF:Lab_Example',2,'sym5', ...
            tab,getWavMSG('Wavelet:moreMSGRF:EX2D_Name_Gatlin'));                               
        cba     = [beg_call_str ',''wgatlin'',''sym5'',2,''BW'');'];
        uimenu(m_demoIDX,'Label',lab_str,'Callback',cba);

        lab_str = getWavMSG('Wavelet:dw2dRF:Lab_Example',3,'bior1.5', ...
            tab,getWavMSG('Wavelet:moreMSGRF:EX2D_Name_Belmont_1'));                               
        cba     = [beg_call_str ',''belmont1'',''bior1.5'',3);'];
        uimenu(m_demoIDX,'Label',lab_str,'Callback',cba);

        lab_str = getWavMSG('Wavelet:dw2dRF:Lab_Example',3,'coif1', ...
            tab,getWavMSG('Wavelet:moreMSGRF:EX2D_Name_Belmont_2'));                               
        cba     = [beg_call_str ',''belmont2'',''coif1'',3);'];
        uimenu(m_demoIDX,'Label',lab_str,'Callback',cba);

        lab_str = getWavMSG('Wavelet:dw2dRF:Lab_Example',4,'haar', ...
            tab,getWavMSG('Wavelet:moreMSGRF:EX2D_Name_jellyfish'));                               
        cba     = [beg_call_str ',''jellyfish256'',''haar'',4,''COL'');'];
        uimenu(m_demoCOL,'Label',lab_str,'Callback',cba);

        lab_str = getWavMSG('Wavelet:dw2dRF:Lab_Example',4,'haar', ...
            tab,getWavMSG('Wavelet:moreMSGRF:EX2D_Name_woodsculp'));                               
        cba     = [beg_call_str ',''woodsculp256.jpg'',''haar'',4,''COL'');'];
        uimenu(m_demoCOL,'Label',lab_str,'Callback',cba);
        
        lab_str = getWavMSG('Wavelet:dw2dRF:Lab_Example',4,'sym4', ...
            tab,getWavMSG('Wavelet:moreMSGRF:EX2D_Name_ArmsCOL'));                               
        cba     = [beg_call_str ',''arms.jpg'',''sym4'',4,''COL'');'];
        uimenu(m_demoCOL,'Label',lab_str,'Callback',cba);
        
        uimenu(m_imp_wrks,...
            'Label',getWavMSG('Wavelet:dw2dRF:Lab_ImportImg'), ...
            'Tag','Import_Img', ...
            'Callback',['dw2dmngr(''import_img'',' str_numwin ');'] ...
            );
        uimenu(m_imp_wrks,...
            'Label',getWavMSG('Wavelet:dw2dRF:Lab_ImportCfs'), ...
            'Tag','Import_Cfs', ...                        
            'Callback',['dw2dmngr(''import_cfs'',' str_numwin ');'] ...            
            );
        uimenu(m_imp_wrks,...
            'Label',getWavMSG('Wavelet:dw2dRF:Lab_ImportDec'),...
            'Tag','Import_Dec', ...            
            'Callback',['dw2dmngr(''import_dec'',' str_numwin ');'] ...
            );
        
        cb_beg = ['dw2dmngr(''exp_wrks'',' str_numwin];
        uimenu(m_exp_wrks,...
            'Label',getWavMSG('Wavelet:dw2dRF:Lab_ExportImg'), ...
            'Tag','Export_Img','Callback',[cb_beg ',''sig'');'] ...
            );
        uimenu(m_exp_wrks,...
            'Label',getWavMSG('Wavelet:dw2dRF:Lab_ExportCfs'), ...
            'Tag','Export_Cfs','Callback',[cb_beg ',''cfs'');'] ...
            );
        uimenu(m_exp_wrks,...
            'Label',getWavMSG('Wavelet:dw2dRF:Lab_ExportDec'),...
            'Tag','Export_Dec','Callback',[cb_beg ',''dec'');'] ...
            );        

        % Begin waiting.
        %---------------
        wwaiting('msg',win_dw2dtool,getWavMSG('Wavelet:commongui:WaitInit'));

        % General parameters initialization.
        %-----------------------------------
        dx = X_Spacing; dx2 = 2*dx;
        dy = Y_Spacing; dy2 = 2*dy;
        d_txt = (Def_Btn_Height-Def_Txt_Height);
        x_frame0   = pos_frame0(1);
        btn_width  = Def_Btn_Width;
        w_subframe = pos_frame0(3)-4*dx;
        w_util     = (pos_frame0(3)-9*dx)/3;
        pop_width  = Pop_Min_Width;

        % Position property of objects.
        %------------------------------
        xlocINI    = pos_frame0([1 3]);
        ybottomINI = pos_win(4)-3.5*Def_Btn_Height-dy2;
        y_low      = ybottomINI;

        bdx        = (pos_frame0(3)-1.5*btn_width)/2;
        x_left     = x_frame0+bdx;
        y_low      = y_low-Def_Btn_Height-3*dy2;
        h_btn      = 1.5*Def_Btn_Height;
        pos_anal   = [x_left, y_low, 1.5*btn_width, h_btn];

        x_left     = x_frame0+dx2;
        y_low      = pos_anal(2)-1.5*Def_Btn_Height-2*dy2;
        push_width = (pos_frame0(3)-3*dx2)/2;
        pos_stat   = [x_left, y_low, push_width, h_btn];

        pos_comp    = pos_stat;
        pos_comp(1) = pos_stat(1)+pos_stat(3)+dx2;

        y_low       = y_low-1.5*Def_Btn_Height-dy;
        pos_hist    = [x_left, y_low, push_width, h_btn];
        pos_deno    = pos_hist;
        pos_deno(1) = pos_hist(1)+pos_hist(3)+dx2;

        y_low          = y_low-Def_Btn_Height-2*dy2;
        w_txt          = 2*push_width+dx2-pop_width;
        pos_txt_declev = [x_left, y_low+d_txt/2, w_txt, Def_Txt_Height];

        pos_declev     = [x_left+4*dx2, y_low, pop_width, Def_Btn_Height];
        pos_declev(1)  = pos_txt_declev(1)+pos_txt_declev(3);

        y_low          = y_low-Def_Btn_Height-2*dy2;
        pos_pop_viewm  = [x_left, y_low, 2*push_width+dx2, Def_Btn_Height];

        y_low          = y_low-Def_Btn_Height-dy2;
        w_txt          = dx2+w_util;
        xl             = x_left+dx;
        yl             = y_low-Def_Txt_Height/2;
        pos_txt_full   = [xl, yl, w_txt, Def_Txt_Height];

        pos_pus_full      = zeros(4,4);
        xl                = pos_txt_full(1)+pos_txt_full(3)+dx;
        pos_pus_full(1,:) = [xl, y_low, w_util, Def_Btn_Height];

        pos_pus_full(2,:) = pos_pus_full(1,:);
        pos_pus_full(2,2) = pos_pus_full(2,2)-Def_Btn_Height;

        pos_pus_full(3,:) = pos_pus_full(1,:);
        pos_pus_full(3,1) = pos_pus_full(3,1)+pos_pus_full(3,3);

        pos_pus_full(4,:) = pos_pus_full(3,:);
        pos_pus_full(4,2) = pos_pus_full(4,2)-pos_pus_full(4,4);

        y_low       = y_low-2*Def_Btn_Height-3*dy2;
        pos_tdrag   = [x_left, y_low+10, w_subframe, Def_Btn_Height]; %High DPI y_low+10

        h_btn       = 1*Def_Btn_Height; %High DPI 1.25 to 1
        w_btn       = 1.25*push_width;
        x_btn       = x_frame0+(pos_frame0(3)-w_btn)/2;
        y_low       = y_low-h_btn;
        pos_pus_vis = [x_btn , y_low+10 , w_btn ,h_btn]; %High DPI y_low+10
        y_low       = y_low-h_btn;
        pos_pus_big = [x_btn , y_low+10 , w_btn ,h_btn];
        y_low       = y_low-h_btn;
        pos_pus_rec = [x_btn , y_low+10 , w_btn ,h_btn];

        % String property of objects.
        %----------------------------
        str_anal       = getWavMSG('Wavelet:commongui:Str_Anal');
        str_stat       = getWavMSG('Wavelet:commongui:Str_STAT');
        str_comp       = getWavMSG('Wavelet:commongui:Str_COMP');
        str_hist       = getWavMSG('Wavelet:commongui:Str_HIST');
        str_deno       = getWavMSG('Wavelet:commongui:Str_DENO');
        str_txtlev_dec = getWavMSG('Wavelet:dw2dRF:Str_Dec_At','');
        str_vallev_dec = int2str((1:max_lev_anal)');
        str_txt_drag   = getWavMSG('Wavelet:dw2dRF:Str_Txt_Drag');
        str_rec        = getWavMSG('Wavelet:dw2dRF:Str_Rec');
        str_vis        = getWavMSG('Wavelet:dw2dRF:Str_Vis');
        str_big        = getWavMSG('Wavelet:dw2dRF:Str_Big');
        str_pop_viewm  = { ...
            getWavMSG('Wavelet:dw2dRF:Str_VM_Sq'); ...
            getWavMSG('Wavelet:dw2dRF:Str_VM_Tr')
            };
        str_txt_full   = getWavMSG('Wavelet:dw2dRF:Str_Big');

        % Callback property of objects.
        %------------------------------
        cba_pus_anal   = ['dw2dmngr(''analyze'',' str_numwin ');'];
        cba_stat       = ['dw2dmngr(''stat'',' str_numwin ');'];
        cba_comp       = ['dw2dmngr(''comp'',' str_numwin ');'];
        cba_hist       = ['dw2dmngr(''hist'',' str_numwin ');'];
        cba_deno       = ['dw2dmngr(''deno'',' str_numwin ');'];
        cba_pop_declev = ['dw2dmngr(''view_dec'',' str_numwin ');'];
        cba_pop_viewm  = ['dw2dmngr(''view_mode'',' str_numwin ');'];

        % Command part of the window.
        %============================
        % Data, Wavelet and Level parameters.
        %------------------------------------
        utanapar('create',win_dw2dtool, ...
                 'xloc',xlocINI,'bottom',ybottomINI,...
                 'Enable','off',        ...
                 'wtype','dwt',         ...
                 'maxlev',max_lev_anal  ...
                 );
        comFigProp = {'Parent',win_dw2dtool,'Units',win_units};
        comPusProp = [comFigProp,'Style','pushbutton','Enable','off'];
        comPopProp = [comFigProp,'Style','Popupmenu','Enable','off'];
        comTxtProp = [comFigProp,'Style','Text', ...
                      'BackgroundColor',Def_FraBkColor];
        pus_anal = uicontrol(comPusProp{:},...
            'Position',pos_anal,...
            'String',str_anal,...
            'Tag',tag_pus_anal,...
            'Interruptible','On',...
            'Callback',cba_pus_anal...
            );

        uicontrol(comPusProp{:},...
            'Position',pos_stat,...
            'String',str_stat,...
            'Tag',tag_pus_stat,...
            'Callback',cba_stat...
            );

        uicontrol(comPusProp{:},...
            'Position',pos_comp,...
            'String',str_comp,...
            'Tag',tag_pus_comp,...
            'Callback',cba_comp...
            );

        uicontrol(comPusProp{:},...
            'Position',pos_hist,...
            'String',str_hist,...
            'Tag',tag_pus_hist,...
            'Callback',cba_hist...
            );

        uicontrol(comPusProp{:},...
            'Position',pos_deno,...
            'String',str_deno,...
            'Tag',tag_pus_deno,...
            'Callback',cba_deno...
            );

        uicontrol(comTxtProp{:},...
            'Position',pos_txt_declev,...
            'HorizontalAlignment','left',...
            'String',str_txtlev_dec...
            );

        uicontrol(comPopProp{:},...
            'Position',pos_declev,...
            'String',str_vallev_dec,...
            'Tag',tag_pop_declev,...
            'Callback',cba_pop_declev...
            );

        uicontrol(comPopProp{:},...
            'Position',pos_pop_viewm,...
            'String',str_pop_viewm,...
            'Callback',cba_pop_viewm,...
            'Tag',tag_pop_viewm...
            );

        txt_full = uicontrol(comTxtProp{:},...
            'HorizontalAlignment','Center',...
            'Position',pos_txt_full,...
            'String',str_txt_full,...
            'Tag',tag_txt_full...
            );

        tooltip = {...
            getWavMSG('Wavelet:dw2dRF:ViewOri'), ...
            getWavMSG('Wavelet:dw2dRF:ViewSyn'), ...
            getWavMSG('Wavelet:dw2dRF:ViewSel'), ...
            getWavMSG('Wavelet:dw2dRF:ViewDec') ...
            }; 
        pus_full = zeros(1,4);
        for k=1:4
            pus_full(k) = uicontrol(comPusProp{:},...
                'Position',pos_pus_full(k,:),...
                'String',sprintf('%.0f',k),...
                'UserData',0,...
                'TooltipString',deblank(tooltip{k}), ...
                'Tag',tag_pus_full(k,:)...
                );
        end


        txt_drag = uicontrol(comTxtProp{:},...
            'Position',pos_tdrag,...
            'HorizontalAlignment','left',...
            'String',str_txt_drag...
            );

        pus_vis = uicontrol(comPusProp{:},...
            'Position',pos_pus_vis,...
            'String',str_vis,...
            'Tag',tag_pus_vis...
            );

        pus_big  = uicontrol(comPusProp{:},...
            'Position',pos_pus_big,...
            'String',str_big,...
            'Tag',tag_pus_big...
            );

        pus_rec  = uicontrol(comPusProp{:},...
            'Position',pos_pus_rec,...
            'String',str_rec,...
            'Tag',tag_pus_rec...
            );


        % Adding colormap GUI.
        %---------------------
        utcolmap('create',win_dw2dtool, ...
                 'xloc',xlocINI,'bkcolor',Def_FraBkColor);

        %  Normalization.
        %----------------
        Pos_Graphic_Area = wfigmngr('normalize',win_dw2dtool, ...
            Pos_Graphic_Area,'On');

        % Callbacks update.
        %------------------
        utanapar('set_cba_num',win_dw2dtool,[m_files;pus_anal]);
        beg_cba     = ['dw2dmngr(''select'','  str_numwin ','];
        cba_pus_vis = [beg_cba , num2mstr(pus_vis) ');'];
        cba_pus_big = [beg_cba , num2mstr(pus_big) ');'];
        cba_pus_rec = [beg_cba , num2mstr(pus_rec) ');'];
        set(pus_vis,'Callback',cba_pus_vis);
        set(pus_big,'Callback',cba_pus_big);
        set(pus_rec,'Callback',cba_pus_rec);
        beg_cba = ['dw2dmngr(''fullsize'',' str_numwin ','];
        for k=1:4
            pus = pus_full(k);
            cba_pus_full = [beg_cba sprintf('%.0f',k) ');'];
            set(pus,'Callback',cba_pus_full);
        end

		% Add Context Sensitive Help (CSHelp).
		%-------------------------------------
		hdl_DW2D_FULLSIZE = [txt_full,pus_full(:)'];
		hdl_DW2D_SELECT = [txt_drag,pus_vis,pus_big,pus_rec];
		wfighelp('add_ContextMenu',win_dw2dtool,hdl_DW2D_FULLSIZE,'DW2D_FULLSIZE');		
		wfighelp('add_ContextMenu',win_dw2dtool,hdl_DW2D_SELECT,'DW2D_SELECT');
		%-------------------------------------

        % Memory for stored values.
        %--------------------------
        wmemtool('ini',win_dw2dtool,n_InfoInit,nb0_stored);
        wmemtool('ini',win_dw2dtool,n_param_anal,nb1_stored);
        wmemtool('ini',win_dw2dtool,n_miscella,nb3_stored);
        wmemtool('wmb',win_dw2dtool,n_param_anal,ind_act_option,option);
        wmemtool('wmb',win_dw2dtool,n_miscella,                 ...
                 ind_graph_area,Pos_Graphic_Area,ind_view_status,'none', ...
                 ind_save_status,'none',ind_sel_funct,select_funct       ...
                 );

        % Setting Colormap.
        %------------------
        cbcolmap('set',win_dw2dtool,'pal',{'pink',default_nbcolors});

        % Creating Axes: ImgIni, ImgBig, ImgVis, ImgSyn & ImgSel.
        %--------------------------------------------------------
        commonProp = {...
            'Parent',win_dw2dtool,     ...
            'Position',[0 0 1 1],      ...
            'Visible','off',           ...
            'XTickLabelMode','manual', ...
            'YTickLabelMode','manual', ...
            'XTicklabel',[],           ...
            'YTickLabel',[],           ...
            'Box','On',                ...
            'XGrid','off',             ...
            'YGrid','off'              ...
            };
        axes(commonProp{:},'Tag',tag_axeimgini);
        axes(commonProp{:},'Tag',tag_axeimgvis);
        axes(commonProp{:},'Tag',tag_axeimgbig);
        axes(commonProp{:},'Tag',tag_axeimgsyn);
        Axe_ImgSel = axes(commonProp{:},'Tag',tag_axeimgsel);

        % Creating AxeImgDec.
        %--------------------
        locProp = {commonProp{:}, ...
            'XTick',[],'YTick',[],'Tag',tag_axeimgdec};  %#ok<CCAT>
        Axe_ImgDec = zeros(1,4*max_lev_anal);
        for ind=1:4*max_lev_anal
            Axe_ImgDec(ind) = axes(locProp{:});
        end
        wmemtool('wmb',win_dw2dtool,tag_axeimgdec,1,Axe_ImgDec);

        % Creating Tree axes.
        %-------------------
        axe_figutil = axes(...
            'Parent',win_dw2dtool,    ...
            'Position',[0 0 1 1],     ...
            'XLim',[0 1],'YLim',[0 1],...
            'Visible','off',          ...
            'Tag',tag_axefigutil      ...
            );

        for k = 1:max_lev_anal+1
            line(...
                 'Parent',axe_figutil,        ...
                 'XData',[0 0],'YData',[0 0], ...
                 'LineWidth',2,'Color',linDW2D_Color, ...
                 'Visible','off',   ...
                 'UserData',k,      ...
                 'Tag',tag_linetree ...
                 );
        end
        fontsize = wmachdep('FontSize','normal');
        axeXColor = get(win_dw2dtool,'DefaultAxesXColor');
        for k = 1:max_lev_anal
            text(...
                 'Parent',axe_figutil,           ...
                 'String',['L_' sprintf('%.0f',k)],...
                 'FontSize',fontsize,            ...
                 'FontWeight','bold',            ...
                 'HorizontalAlignment','left',   ...
                 'Visible','off',                ...
                 'UserData',k,                   ...
                 'Color',axeXColor,              ...
                 'Tag',tag_txttree               ...
                 );
        end
        dw2darro('ini_arrow',win_dw2dtool);
        wboxtitl('create',axe_figutil,Axe_ImgSel,Col_BoxTitleSel,...
                          getWavMSG('Wavelet:dw2dRF:Img_SEL'),'off');

        % End waiting.
        %---------------
        wwaiting('off',win_dw2dtool);
		
    case 'close'
        called_win = wfindobj('figure','UserData',varargin{1});
        delete(called_win);

    otherwise
        errargt(mfilename,getWavMSG('Wavelet:moreMSGRF:Unknown_Opt'),'msg');
        error(message('Wavelet:FunctionArgVal:Invalid_Input'));
end
