% Author: Matthew Frehlich, Ye Mei, Luis Garcia Dominguez, Faranak Farzan
%         2016
%         Ben Schwartzmann 
%         2017

% eegdatapro_show() - Displays EEG data in a butterfly plot after the processing 
% step given by afterstep

% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

function [] = eegdatapro_show(afterstep, laststep)

%Check if previous step was done
if tmseeg_previous_step(afterstep+1) 
    return 
end

global backcolor VARS 

hfig = figure('Units','normalized',...
              'name',['EEG After Step ' num2str(afterstep)],...
              'numbertitle','off',...
              'resize','on',...
              'color',backcolor,...
              'Position',[0 0 0.7 0.7],...
              'DockControls','off');
axes('Position',[0.2 0.1 0.75 0.9]);   

data_button = uicontrol('style','pushbutton',...
     'Units','normalized',...
     'Position',[0.02 0.90 0.1 0.05],...
     'Parent',hfig,...
     'string','View Data',...
     'Callback',{@view_data,afterstep,laststep}); %#ok
spectrum_button = uicontrol('style','pushbutton',...
     'Units','normalized',...
     'Position',[0.02 0.85 0.1 0.05],...
     'Parent',hfig,...
     'string','View Spectrum',...
     'Callback',{@view_spectrum,afterstep}); %#ok
zoom_button = uicontrol('style','pushbutton',...
     'Units','normalized',...
     'Position',[0.02 0.75 0.1 0.05],...
     'Parent',hfig,...
     'string','Zoom on Data',...
     'Callback',{@view_zoom,afterstep,laststep}); %#ok
          
pre_pulse_deletion = 0;
 
xshowmin = VARS.XSHOWMIN;
xshowmax = VARS.XSHOWMAX;
yshowlimit = VARS.YSHOWLIMIT;


%Adjust display limits based on epoching
if VARS.XSHOWMIN < VARS.EPCH_STRT
    xshowmin = VARS.EPCH_STRT;
end

if VARS.XSHOWMAX > VARS.EPCH_END
    xshowmax = VARS.EPCH_END;
end


%Load proper dataset
[~, EEG] = eegdatapro_load_step(afterstep + 1);

if afterstep == laststep
    pre_pulse_deletion = 1;
end


%Create topo plot with EEGLAB timtopo() command  
if(pre_pulse_deletion)
    data_temp = squeeze(nanmean(EEG.data,3));
else
    data_temp = squeeze(nanmean(EEG.data,3));

    if(isfield(EEG,{'TMS_period2remove_b'}))
        ix = min(EEG.TMS_period2remove_b);
        rm_pulse_fill = NaN(size(data_temp,1),length(EEG.TMS_period2remove_b));
        data_temp = cat(2,data_temp(:,1:ix-1),rm_pulse_fill,data_temp(:,ix:end));
    end

    if(isfield(EEG,{'TMS_period2remove_1'}))
        %Insert NaN values to fill space where TMS pulse was removed
        ix       = min(EEG.TMS_period2remove_1);
        rm_pulse_fill = NaN(size(data_temp,1),length(EEG.TMS_period2remove_1));
        data_temp = cat(2,data_temp(:,1:ix-1),rm_pulse_fill,data_temp(:,ix:end));
    end

end

EEGtimes =  min(EEG.times):1000/EEG.srate:max(EEG.times);
data = data_temp(:,EEGtimes>=xshowmin & EEGtimes<=xshowmax);

timtopo(data, EEG.chanlocs,'limits',[xshowmin xshowmax -yshowlimit yshowlimit]);

end

function view_spectrum(varargin)

h = findall(gcf,'Type','axes');
delete(h);

global VARS 
afterstep=varargin{3};

[~, EEG] = eegdatapro_load_step(afterstep + 1);          

axes('Position',[0.2 0.15 0.70 0.75]);          
pop_spectopo( EEG,1,[EEG.xmin*1000 EEG.xmax*1000],'EEG','percent',100,'freq',[2 6 20],'freqrange',VARS.SPECTRUM_RNG,'electrodes','on');

end


function view_data(varargin)

h = findall(gcf,'Type','axes');
delete(h);

global VARS 

afterstep=varargin{3};
laststep=varargin{4};

pre_pulse_deletion = 0;
 
xshowmin = VARS.XSHOWMIN;
xshowmax = VARS.XSHOWMAX;
yshowlimit = VARS.YSHOWLIMIT;


%Adjust display limits based on epoching
if VARS.XSHOWMIN < VARS.EPCH_STRT
    xshowmin = VARS.EPCH_STRT;
end

if VARS.XSHOWMAX > VARS.EPCH_END
    xshowmax = VARS.EPCH_END;
end


%Load proper dataset
[~, EEG] = eegdatapro_load_step(afterstep + 1);

if afterstep == laststep
    pre_pulse_deletion = 1;
end


%Create topo plot with EEGLAB timtopo() command  
if(pre_pulse_deletion)
    data_temp = squeeze(nanmean(EEG.data,3));
else
    data_temp = squeeze(nanmean(EEG.data,3));

    if(isfield(EEG,{'TMS_period2remove_b'}))
        ix = min(EEG.TMS_period2remove_b);
        rm_pulse_fill = NaN(size(data_temp,1),length(EEG.TMS_period2remove_b));
        data_temp = cat(2,data_temp(:,1:ix-1),rm_pulse_fill,data_temp(:,ix:end));
    end

    if(isfield(EEG,{'TMS_period2remove_1'}))
        %Insert NaN values to fill space where TMS pulse was removed
        ix       = min(EEG.TMS_period2remove_1);
        rm_pulse_fill = NaN(size(data_temp,1),length(EEG.TMS_period2remove_1));
        data_temp = cat(2,data_temp(:,1:ix-1),rm_pulse_fill,data_temp(:,ix:end));
    end
    
end

EEGtimes =  min(EEG.times):1000/EEG.srate:max(EEG.times);
data = data_temp(:,EEGtimes>=xshowmin & EEGtimes<=xshowmax);

axes('Position',[0.2 0.1 0.75 0.9]);   
timtopo(data, EEG.chanlocs,'limits',[xshowmin xshowmax -yshowlimit yshowlimit]);       

end


function view_zoom(varargin)

%Load EEG Data, display in custom plot allowing zoom feature
global VARS

afterstep=varargin{3};
laststep=varargin{4};

figure('units','normalized',...
        'menubar','none',...
        'numbertitle','off',...
        'toolbar','none',...
        'name',['After Step ' num2str(afterstep)],...
        'position',[0 0 .9 .9 ]);

pre_pulse_deletion = 0;
 
xshowmin = VARS.XSHOWMIN;
xshowmax = VARS.XSHOWMAX;

%Adjust display limits based on epoching
if VARS.XSHOWMIN < VARS.EPCH_STRT
    xshowmin = VARS.EPCH_STRT;
end

if VARS.XSHOWMAX > VARS.EPCH_END
    xshowmax = VARS.EPCH_END;
end


%Load proper dataset
[~, EEG] = eegdatapro_load_step(afterstep + 1);

if afterstep == laststep
    pre_pulse_deletion = 1;
end


%Create topo plot with EEGLAB timtopo() command  
if(pre_pulse_deletion)
    data_temp = squeeze(nanmean(EEG.data,3));
else
    data_temp = squeeze(nanmean(EEG.data,3));

    if(isfield(EEG,{'TMS_period2remove_b'}))
        ix = min(EEG.TMS_period2remove_b);
        rm_pulse_fill = NaN(size(data_temp,1),length(EEG.TMS_period2remove_b));
        data_temp = cat(2,data_temp(:,1:ix-1),rm_pulse_fill,data_temp(:,ix:end));
    end

    if(isfield(EEG,{'TMS_period2remove_1'}))
        %Insert NaN values to fill space where TMS pulse was removed
        ix       = min(EEG.TMS_period2remove_1);
        rm_pulse_fill = NaN(size(data_temp,1),length(EEG.TMS_period2remove_1));
        data_temp = cat(2,data_temp(:,1:ix-1),rm_pulse_fill,data_temp(:,ix:end));
    end

end

EEGtimes =  min(EEG.times):1000/EEG.srate:max(EEG.times);
data = data_temp(:,EEGtimes>=xshowmin & EEGtimes<=xshowmax);

x = EEGtimes(EEGtimes>=xshowmin & EEGtimes<=xshowmax);
y = squeeze(nanmean(data,3));
plot(x,y)

titlestr  = ['Data after processing step ' num2str(afterstep) ...
    ', use cursor to zoom in, shift + click to zoom out'];
title(titlestr)
xlabel('Time(ms)');
ylabel(['Amplitude (' char(0181) 'V)']);
zoom on

end