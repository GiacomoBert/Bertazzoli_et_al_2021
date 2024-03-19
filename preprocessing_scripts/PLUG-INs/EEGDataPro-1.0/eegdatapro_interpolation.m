% Author: Ye Mei, Luis Garcia Dominguez,Faranak Farzan   
%         2015
%         Ben Schwartzmann
%         2017

% eegdatapro_interpolation() - interpolated removed channels based off original
% channel configuration from data loading.  Also rereferences data and pads
% for time segment removed in tmseeg_rm_TMS_artifact.
%
% Inputs: S        - parent GUI information (structure)
%         step_num - Step number for current cleaning step in workflow

% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

function eegdatapro_interpolation(S,step_num)

%Check if previous steps were done
if tmseeg_previous_step(step_num)
    return 
end

global basepath basefile VARS %#ok
h = msgbox('Interpolating now!');
[files, EEG] = eegdatapro_load_step(step_num);

%Compare current and original channel configurations
badchan = ~ismember({EEG.chanloc_orig.labels},{EEG.chanlocs.labels});

% Interpolate channels if they have been removed
if isfield(EEG,'nbchan_o')
    if EEG.nbchan < EEG.nbchan_o
        EEG = pop_interp(EEG, EEG.chanloc_orig, 'spherical');
    end
end
EEG.channels_interpolated = find(badchan);
EEG = eeg_checkset(EEG);

re_ref = questdlg('Re-reference to average?');
if strcmp(re_ref,'Yes')
    EEG = pop_reref(EEG, []);
    EEG = eeg_checkset( EEG );
end

%Add back removed time segment
r=S.routine;

if sum(strcmp(r.option_name,'REMOVE TMS ARTIFACT'))~=0
    [~,ind]=max((strcmp(r.option_name,'REMOVE TMS ARTIFACT'))==1);
    add_buffer = questdlg(strcat(['Add buffer time to for time periods deleted in step ' num2str(ind) '?']));
    if strcmp(add_buffer,'Yes')
        EpochSecs = EEG.epoch_length;
        EEG = tmseeg_addTMSTimeBack(EEG, EpochSecs); 
    end
end

save(fullfile(basepath,[basefile '_eegdatapro_settings.mat']), 'VARS');
tmseeg_step_check(files, EEG, S, step_num)

close(h)
