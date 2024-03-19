% eegplugin_bidsvalidator() - EEGLAB plugin for validating BIDS dataset
%
% Usage:
%   >> eegplugin_bidsvalidator(fig, trystrs, catchstrs);
%
% Inputs:
%   fig        - [integer]  EEGLAB figure
%   trystrs    - [struct] "try" strings for menu callbacks.
%   catchstrs  - [struct] "catch" strings for menu callbacks.

function vers = eegplugin_bidsvalidator(fig, trystrs, catchstrs)

    vers = '1.1';
    if nargin < 3
        error('eegplugin_bidsvalidator requires 3 arguments');
    end
    
%     % add folder to path
%     % ------------------
%     p = which('pop_importbids.m');
%     p = p(1:findstr(p,'pop_importbids.m')-1);
%     if ~exist('pop_importbids')
%         addpath( p );
%     end
    
    % find data menu
    % ---------------------
    bids = findobj(fig, 'label', 'BIDS tools');
    if isempty(bids)
        menui3 = findobj(fig, 'label', 'File');
        bids = uimenu( menui3, 'label', 'BIDS tools', 'separator', 'on', 'position', 5, 'userdata', 'startup:on;study:on');
    end
    
    % create BIDS menus
    % -----------------
    if isempty(get(bids, 'children'))
        comvalidatebids = [ trystrs.no_check 'pop_validatebids();' catchstrs.add_to_hist ];
        uimenu( bids, 'label', 'Validate BIDS dataset', 'callback', comvalidatebids, 'userdata', 'startup:on;study:on');
    end
