
% plot topographies at selected latencies

clear
close all
clc

%% add toolboxes
%fieldtrip
addpath('Z:\Giacomo_20190828\Analysis\Plug-in_extensions\fieldtrip-20190905');
%eeglab
addpath('Z:\Giacomo_20190828\Analysis\Plug-in_extensions\eeglab2019_0');
%colormaps
addpath('Z:\Giacomo_20190828\Analysis\Plug-in_extensions\Colormaps-20200430T170805Z-001\Colormaps');

%% define stat dir
TEPs_dir='Z:\Giacomo_20190828\Analysis\Post_processing\TEPs\NEW_postprocessing\Correlation_dot_product_NEW\T1vsT2\illustrative_topographies';
pipeline={'ARTIST', 'TMSEEG', 'TESA', 'SOUND'};
pipeline_short={'AR', 'TM', 'TE', 'SO'};
pipeline_mod_dir1={'', '_newICsel', '_newICsel', '\SOUND_normal_filt_check'};
pipeline_mod_dir2={'', '_newICsel', '_newICsel', '_normal_filt_check'};

area_full={'parietal','prefrontal'};
area={'lP', 'lPF'};
session={'T1', 'T2'};
subnum=16;
%% load grand avg TEPs
data_gavg=cell(4,2,2);
data_avg=cell(4,2,2);

for ses_count=1:length(session)
    for area_count=1:length(area)
        for pip_count=1:length(pipeline)
            data_gavg(pip_count, area_count, ses_count)= {load(['Z:\Giacomo_20190828\Analysis\Post_processing\NEW_postprocessing\' pipeline{pip_count} '_postprocessing' pipeline_mod_dir1{pip_count} '\' pipeline{pip_count} pipeline_mod_dir2{pip_count} '_left_' area_full{area_count} '_' session{ses_count} '\' pipeline{pip_count} pipeline_mod_dir2{pip_count} '_gr_avg_' area{area_count} '_' session{ses_count} '.mat'])};
            data_avg(pip_count, area_count, ses_count)= {load(['Z:\Giacomo_20190828\Analysis\Post_processing\NEW_postprocessing\' pipeline{pip_count} '_postprocessing' pipeline_mod_dir1{pip_count} '\' pipeline{pip_count} pipeline_mod_dir2{pip_count} '_left_' area_full{area_count} '_' session{ses_count} '\' pipeline{pip_count} pipeline_mod_dir2{pip_count} '_all_avg_' area{area_count} '_' session{ses_count} '.mat'])};
        end
    end
end

%% define latency

latency=[0.180 0.180];


%% Plot grand grand average for all conditions for each area, averaged
for ses_count=1:2
    for area_count=2
        for pip_count=1:4
            
            selected_data1=struct2array(data_gavg{pip_count, area_count, ses_count}); %struct2array to remove usless intermediate
            
            %%
            data_gavg{pip_count, area_count, ses_count} %just to show the pip
            
            %% make topoplots with stats
            f= figure('Renderer', 'painters', 'Position', [10 10 1200 600]); %open a bigger window
            
            cfg=[];
            %                 cfg.parameter='stat';
            cfg.layout = 'EEG1010_mod_GIA.lay';
            cfg.xlim= latency;

            %lP
            %    cfg.zlim= [-3 4]; % 10 30 ms
            %    cfg.zlim= [-2 4]; % 40 60 ms
            %    cfg.zlim= [-5 5]; % 100 130 ms
            %    cfg.zlim= [-8 12]; % 180 210 ms
            %    cfg.zlim= [-3 4]; % 280 310 ms
            
            %lPF8
            %  cfg.zlim= [-2 4]; % 10 30 ms
            %  cfg.zlim= [-2 4]; % 40 60 ms
            %  cfg.zlim= [-6 6]; % 100 130 ms
            %  cfg.zlim= [-10 12]; % 180 210 ms
            % cfg.zlim= [-3 4]; % 280 310 ms
            
%             if ses_count==1 && area_count==1
%                 
%                 cfg.zlim=[-9 12];
%             elseif ses_count==1 && area_count==2
%                 cfg.zlim=[-11 13];
%                 
%             elseif ses_count==2 && area_count==1
%                 cfg.zlim=[-9 11];
%                 
%                 
%             elseif ses_count==1 && area_count==2
%                 cfg.zlim=[-10 13];
%                 
%             end
                
    cfg.zlim=[-11 11];
%    colorbar

            cfg.gridscale=300;
            cfg.colormap='inferno';
            cfg.comment='no';
            
            t= ft_topoplotER(cfg, selected_data1);
            hold on
            %             title([pipeline{pip_count} ' ',...
            %                 area{area_count} ' '...
            %                 session{ses_count}]);
            current_obj=gcf;
            H=findobj(gcf);
            H(12).LineWidth=1; %topo line bold
            
%             hold off
            
            %% save
            
            savefig([TEPs_dir '\' pipeline{pip_count} '_'...
                area{area_count} '_'...
                session{ses_count} '_'...
                num2str(latency(1)*1000) '_'...
                num2str(latency(2)*1000) '_TOPO' ]);
            saveas(gca,...
                [TEPs_dir '\' pipeline{pip_count} '_'...
                area{area_count} '_'...
                session{ses_count} '_'...
                num2str(latency(1)*1000) '_'...
                num2str(latency(2)*1000) '_TOPO' ],...
                'tif');
            
            
            
            %                            keyboard
            close all
            
            
        end
    end
end

