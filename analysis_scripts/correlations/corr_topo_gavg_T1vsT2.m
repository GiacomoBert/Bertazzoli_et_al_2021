% calculate correlation via both  spearman correlation
% across pipelines in all conditions
clear
close all
clc

%% add toolboxes
%fieldtrip
addpath('Z:\Giacomo_20190828\Analysis\Plug-in_extensions\fieldtrip-20190905');
%Multiple testing toolbox
addpath('Z:\Giacomo_20190828\Analysis\Plug-in_extensions\MultipleTestingToolbox');
%eeglab
addpath('Z:\Giacomo_20190828\Analysis\Plug-in_extensions\eeglab2019_0');
%shaded error bar
addpath('Z:\Giacomo_20190828\Analysis\Plug-in_extensions\raacampbell-shadedErrorBar-67b0080');

%% define stat dir
TEPs_stat_dir='Z:\Giacomo_20190828\Analysis\Post_processing\TEPs\NEW_postprocessing\Correlation_dot_product_NEW\T1vsT2\topo';
pipeline={'ARTIST', 'TMSEEG', 'TESA', 'SOUND'};
pipeline_short={'AR', 'TM', 'TE', 'SO'};
pipeline_mod_dir1={'', '_newICsel', '_newICsel', '\SOUND_normal_filt_check'};
pipeline_mod_dir2={'', '_newICsel', '_newICsel', '_normal_filt_check'};

area_full={'parietal','prefrontal'};
area={'lP', 'lPF'};
session={'T1', 'T2'};
%% load grand avg TEPs

data_gavg=cell(4,2,2);
data_gavg_ind=cell(4,2,2);
for ses_count=1:length(session)
    for area_count=1:length(area)
        
        for pip_count=1:length(pipeline)
            data_gavg(pip_count, area_count, ses_count)= {load(['Z:\Giacomo_20190828\Analysis\Post_processing\NEW_postprocessing\' pipeline{pip_count} '_postprocessing' pipeline_mod_dir1{pip_count} '\' pipeline{pip_count} pipeline_mod_dir2{pip_count} '_left_' area_full{area_count} '_' session{ses_count} '\' pipeline{pip_count} pipeline_mod_dir2{pip_count} '_gr_avg_' area{area_count} '_' session{ses_count} '.mat'])};
            data_gavg_ind(pip_count, area_count, ses_count)= {load(['Z:\Giacomo_20190828\Analysis\Post_processing\NEW_postprocessing\' pipeline{pip_count} '_postprocessing' pipeline_mod_dir1{pip_count} '\' pipeline{pip_count} pipeline_mod_dir2{pip_count} '_left_' area_full{area_count} '_' session{ses_count} '\' pipeline{pip_count} pipeline_mod_dir2{pip_count} '_gr_avg_' area{area_count} '_' session{ses_count} '_KeepSubj.mat'])};
            
        end
    end
end


%% spearman corr
time=struct2array(data_gavg_ind{1,1,1});
time=time.time;


% for ses_count=1:length(session)
for area_count=1:length(area)
    figure('Renderer', 'painters', 'Position', [10 10 1200 600]); %open a bigger window
    %         comp_count=1;
    for pip_count=1:length(pipeline)
        %             cond_count=pip_count+1;
        %             while cond_count<=4
        
        for sub_count=1:16
            
            for time_count=1:600
                selected_data1= struct2array(data_gavg_ind{pip_count,area_count,1}); %session1
                selected_data1= selected_data1.individual(sub_count,:,time_count);
                
                selected_data2= struct2array(data_gavg_ind{pip_count,area_count,2}); %session2
                selected_data2= selected_data2.individual(sub_count,:,time_count);
                
                x1=selected_data1';
                x2=selected_data2';
                
                RHO(sub_count,time_count) = ...
                    corr(x1,x2,'Type','Spearman', 'tail', 'both');
                
                
                %                     scorr = @(x1,x2)(corr(x1,x2,'type','Spearman'));
                %                     stat_set=statset('UseParallel',1);
                %                     [CI_bootci(sub_count).subject(:,time_count), boorstat]=bootci(10000, {scorr, x1, x2}, 'type', 'norm', 'Options', stat_set);
                
            end
        end
        z_RHO=atanh(RHO); %fisher z-transform
        
        
        %test the null hypothesis that the correlation are
        %different RHO = 0
        [H_spear(pip_count,:),PvalSPEAR(pip_count,:)]= ...
            ttest(z_RHO, 0, 'alpha', 0.025, 'dim', 1, 'tail', 'both');
        
        % correct for multiple comparisons
        [FDR_PvalSPEAR(pip_count,:), FDR_alphaSPEAR(pip_count,:), FDR_H_spear(pip_count,:)]=...
            fdr_BH(PvalSPEAR(pip_count,:), 0.025);
        
        mean_z_RHO=mean(z_RHO,1); %mean z RHO
        std_z_RHO=std(z_RHO,1); %std z RHO
        
        for inverse_count=1:length(mean_z_RHO) %inverse fisher transform
            mean_RHO(inverse_count)=(exp(2*mean_z_RHO(inverse_count))-1)/ (exp(2*mean_z_RHO(inverse_count))+1);
            std_RHO(inverse_count)=(exp(2*std_z_RHO(inverse_count))-1)/ (exp(2*std_z_RHO(inverse_count))+1);
        end
        
        %inverse fisher transform for the signle subject RHO value.
        
        for sub_count=1:16
            for inverse_count=1:length(z_RHO)
                all_RHO(sub_count, inverse_count)=(exp(2*z_RHO(sub_count, inverse_count))-1)/ (exp(2*z_RHO(sub_count, inverse_count))+1);
                all_std_RHO(sub_count, inverse_count)=(exp(2*z_RHO(sub_count, inverse_count))-1)/ (exp(2*z_RHO(sub_count, inverse_count))+1);
            end
        end
        
        % %                 % compute the range around the mean, because SEM or
        % %                 %CI exeed 1
        % %                 SEM=std_RHO / (sqrt(16)*sqrt(6/(6-1)));
        % %                 percentile95(1,:)=prctile(all_RHO, 95);
        
        %%
        s=plot(time, mean_RHO, 'LineWidth',2, 'LineStyle', '-');
        hold on
        
        h1=yline(0.6, 'LineStyle', '--', 'LineWidth',1.5);
        set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        
        h2=yline(0.8, 'LineStyle', '--', 'LineWidth',1.5);
        set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        
        %draw a gray rectangle in the interpolated area for the TMS pulse
        timeRect = [-0.001 0.006];
        ymax = 80;
        ymin = -40;
        rectangle('Position',[timeRect(1) ymin (timeRect(2)-timeRect(1)) ymax],'FaceColor',[0.7 0.7 0.7 0.3], 'LineStyle','none');
        
        ylim([-0.2 1]);
        xlim([-0.05, 0.35]);
        xticklabels({ '-50', '0', '50' ,'100', '150', '200', '250', '300', '350'})
        yticks([-0.2 0 0.2 0.4 0.6 0.8 1]);
        
        xlabel('Time (ms)');
        ylabel('Spearman \rho');
        % title(title_data);
        set(gca,'TickDir','out'); % The only other option is 'in'
        
        set(gca,'box','off') %remove ticks on the figure
        
        ax = gca;
        c = ax.Color;
        ax.FontWeight  = 'bold';
        ax.FontSize    = 20;
        
        
        % disp conditions
        data_gavg_ind{pip_count,area_count,1}
        data_gavg_ind{pip_count,area_count,2}
        %
        %                 cond_count=cond_count+1;
        %                 comp_count=comp_count+1;
        %             end
    end
    % legend({'AR_TM' 'AR_TE', 'AR_SO', 'TM_TE', 'TM_SO', 'TE_SO'})
%    error('remove this line to run the script'); %just to prevent that I run the script accidentaly and I overwrite results while I do some checks
    hold off
    
    savefig([TEPs_stat_dir '\corrspearPIP_correct_NEW_def'...
        area{area_count}]);
    saveas(gcf,...
        [TEPs_stat_dir '\corrspearPIP_correct_NEW_def'...
        area{area_count}],...
        'tif');
    % save
    save([TEPs_stat_dir '\stat_corrspearPIP_correct_NEW_def'...
        area{area_count}],...
        'FDR_PvalSPEAR', 'FDR_alphaSPEAR','FDR_H_spear')
    
    close all
    
end
% end
%%


