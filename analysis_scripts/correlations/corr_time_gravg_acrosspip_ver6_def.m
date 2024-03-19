% calculate correlation via spearman correlation
% across pipelines in all conditions
% corraltion in time

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
%colormaps
addpath('Z:\Giacomo_20190828\Analysis\Plug-in_extensions\Colormaps-20200430T170805Z-001\Colormaps');

%% define stat dir
TEPs_stat_dir='Z:\Giacomo_20190828\Analysis\Post_processing\TEPs\NEW_postprocessing\Correlation_dot_product_NEW\corr_TIME';
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

% initialize eeglab
eeglab
close
%% spearman corr
time=struct2array(data_gavg_ind{1,1,1}).time;
chan_num=length(data_gavg{1,1,1}.ARTIST_gr_avg_lP_T1.label);

interval=([0.150 0.350]);
index_interval=[find(round(time*1000)==round(interval(1)*1000)) ...
    find(round(time*1000)==round(interval(2)*1000))];

% dummy fieldtrip structure
dummy_ft=struct2array(data_gavg{1,1,1});

%% load dummy chanloc
% % EEG = pop_loadset(...
% %     'filename', 'ANAD_lP_T1_TMSEEG_Processed.set',...
% %     'filepath', 'Z:\Giacomo_20190828\Analysis\TMSEEG_outputs_newICsel');
% %
% % % generate channel Tp9 (original online ref)
% % EEG.data(63,length(EEG.data)) = 0; %adds one 0 microV data channel (ref TP9)
% % EEG.nbchan = 63; %change nbchan to keep consistency in number of data channels
% % EEG.chanlocs(63).type = 'EEG';
% % EEG.chanlocs(63).labels = 'Tp9';
% % EEG= pop_chanedit(EEG,  'lookup', 'standard-10-5-cap385.elp'); %adds coordinates
% % EEG = pop_interp(EEG, 63, 'spherical');

%%

for ses_count=1:length(session)
    for area_count=1:length(area)
        comp_count=1;
        for pip_count=1:length(pipeline)
            cond_count=pip_count+1;
            while cond_count<=4
                
                parfor sub_count=1:16
                    
                    for chan_count=1:chan_num
                        selected_data1 =squeeze(struct2array(data_gavg_ind{pip_count,area_count,ses_count}).individual(sub_count,chan_count,:));
                        selected_data2 =squeeze(struct2array(data_gavg_ind{cond_count,area_count,ses_count}).individual(sub_count,chan_count,:));
                        
                        x1=selected_data1(index_interval(1):index_interval(2)); %do not include baseline and interp period
                        x2=selected_data2(index_interval(1):index_interval(2));
                        
                        RHO(sub_count,chan_count) = ...
                            corr(x1,x2,'Type','Spearman', 'tail', 'both');
                    end
                end
                z_RHO=atanh(RHO); %fisher z-transform
                if find(z_RHO==Inf) %in there are Inf value due to 1 correlation, make it a NaN
                    z_RHO(find(z_RHO==Inf))=NaN;
                end
                
                % test significance against 0
                % test significance
                [H_spear(comp_count,:),PvalSPEAR(comp_count,:)]= ttest(z_RHO, 0, 'alpha', 0.025, 'dim', 1, 'tail', 'both');
                % correct for multiple comparisons
                [FDR_PvalSPEAR(comp_count,:), FDR_alphaSPEAR(comp_count,:), FDR_H_spear(comp_count,:)]=...
                    fdr_BH(PvalSPEAR(comp_count,:), 0.025);
                
                mean_z_RHO=mean(z_RHO,1, 'omitnan'); %mean z RHO
                std_z_RHO=std(z_RHO,1); %std z RHO
                
                for inverse_count=1:length(mean_z_RHO) %inverse fisher transform
                    mean_RHO(inverse_count)=(exp(2*mean_z_RHO(inverse_count))-1)/ (exp(2*mean_z_RHO(inverse_count))+1);
                    std_RHO(inverse_count)=(exp(2*std_z_RHO(inverse_count))-1)/ (exp(2*std_z_RHO(inverse_count))+1);
                end
                
                %%
                mean_RHO_round=round(mean_RHO', 1); %trasform value to be plotted in the topoplot
                mean_RHO_chart=num2str(mean_RHO_round);
                for x=1:chan_num
                    EEG.chanlocs(x).labels = mean_RHO_chart(x,:);
                    %   dummy_ft.label{x}=  num2str(mean_RHO_chart(x,:));
                end
                
                %% time ft
                figure('Renderer', 'painters', 'Position', [10 10 1200 600]); %open a bigger window
                
                % prepare dummy structure
                dummy_ft.avg=mean_RHO';
                dummy_ft.time=0.1;
                
                mask=FDR_H_spear(comp_count,:);
                
                cfg=[];
                %                 cfg.parameter='stat';
                cfg.layout = 'EEG1010_mod_GIA.lay';
                
                cfg.zlim= [0 1];
                
                cfg.highlight          = 'on';
                cfg.highlightchannel   =  find(mask)';
                cfg.highlightsymbol    = '*';
                cfg.gridscale  =300;
                cfg.colormap='jet';
                cfg.highlightcolor     = [255, 250, 250]/255;
                %cfg.highlightsize      = highlight marker size (default = 6)
                cfg.highlightfontsize  = 20;
                cfg.comment='no';
                
                t= ft_topoplotER(cfg, dummy_ft);
                
                current_obj=gcf;
                H=findobj(gcf);
                if sum(mask)>0 && sum(mask) ~=	63 %if and only if there are significant channels
                    H(17).LineWidth=5; %significant cahnnels bold
                    H(18).LineWidth=1; % topo lines bold
                    H(11).LineWidth=5; % cahnnels bo
                    
                elseif sum(mask)==0
                    H(17).LineWidth=1; % topo lines bold
                    H(11).LineWidth=5; % cahnnels bo
                    
                elseif sum(mask)==63
                    H(16).LineWidth=5; %significant cahnnels bold
                    H(17).LineWidth=1; % topo lines bold
                end
                
                %%
                % %                 topoplot(mean_RHO, EEG.chanlocs,...
                % %                     'style', 'blank',...
                % %                     'plotrad' ,0.8,...
                % %                     'maplimits', [0 1],...
                % %                     'headrad'  , 0.7,...
                % %                     'conv'    , 'on',...
                % %                     'electrodes' ,'labels',...
                % %                     'emarker', {[],'k',10,1});
                %%
                hold on
                %                 c=colorbar;
                %                 c.Label.String = 'Dot Product';
                %                 c.Label.FontSize = 12;
                
                set(gca,'box','off') %remove ticks on the figure
                
                % legend({'AR_TM' 'AR_TE', 'AR_SO', 'TM_TE', 'TM_SO', 'TE_SO'})
                
                hold off
                
                savefig([TEPs_stat_dir '\corrspearPIP_TIME_correct_NEW_'...
                    pipeline_short{pip_count} '_'...
                    pipeline_short{cond_count} '_'...
                    area{area_count} '_'...
                    session{ses_count} '_'...
                    num2str(interval(1)*1000) '_'...
                    num2str(interval(2)*1000)]);
                saveas(gcf,...
                    [TEPs_stat_dir '\corrspearPIP_TIME_correct_NEW_'...
                    pipeline_short{pip_count} '_'...
                    pipeline_short{cond_count} '_'...
                    area{area_count} '_'...
                    session{ses_count} '_'...
                    num2str(interval(1)*1000) '_'...
                    num2str(interval(2)*1000)], ...
                    'tif');
                % save
                save([TEPs_stat_dir '\stat_corrspearPIP_TIME_correct_NEW_'...
                    pipeline_short{pip_count} '_'...
                    pipeline_short{cond_count} '_'...
                    area{area_count} '_'...
                    session{ses_count} '_' ...
                    num2str(interval(1)*1000) '_'...
                    num2str(interval(2)*1000)],...
                    'FDR_PvalSPEAR', 'FDR_alphaSPEAR','FDR_H_spear')
                
                close all
                
                
                cond_count=cond_count+1;
                comp_count=comp_count+1;
            end
        end
        
    end
end
%%


