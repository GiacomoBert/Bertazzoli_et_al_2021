% one way anova across pipeline. Performa a rmANOVA with cluster-based
% multiple comparisons correction as implemented in fieldtrip. Do it
% separately for all conditions lP_T1 lP_T2 lPF_T1 lPF_T2

clear
close all
clc

%% define stat dir
TEPs_stat_dir='Z:\Giacomo_20190828\Analysis\Post_processing\TEPs\NEW_postprocessing\TEPs_stat_output_NEW\Cluster_based_analysis_NEW\6_350';

%% add toolboxes
%fieldtrip
addpath('Z:\Giacomo_20190828\Analysis\Plug-in_extensions\fieldtrip-20190905');

%% load grand avg TEPs
% lP_T1
load('Z:\Giacomo_20190828\Analysis\Post_processing\NEW_postprocessing\ARTIST_postprocessing\ARTIST_left_parietal_T1\ARTIST_gr_avg_lP_T1.mat');
load('Z:\Giacomo_20190828\Analysis\Post_processing\NEW_postprocessing\TMSEEG_postprocessing_newICsel\TMSEEG_newICsel_left_parietal_T1\TMSEEG_newICsel_gr_avg_lP_T1.mat');
load('Z:\Giacomo_20190828\Analysis\Post_processing\NEW_postprocessing\TESA_postprocessing_newICsel\TESA_newICsel_left_parietal_T1\TESA_newICsel_gr_avg_lP_T1.mat');
load('Z:\Giacomo_20190828\Analysis\Post_processing\NEW_postprocessing\SOUND_postprocessing\SOUND_normal_filt_check\SOUND_normal_filt_check_left_parietal_T1\SOUND_normal_filt_check_gr_avg_lP_T1.mat');
% lP_T2
load('Z:\Giacomo_20190828\Analysis\Post_processing\NEW_postprocessing\ARTIST_postprocessing\ARTIST_left_parietal_T2\ARTIST_gr_avg_lP_T2.mat');
load('Z:\Giacomo_20190828\Analysis\Post_processing\NEW_postprocessing\TMSEEG_postprocessing_newICsel\TMSEEG_newICsel_left_parietal_T2\TMSEEG_newICsel_gr_avg_lP_T2.mat');
load('Z:\Giacomo_20190828\Analysis\Post_processing\NEW_postprocessing\TESA_postprocessing_newICsel\TESA_newICsel_left_parietal_T2\TESA_newICsel_gr_avg_lP_T2.mat');
load('Z:\Giacomo_20190828\Analysis\Post_processing\NEW_postprocessing\SOUND_postprocessing\SOUND_normal_filt_check\SOUND_normal_filt_check_left_parietal_T2\SOUND_normal_filt_check_gr_avg_lP_T2.mat');
% lPF_T1
load('Z:\Giacomo_20190828\Analysis\Post_processing\NEW_postprocessing\ARTIST_postprocessing\ARTIST_left_prefrontal_T1\ARTIST_gr_avg_lPF_T1.mat');
load('Z:\Giacomo_20190828\Analysis\Post_processing\NEW_postprocessing\TMSEEG_postprocessing_newICsel\TMSEEG_newICsel_left_prefrontal_T1\TMSEEG_newICsel_gr_avg_lPF_T1.mat');
load('Z:\Giacomo_20190828\Analysis\Post_processing\NEW_postprocessing\TESA_postprocessing_newICsel\TESA_newICsel_left_prefrontal_T1\TESA_newICsel_gr_avg_lPF_T1.mat');
load('Z:\Giacomo_20190828\Analysis\Post_processing\NEW_postprocessing\SOUND_postprocessing\SOUND_normal_filt_check\SOUND_normal_filt_check_left_prefrontal_T1\SOUND_normal_filt_check_gr_avg_lPF_T1.mat');
% lPF_T2
load('Z:\Giacomo_20190828\Analysis\Post_processing\NEW_postprocessing\ARTIST_postprocessing\ARTIST_left_prefrontal_T2\ARTIST_gr_avg_lPF_T2.mat');
load('Z:\Giacomo_20190828\Analysis\Post_processing\NEW_postprocessing\TMSEEG_postprocessing_newICsel\TMSEEG_newICsel_left_prefrontal_T2\TMSEEG_newICsel_gr_avg_lPF_T2.mat');
load('Z:\Giacomo_20190828\Analysis\Post_processing\NEW_postprocessing\TESA_postprocessing_newICsel\TESA_newICsel_left_prefrontal_T2\TESA_newICsel_gr_avg_lPF_T2.mat');
load('Z:\Giacomo_20190828\Analysis\Post_processing\NEW_postprocessing\SOUND_postprocessing\SOUND_normal_filt_check\SOUND_normal_filt_check_left_prefrontal_T2\SOUND_normal_filt_check_gr_avg_lPF_T2.mat');


%% load avg all subj
% lP_T1
load('Z:\Giacomo_20190828\Analysis\Post_processing\NEW_postprocessing\ARTIST_postprocessing\ARTIST_left_parietal_T1\ARTIST_all_avg_lP_T1.mat');
load('Z:\Giacomo_20190828\Analysis\Post_processing\NEW_postprocessing\TMSEEG_postprocessing_newICsel\TMSEEG_newICsel_left_parietal_T1\TMSEEG_newICsel_all_avg_lP_T1.mat');
load('Z:\Giacomo_20190828\Analysis\Post_processing\NEW_postprocessing\TESA_postprocessing_newICsel\TESA_newICsel_left_parietal_T1\TESA_newICsel_all_avg_lP_T1.mat');
load('Z:\Giacomo_20190828\Analysis\Post_processing\NEW_postprocessing\SOUND_postprocessing\SOUND_normal_filt_check\SOUND_normal_filt_check_left_parietal_T1\SOUND_normal_filt_check_all_avg_lP_T1.mat');
% lP_T2
load('Z:\Giacomo_20190828\Analysis\Post_processing\NEW_postprocessing\ARTIST_postprocessing\ARTIST_left_parietal_T2\ARTIST_all_avg_lP_T2.mat');
load('Z:\Giacomo_20190828\Analysis\Post_processing\NEW_postprocessing\TMSEEG_postprocessing_newICsel\TMSEEG_newICsel_left_parietal_T2\TMSEEG_newICsel_all_avg_lP_T2.mat');
load('Z:\Giacomo_20190828\Analysis\Post_processing\NEW_postprocessing\TESA_postprocessing_newICsel\TESA_newICsel_left_parietal_T2\TESA_newICsel_all_avg_lP_T2.mat');
load('Z:\Giacomo_20190828\Analysis\Post_processing\NEW_postprocessing\SOUND_postprocessing\SOUND_normal_filt_check\SOUND_normal_filt_check_left_parietal_T2\SOUND_normal_filt_check_all_avg_lP_T2.mat');
% lPF_T1
load('Z:\Giacomo_20190828\Analysis\Post_processing\NEW_postprocessing\ARTIST_postprocessing\ARTIST_left_prefrontal_T1\ARTIST_all_avg_lPF_T1.mat');
load('Z:\Giacomo_20190828\Analysis\Post_processing\NEW_postprocessing\TMSEEG_postprocessing_newICsel\TMSEEG_newICsel_left_prefrontal_T1\TMSEEG_newICsel_all_avg_lPF_T1.mat');
load('Z:\Giacomo_20190828\Analysis\Post_processing\NEW_postprocessing\TESA_postprocessing_newICsel\TESA_newICsel_left_prefrontal_T1\TESA_newICsel_all_avg_lPF_T1.mat');
load('Z:\Giacomo_20190828\Analysis\Post_processing\NEW_postprocessing\SOUND_postprocessing\SOUND_normal_filt_check\SOUND_normal_filt_check_left_prefrontal_T1\SOUND_normal_filt_check_all_avg_lPF_T1.mat');
% lPF_T2
load('Z:\Giacomo_20190828\Analysis\Post_processing\NEW_postprocessing\ARTIST_postprocessing\ARTIST_left_prefrontal_T2\ARTIST_all_avg_lPF_T2.mat');
load('Z:\Giacomo_20190828\Analysis\Post_processing\NEW_postprocessing\TMSEEG_postprocessing_newICsel\TMSEEG_newICsel_left_prefrontal_T2\TMSEEG_newICsel_all_avg_lPF_T2.mat');
load('Z:\Giacomo_20190828\Analysis\Post_processing\NEW_postprocessing\TESA_postprocessing_newICsel\TESA_newICsel_left_prefrontal_T2\TESA_newICsel_all_avg_lPF_T2.mat');
load('Z:\Giacomo_20190828\Analysis\Post_processing\NEW_postprocessing\SOUND_postprocessing\SOUND_normal_filt_check\SOUND_normal_filt_check_left_prefrontal_T2\SOUND_normal_filt_check_all_avg_lPF_T2.mat');


%% Prepare neighbours

cfg = [];
cfg.layout = 'EEG1010_mod_GIA.lay';
cfg.method = 'triangulation';
%cfg.elec = elec;
cfg.feedback='yes';

neighbours = ft_prepare_neighbours(cfg, ARTIST_gr_avg_lP_T1);
close
%% ANOVA

% parameters
cfg=[];
cfg.channel     = 'all';
cfg.latency     = [0.006 0.350]; %no baseline
cfg.avgoverchan = 'no';
cfg.avgovertime = 'no';
cfg.parameter   = 'avg';

cfg.method       = 'montecarlo'; %get Monte-Carlo estimates of the significance probabilities and/or critical values from the permutation distribution,
cfg.statistic = 'ft_statfun_depsamplesFunivariate'; %rm anova 
cfg.correctm = 'cluster';
cfg.alpha=0.05; 
cfg.clusteralpha = 0.05;
cfg.tail  = 1; % For the Fstatistic only cfg.tail = 1 makes sense.
cfg.numrandomization = 1000;% number of draws from the permutation distribution

cfg.clusterstatistic = 'maxsum';
cfg.clustertail = 1;
cfg.minnbchan = 2;               % minimum number of neighborhood channels that is
% required for a selected sample to be included
% in the clustering algorithm (default=0).
cfg.neighbours    = neighbours;

nsubj = 16;
design = zeros(2,4*nsubj); %construct design matrix
design(1,1:nsubj)=1:nsubj; %num subjects repeated for 4 condition in the first row of design
design(1,nsubj+1:nsubj*2)=1:nsubj;
design(1,nsubj*2+1:nsubj*3)=1:nsubj;
design(1,nsubj*3+1:nsubj*4)=1:nsubj;

design(2,1:nsubj)=1; % in the second row the number of condition per each subj, so 1 till nsubj 2 till nsubj*2 etc...
design(2,nsubj+1:nsubj*2)=2;
design(2,nsubj*2+1:nsubj*3)=3;
design(2,nsubj*3+1:nsubj*4)=4;

%The design matrix in a within-UO design is different from the design matrix
%in a between-UO design. In the design matix for a within-UO design, you have
%to specify the unit variable. The unit variable specifies the units that have
%produced the different condition-specific data structures. For example, consider a
% hypothetical study with 4 subjects and 2 experimental conditions. The design matrix may
% then look like this: design = [1 2 3 4 1 2 3 4; 1 1 1 1 2 2 2 2 ]. The first row of this
% matrix is the unit variable: it specifies that the first subject produced the first and
% the fifth data structure, the second subject produced the second and the sixth data structure,
% etc. The second row of the design matrix is the independent variable.

cfg.design = design;

cfg.uvar  = 1;
cfg.ivar  = 2;

%   If you use a cluster-based statistic, you can specify the following
%   options that determine how the single-sample or single-voxel
%   statistics will be thresholded and combined into one statistical
%   value per cluster.

%lP T1
[stat_TEPs_rANOVA_lP_T1_NEW] = ft_timelockstatistics(cfg,... %left parietal
    ARTIST_all_avg_lP_T1{:},...
    TESA_newICsel_all_avg_lP_T1{:},...
    TMSEEG_newICsel_all_avg_lP_T1{:},...
    SOUND_normal_filt_check_all_avg_lP_T1{:});

save([TEPs_stat_dir '\stat_TEPs_rANOVA_lP_T1_NEW'], 'stat_TEPs_rANOVA_lP_T1_NEW');


%lP T2
[stat_TEPs_rANOVA_lP_T2_NEW] = ft_timelockstatistics(cfg,... %left parietal
    ARTIST_all_avg_lP_T2{:},...
    TESA_newICsel_all_avg_lP_T2{:},...
    TMSEEG_newICsel_all_avg_lP_T2{:},...
    SOUND_normal_filt_check_all_avg_lP_T2{:});

save([TEPs_stat_dir '\stat_TEPs_rANOVA_lP_T2_NEW'], 'stat_TEPs_rANOVA_lP_T2_NEW');

%lPF T1
[stat_TEPs_rANOVA_lPF_T1_NEW] = ft_timelockstatistics(cfg,... %left parietal
    ARTIST_all_avg_lPF_T1{:},...
    TESA_newICsel_all_avg_lPF_T1{:},...
    TMSEEG_newICsel_all_avg_lPF_T1{:},...
    SOUND_normal_filt_check_all_avg_lPF_T1{:});

save([TEPs_stat_dir '\stat_TEPs_rANOVA_lPF_T1_NEW'], 'stat_TEPs_rANOVA_lPF_T1_NEW');

%lPF T2
[stat_TEPs_rANOVA_lPF_T2_NEW] = ft_timelockstatistics(cfg,... %left parietal
    ARTIST_all_avg_lPF_T2{:},...
    TESA_newICsel_all_avg_lPF_T2{:},...
    TMSEEG_newICsel_all_avg_lPF_T2{:},...
    SOUND_normal_filt_check_all_avg_lPF_T2{:});

save([TEPs_stat_dir '\stat_TEPs_rANOVA_lPF_T2_NEW'], 'stat_TEPs_rANOVA_lPF_T2_NEW');

