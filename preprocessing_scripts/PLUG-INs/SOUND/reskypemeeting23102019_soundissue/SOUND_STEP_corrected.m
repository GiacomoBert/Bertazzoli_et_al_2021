% AN EXAMPLE MATLAB SCRIPT FOR USING DDWIENER AND SOUND
%
% This is an EEGLAB-compatible MATLAB script that demonstrates how DDWiener
% and SOUND can be used to clean a partially noisy EEG dataset.
%
% Required toolbox to run this script: BioSig, TESA, bva-io (plugin to import BrainVision data)
% 18 September 2017 : Tuomas Mutanen and Marta Bortoletto.
% Department of Neuroscience and Biomedical Engineering (NBE), School of Science, Aalto University
% .........................................................................


% Relevant Variables that are created in the script
% EEG_orig %original continous data
% EEG_janela %continous after janelacut
% EEG_artrej %epoched after rejection (148 epochs)
% EEG_clean %after sound
% EEG_sICA %after ICA (performed after sound)
% EEG_clean_muscle %after SSP-SIR

% .........................................................................


%% Define directories and subjects
% STEP 1: Data pre-processing up to artifact rejection

clear all;		% prepare a clean workspace with no variables
close all;		% close all open events

save_everything = 1;

% Add eeglab;
addpath('\\home.org.aalto.fi\tpmutane\data\Documents\TMS-EEG preprcosessing\eeglab2019_0');
eeglab

% Add SOUND:
addpath('\\home.org.aalto.fi\tpmutane\data\Documents\TMS-EEG preprcosessing\sound');

% define directories
% all directories must have been created
locpath=('\\home.org.aalto.fi\tpmutane\data\Documents\TMS-EEG preprcosessing\eeglab2019_0\plugins\dipfit\standard_BESA\standard-10-5-cap385.elp'); %R

datalocation =('\\home.org.aalto.fi\tpmutane\data\Documents\TMS-EEG preprcosessing\TMS-EEG_preprocessing\GRGS');

datafolder = '\\home.org.aalto.fi\tpmutane\data\Documents\TMS-EEG preprcosessing\SSPSIR_SOUND\SOUND_problem';

stepfolder='\\home.org.aalto.fi\tpmutane\data\Documents\TMS-EEG preprcosessing\SSPSIR_SOUND\Steps';

figurefolder = '\\home.org.aalto.fi\tpmutane\data\Documents\TMS-EEG preprcosessing\SSPSIR_SOUND\Figures';

outputfolder=' \\home.org.aalto.fi\tpmutane\data\Documents\TMS-EEG preprcosessing\SSPSIR_SOUND\Final';

SSP_SIRfolder='\\home.org.aalto.fi\tpmutane\data\Documents\TMS-EEG preprcosessing\SSPSIR_SOUND\SSPSIR';

subject = 'MLEN_lP_T1';


%Event selection for last steps from 6
trigger= [16];

    %% 
    disp('Shell 2 started')
    
    EEG = pop_loadset('filename',[subject '_2PreProc.set'],'filepath', datafolder);
    
    %For testing sound and speeding up things let's study only a subsection
    EEG = pop_select( EEG, 'time',[-0.005 0.031] );
   %  DownSample data (5000 Hz to 1024 Hz)
    %   EEG = pop_resample( EEG, 1024);
    figure; pop_timtopo(EEG, [-5  30], [10  15  20  25 ])

    LFM_sphere= load([stepfolder '\','GRGS_rP_LFM_sphere.txt']);

    
      %% 4. Use SOUND to clean the channel-specific noise
    tic
    disp('Shell 4 started')
    
    % Set the SOUND parameters:
    % Set the reqularization parameter of the minimum-norm estimate when
    % finding the noise-free current estimates. This can be adjusted if clear
    % under- or overcorrection is observed. lambda_value = 0.1 was used in the original
    % article.
    lambda_value = 0.1;
    
    % Set the number of iterations when estimating the channel-specific noise
    % estimates. Iter = 5 was found to be sufficient in the datasets studied in
    % the original article.
    iter = 5;
    
    EEG_clean = EEG;
    
   chanN = size(EEG.data,1);
   
   % For loop for estimating the noise at each time sample separately:
   
   step = 1;
   
   all_sigmas = zeros(size(EEG.data,1),size(EEG.data,2));
   
   for i = 1:(size(EEG.data,2))% 1:ms_samples:(size(EEG.data,2) - ms_samples + 1 )

        % Find the number of trials:
        N_tr= size(EEG.data,3);
                
        % Run the SOUND algorithm:
        
        % Inform about the progress of SOUND:
        EEG.times(i)
        
        % Estimate the noise for each sample
  
           tmp_data = reshape(EEG.data(1:61, i ,:), [61, N_tr] );

            if i ==1
                sigmas = ones(size(tmp_data,1),1);
                [y_solved, sigmas] = DDWiener(tmp_data);
              
            end
            
            [~, ~, sigmas] = SOUND(tmp_data, LFM_sphere, iter, lambda_value,sigmas,1);
            
            all_sigmas(:,step) = sigmas;
            step = step + 1;%         
            
%        disp(['processing ' all_datasets{maincounter} ' _ dataset N° ' num2str(maincounter)])
       
   end
   
   % Now we clean each sample by smoothly estimating the sample-specific
   % noises with a smooth time-window
   
   scaling_w = 2;
   smoothness = 1;
   all_sigmas_ext = [repmat(all_sigmas(:,1),[1  scaling_w]), all_sigmas,  repmat(all_sigmas(:,end),[1  scaling_w])];
   step =  scaling_w + 1;
   
   smoothing_func = gaussmf(1:(2*scaling_w + 1),[smoothness  scaling_w + 1]);
   smoothing_func = smoothing_func/sum(smoothing_func); figure; plot(smoothing_func)
     smoothing_func = repmat(  smoothing_func, [size(EEG.data,1), 1]);

  
   for i = 1:size(EEG.data,2)% 1:ms_samples:(size(EEG.data,2) - ms_samples + 1 )

        % Find the number of trials:
        N_tr= size(EEG.data,3);
                
        % Run the SOUND algorithm:
        
        % Inform about the progress of SOUND:
        EEG.times(i)
        
            
            % Make the final estimate for the noise distribution for this
            % sample, by using the noise estimate and its neigbourhood:
            
            sigmas =  sum(all_sigmas_ext(:,(step-scaling_w):(step +scaling_w) ).*smoothing_func,2);
                     
            % Final cleaning step:
            
            W = diag(1./sigmas);
            WL = W*LFM_sphere;
            WLLW = WL*WL';
            
            data = reshape(EEG.data(1:61, i , :), [chanN, N_tr] );
            
            x = WL'*((WLLW + lambda_value*trace(WLLW)/chanN*eye(chanN))\(W*data));
            corrected_data = LFM_sphere*x;
            EEG_clean.data(1:61, i, :) = reshape(corrected_data,[size(corrected_data,1), 1, N_tr]);
        
            step = step + 1;
              
   end


    %----------------------------------------------------------------
    %save data after SOUND:
     EEG.setname= 'After SOUND';
    EEG = pop_saveset( EEG_clean, 'filename',[subject '_3Sound_TM.set'],'filepath',stepfolder);
    %---------------------------------------------------------------------
    
    disp('Shell 4 ended')
    
     figure; pop_timtopo(EEG, [-5  30], [10  15  20  25 ])
     toc
     

    
   