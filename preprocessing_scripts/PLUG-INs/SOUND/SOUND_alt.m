function [corrected_data, x, sigmas,dn] = SOUND_alt(data, LFM, iter,lambda0,sigmas, estimate_just_noise)
% This function performs the SOUND algorithm for a given data.
%
% .........................................................................
% 24 September 2017: Tuomas Mutanen, NBE, Aalto university  
% .........................................................................

chanN = size(data,1);
if nargin < 5
sigmas = ones(chanN,1);
[y_solved, sigmas] = DDWiener(data);
end
if nargin<4
    lambda0 = 1;
end

if nargin < 6
    estimate_just_noise = 0;
end

% Number of time points
T = size(data,2);
LL = LFM*LFM';

chanPerms = zeros(size(data,1),size(data,1)-1);
for i = 1:size(data,1)
    chanPerms(i,:) = setdiff(1:chanN,i);
end

% Going through all the channels as many times as requested
for k=1:iter
    sigmas_old = sigmas;
    
    %Evaluating each channel in a random order
    for i=randperm(chanN)
        chan = chanPerms(i,:);
            % Defining the whitening operator with the latest noise
            % estimates
            W = diag(1./sigmas);
            
            % Computing the whitened version of the lead field
            %WL = (W(chan,chan))*(LFM(chan,:));
            WL =  (1./sigmas(chan)).* LFM(chan,:);
           % WLLW = W(chan,chan)*( LL(chan,chan)*(W(chan,chan))' );
             WLLW =  (1./sigmas(chan)).*( LL(chan,chan)*(W(chan,chan))' );
             
            % Computing the MNE, the Wiener estimate in the
            % studied channel, as well as the corresponding noise estimate
            %x = (WL)'*((WLLW + lambda0*trace(WLLW)/(chanN-1)*eye(chanN-1))\((W(chan,chan))*(data(chan,:))));
            %y_solved = LFM*x;

        %    y_solved = LFM*( (WL)'*((WLLW + lambda0*trace(WLLW)/(chanN-1)*eye(chanN-1))\((W(chan,chan))*(data(chan,:)))));
             y_solved = LFM*( (WL)'*((WLLW + lambda0*trace(WLLW)/(chanN-1)*eye(chanN-1))\((1./sigmas(chan)).*(data(chan,:)))));
            sigmas(i) = sqrt((y_solved(i,:)-data(i,:))*(y_solved(i,:)-data(i,:))')/sqrt(T);
    end
    
    % Following and storing the convergence of the algorithm
    dn(k) = max(abs(sigmas_old - sigmas)./sigmas_old);

end

if estimate_just_noise
    x = [];
    corrected_data = [];

else
% Final data correction based on the final noise-covariance estimate.

            W = diag(1./sigmas);
            WL = W*LFM;
            WLLW = WL*WL';
            x = WL'*((WLLW + lambda0*trace(WLLW)/chanN*eye(chanN))\(W*data));
            corrected_data = LFM*x;
            
            
end
end