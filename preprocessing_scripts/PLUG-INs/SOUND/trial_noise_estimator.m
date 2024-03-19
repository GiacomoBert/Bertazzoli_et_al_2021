function [ERP,est_noise,tris] = trial_noise_estimator(data)
% This function uses DDWiener to estimate the noise level in each trial.
%
% .........................................................................
% 15 January 2016 : Tuomas Mutanen, NBE, Aalto university  
% .........................................................................

disp('Using DDWiener to compute the ERP with least noise. This may take several minutes.')

[CN,T,TriN] = size(data);

tris = zeros(CN,T,TriN);
est_noise = zeros(CN,TriN);
for i=1:CN
    disp(['Evaluating trials of channel ',num2str(i),'/',num2str(CN)]);
    D = squeeze(data(i,:,:))';
    C = D*D';
    gammaV = mean(diag(C));

    
    for j = 1:TriN
        trials = setdiff(1:TriN,j);
        tri_est = C(j,trials)*((C(trials,trials) + gammaV*eye(TriN-1))\D(trials,:));
        tris(i,:,j) = tri_est;
        est_noise(i,j) = sqrt((tri_est-D(j,:))*(tri_est-D(j,:))')/sqrt(T);

    end
end
        scals = permute(repmat(1./est_noise./repmat(sum(1./est_noise,2),[1,TriN]),[1,1,T]),[1,3,2]);

        ERP = sum(scals.*data,3);


end

