function feature = extractsfeatures(Spatial, chanlocs)
warning off;

%% SPATIAL MAP
pattern = Spatial;
C = length(pattern);

%% CONVERT CARTESIAN TO SPHERICAL COORDINATES
[px,py] = cart2sph([chanlocs(1:C).X],[chanlocs(1:C).Y],[chanlocs(1:C).Z]);
theta = px/pi*180; % azimuth coordinate
phi = py/pi*180; % polar coordinate


%% Normalize spatial patterns to have unit variance
pattern = pattern./std(pattern);

%% Compute spatial range
spatialrange = log(max(pattern) - min(pattern));

%% Compute horizontal eye movement artifact
lhorizelec = find((theta >= 10) & (theta <= 60) & (phi <= 0) & (phi >= -40));
rhorizelec = find((theta <= -10) & (theta >= -60) & (phi <= 0) & (phi >= -40));
lrdiff = [];
for ii = 1:length(lhorizelec)
    for jj = 1:length(rhorizelec)
           lrdiff = [lrdiff abs(pattern(lhorizelec(ii)) - pattern(rhorizelec(jj)))];
    end
end

horizmov = log(max(lrdiff));

%% Compute saccade artifact
lhorizelec = (theta >= 40) & (theta <= 60) & (phi <= -10) & (phi >= -40);
rhorizelec = (theta <= -40) & (theta >= -60) & (phi <= -10) & (phi >= -40);
saccade = log(abs(mean(pattern(lhorizelec)) + mean(pattern(rhorizelec))));

%% Compute blink/vertical eye movement artifact
blinkelec = (abs(theta) <= 40) & (phi > -10) & (phi <= 0);
blink = log(abs(mean(pattern(blinkelec))));

%% Compute EKG artifact
Noutermost = 2*(sum((phi == min(phi)) & (abs(theta)>90))); % number of partitions on the border
inttheta = 360/(Noutermost); % interelectrode angle distance on the border
Theta = theta;
Theta(Theta < 0) = Theta(Theta < 0) + 360; % convert -180~0 to 180~360
outelec = [];
for n = 1:Noutermost
    a = (n-1)*inttheta;
    b = n*inttheta;
    I = find((Theta > (a-inttheta/2))& (Theta < (b-inttheta/2)));
    [Y, J] = min(phi(I));
    if (Y <= 0) & (Y >= -60)
       outelec = [outelec I(J)]; % outermost electrodes
    end
end
corrset = zeros(1, length(outelec));
elecset = [];
for m = 1:length(outelec) 
    elecset{m} = ((Theta(outelec)-Theta(outelec(m)) >= 0) & (Theta(outelec)-Theta(outelec(m)) <= 90)); 
%     corrset(m) = mean(ones(1, length(elecset)), pattern(outelec(elecset{m})));
    corrset(m) = mean(pattern(outelec(elecset{m})));
end
[Y, I] = max(corrset);
I2 = find((Theta(outelec)-Theta(outelec(I)) >= 90) & (Theta(outelec)-Theta(outelec(I)) <= 180));
corrtemplate = zeros(1, length(I2));
for l = 1:length(I2)
   template = zeros(1, length(pattern));
   template(elecset{I}) = 1;
   template(elecset{I2(l)}) = -1;
   corrtemplate(l) = corr(template', pattern);
end
if max(corrtemplate) > 0.6
    EKG = 1;
else
    EKG = 0;
end

%% Compute central mean activation
centralelec = (phi > 20);
centralact = log(abs(mean(pattern(centralelec))));

%% Compute frontal mean activation
frontalelec = (abs(theta) <= 60) & (abs(phi) <= 30);
frontalact = log(abs(mean(pattern(frontalelec))));

%% Compute occipital mean activation
occipelec = (abs(theta) <= 180) & (abs(theta) >= 155) & (phi <= 20);
occipact = log(abs(mean(pattern(occipelec))));

%% Compute left temporal mean activation
ltempelec = (theta >= 30) & (theta <= 150) & (phi <= 20);
ltempact = log(abs(mean(pattern([ltempelec]))));

%% Compute right temporal mean activation
rtempelec = (theta <= -30) & (theta >= -150) & (phi <= 20);
rtempact = log(abs(mean(pattern([rtempelec]))));

%% Compute the border activation feature
borderelec = find(phi == min(phi));
[Y, I] = max(abs(pattern));
if ismember(I(1), borderelec)
    borderact = 1;
else borderact = 0;
end

% %% Current source density
% channelfile = 'channel_bpEtkin.mat';
% SurfaceFile = 'tess_cortex_pial_low_3000V.mat';
% SkullFile = 'tess_innerskull_bem_1922V.mat';
% HeadFile = 'tess_head_bem_1922V.mat';
% Kernel = WMNE(channelfile, SurfaceFile, SkullFile, HeadFile);


% load ImagingKernel.mat;
% csdensity = log(sqrt(sum((Kernel*pattern).^2)));


% %current source density
% load inversematrix.mat;
% csdensity = log(sqrt(sum((M100*pattern).^2)));

%% APPEND FEATURES
feature = [spatialrange, blink, horizmov, frontalact, centralact, occipact, ltempact, rtempact, borderact]; % csdensity, frontal, saccade