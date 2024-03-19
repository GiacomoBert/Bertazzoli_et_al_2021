%% Adding User-Defined Wavelets
% This example shows how to add a new wavelet family in 
% Wavelet Toolbox(TM) using the |wavemngr| function. Wavelet Toolbox 
% contains five types of wavelets (WTs):
%
% * WT = 1  orthogonal wavelets
% * WT = 2  biorthogonal wavelets
% * WT = 3  wavelets with scale function
% * WT = 4  wavelets without scale function
% * WT = 5  complex wavelets without scale function
%
% Adding your own wavelet requires two steps: 
%
% * *_Define the parameters._*
% For type 1 or type 2 wavelets, we must define five parameters  
% and pass them to the |wavemngr| function with the *add* option. 
% For types 3, 4, or 5, we must define a sixth parameter.
% 
% * *_Define a file._*
% The wavelets of a family must be built using either a
% MAT-file or a MATLAB(R) file. This file is called by the functions which 
% operate on wavelet families.
%
% In this example, we show how to use the |wavemngr| function to add  
% wavelet families of types 1, 2, 3, and 4. Examples of the file 
% associated with a wavelet family are also shown.

% Copyright 2008-2012 The MathWorks, Inc.

%% Create a Wavelet Family with One Type 1 Wavelet
% First, we create an orthogonal wavelet family containing only 
% one wavelet. 

%%
% Define the five parameters associated with the wavelet family. See the
% |wavemngr| function for more details.
Current_DIR = cd;   % Save the current directory name.
cd(tempdir);        % Work in a temporary directory.
familyName      = 'MyWAVE T1';
familyShortName = 'mywa';
familyWaveType  = 1;
familyNums      = '';
fileWaveName    = 'mywa.mat';

%%
% Create a filter associated with the wavelet and save it in the MAT-file
% *mywa.mat*. This step must be done before using the new wavelet family.
%%
% This filter will be used by the |wfilters| function to build four filters
% used by the Discrete Wavelet Transform (DWT).
sq3 = sqrt(3);
mywa = [(1+sq3) (3+sq3) (3-sq3) (1-sq3)]/8;
save mywa mywa

%%
% Add the new wavelet family to the stack of wavelet families.
wavemngr('add',familyName,familyShortName,familyWaveType, ...
    familyNums,fileWaveName)

%%
% Verify that the wavelet family was created.
wavemngr('read')

%%
% All the available functions for orthogonal wavelets can operate on the  
% new wavelet. Some examples are shown below.
%
% First we can compute the four filters associated with this wavelet 
% using the |wfilters| function.
[LO_D,HI_D,LO_R,HI_R] = wfilters('mywa')

%%
% We can also display the scaling and wavelet functions using the
% |wavefun| function.
wname = 'mywa';
wavefun(wname,'plot',7);
a = findobj(gcf,'Type','axes');
axis(a,'tight')

%%
% When using the Wavelet and the Wavelet Packet Display tools 
% from the |wavemenu| GUI tool, the new wavelet is also available.

%%
% Here is a screen shot of the Wavelet Display Tool.
%%
% <<imgmywa.png>>

%%
% All the discrete analysis functions, including |dwt|, 
% |idwt|, |wavedec|, etc. can operate on the new wavelet.
load freqbrk; % Load a signal
x = freqbrk;
[C,L] = wavedec(x,2,wname);
A2 = wrcoef('a',C,L,wname,2); % Approximation of level 2
D1 = wrcoef('d',C,L,wname,1); % Detail of level 1
D2 = wrcoef('d',C,L,wname,2); % Detail of level 2
subplot(4,1,1); plot(x,'r'); axis tight
title('Signal, Approximation A2 and Details D2 and D1')
ylabel('S','Rotation',0)
subplot(4,1,2); plot(A2,'Color',[0.5 0.5 0.9]); axis tight
ylabel('A2','Rotation',0)
subplot(4,1,3); plot(D2,'Color',[0.5 0.9 0.5]); axis tight
ylabel('D2','Rotation',0)
subplot(4,1,4); plot(D1,'Color',[0.5 0.9 0.5]); axis tight
ylabel('D1','Rotation',0)
xlabel('Time or Space')

%%
% The signals from top to bottom in the figure above are 
% *x* (the signal), *A2* (Approximation of level 2), *D2* 
% (Detail of level 2) and *D1* (Detail of level 1).

%%
% If the wavelet is an orthogonal wavelet, the reconstruction must
% be perfect, which is *x = A2 + D2 + D1*.
% This is the case as we can see below. 
maxdiff = max(abs(x-(A2+D1+D2)))

%%
% All the continuous analysis functions, including |cwt|, |wscalogram|, etc. and 
% the corresponding GUI tools, can also operate on the new wavelet.
scales = 1:1:128;
clf; cwt(x,scales,wname,'3Dplot');

%%
% *NOTE*: There is no guarantee that the four filters generated by the 
% |wfilters| function are associated with an orthogonal wavelet.
% To be an orthogonal wavelet, several mathematical properties must
% be verified. For more details, see Chapter 6 of the User's Guide.
% So it is almost impossible to find a good filter at random.
% Yet, even if the created wavelet is not orthogonal, all the discrete 
% and continuous analysis functions and the corresponding GUI tools, 
% can also operate on the new wavelet. But reconstruction will not be 
% perfect and therefore *x* will not be equal to *A2 + D2 + D1*.

%%
% Delete the new wavelet. 
wavemngr('del',familyShortName);

%%
% Verify that the wavelet family is deleted.
wavemngr('read')

%% Create a Wavelet Family with Several Type 1 Wavelets
% Now, we create an orthogonal wavelet family containing multiple wavelets.
familyName      = 'MyWAVE F1';
familyShortName = 'lemw';
familyWaveType  = 1;
familyNums      = '1 2 3 **';
fileWaveName    = 'lemwavf';

%%
% The variable |familyNums|, which defines the wavelet order, 
% contains three predefined numbers (1, 2 and 3)  
% and the special sequence **** which lets you work use any number.  
% To be available for use, we must define the selected number in the 
% file associated with the wavelet family.
%
%%
% Add a new wavelet family to the stack of wavelet families.
wavemngr('add',familyName,familyShortName,familyWaveType, ...
    familyNums,fileWaveName)

%%
% Verify that the wavelet family is created.
wavemngr('read')

%%
% View the content of the new family.
wavemngr('read',1)

%%
% As we can see above, the *MyWAVE F1* family, whose short name is *lemw*,
% contains three predefined wavelets, and lets you use a "postdefined" 
% wavelet, because we defined the family numbers using ****.

%%
% We now display the scaling and wavelet functions of this new wavelet family.
wname = 'lemw2';
wavefun(wname,'plot',7);
a = findobj(gcf,'Type','axes');
axis(a,'tight')

%%
% We can also display the scaling and wavelet functions for a wavelet
% defined by the special sequence ****. This wavelet has a number greater 
% than three and was defined in the file associated with the wavelet family.
wname = 'lemw5';
wavefun(wname,'plot',7);
a = findobj(gcf,'Type','axes');
axis(a,'tight')

%%
% This new wavelet *lemw5* can also be used to analyze signals.
load noisdopp; % Load a signal
x = noisdopp;
[C,L] = wavedec(x,3,wname);
A = wrcoef('a',C,L,wname,3);  % Approximation of level 3
D = wrmcoef('d',C,L,wname);   % Details of level 1, 2 and 3.
subplot(5,1,1); plot(x,'r'); axis tight;
title('Signal, Approximation A3 and Details D3, D2 and D1')
ylabel('S','Rotation',0)
subplot(5,1,2); plot(A,'Color',[0.5 0.5 0.9]); axis tight
ylabel('A3','Rotation',0)
subplot(5,1,3); plot(D(3,:),'Color',[0.5 0.9 0.5]); axis tight
ylabel('D3','Rotation',0)
subplot(5,1,4); plot(D(2,:),'Color',[0.5 0.9 0.5]); axis tight
ylabel('D2','Rotation',0)
subplot(5,1,5); plot(D(1,:),'Color',[0.5 0.9 0.5]); axis tight
ylabel('D1','Rotation',0)
xlabel('Time or Space')

%%
% Delete the new wavelet. 
wavemngr('del',familyShortName);

%% Create a Wavelet Family of Type 2
% First, we create a biorthogonal wavelet family containing only 
% one wavelet. 
familyName      = 'MyWAVE T2';
familyShortName = 'mywb';
familyWaveType  = 2;
familyNums      = '';
fileWaveName    = 'mywb.mat';

%%
% Create the two filters associated with the biorthogonal wavelet and save 
% them in a MAT-file, *mywb.mat*.
Rf = [1/2 1/2];
Df = [7/8  9/8  1/8  -1/8]/2;
save mywb Rf Df

%% 
% These filters give a true biorthogonal wavelet, because the
% filters are well designed. For example, if you choose:
% Df = [-1/16 1/16 1/2 1/2 1/16 -1/16], you recover the 'bior1.3' wavelet.

%%
% Add the new wavelet family to the stack of wavelet families.
wavemngr('add',familyName,familyShortName,familyWaveType, ...
    familyNums,fileWaveName)

%%
% Verify that the wavelet family is created.
wavemngr('read')

%%
% Display the two pairs of scaling and wavelet functions.
wname = 'mywb';
clf; wavefun(wname,'plot',7);

%%
% We can now use this new biorthogonal wavelet to analyze a signal.
load noisdopp; % Load a signal
x = noisdopp;
scales = 1:1:128;
coefs = cwt(x,scales,wname); 
clf; wscalogram('image',coefs,'scales',scales,'ydata',x);

%%
% This new biorthogonal wavelet can also be used to analyze an image.
load mask; % Load an image
[cA,cH,cV,cD] = dwt2(X,wname); 
clf; imagesc(abs([cA cH ; cV cD]));
clear cA cH cV cD

%%
% Delete the new wavelet. 
wavemngr('del',familyShortName);

%% Create a Wavelet of Type 3
% We now create a wavelet family of type 3, which is a clone of the Meyer 
% wavelet family. This family contains only one wavelet of type 3 and the 
% file used to compute the wavelet and scaling functions is *meyer.m*.
%
% *NOTE* 
% For wavelets of type 3, 4 or 5, we must define the bounds of the
% effective support, which is the set of real numbers where the 
% wavelet function is significantly different from 0. This is the 
% additional parameter required for these wavelet types.

%%
% Six parameters are associated with this wavelet family and must be passed 
% to the |wavemngr| function.
familyName      = 'MyWAVE T3';
familyShortName = 'mywc';
familyWaveType  = 3;
familyNums      = '';
fileWaveName    = 'meyer';
familyBounds    = [-4 4];

%%
% Add this new wavelet family to the stack of wavelet families.
wavemngr('add',familyName,familyShortName,familyWaveType, ...
    familyNums,fileWaveName,familyBounds)

%%
% Verify that the wavelet family is created.
wavemngr('read')

%%
% Display the pair of scaling and wavelet functions.
wname = 'mywc';
wavefun(wname,'plot',7);

%%
% Delete the new wavelet. 
wavemngr('del',familyShortName);

%% Create a Wavelet of Type 4
% Next we create a wavelet family of type 4. These wavelets have no 
% associated scaling function and so, are only available for continuous
% wavelet analysis.
familyName      = 'MyWAVE T4';
familyShortName = 'myw';
familyWaveType  = 4;
familyNums      = '3 4 5';
fileWaveName    = 'mywwavf';
familyBounds    = [0 1];

%%
% Add this new wavelet family to the stack of wavelet families.
wavemngr('add',familyName,familyShortName,familyWaveType, ...
    familyNums,fileWaveName,familyBounds)

%%
% The file *mywwavf.m* is already defined, so we can display its content.
type('mywwavf.m')
%%
% Note that this function is available in the directory 
% |$MATLABHOME/wavelet/wavedemo| of Wavelet Toolbox. 

%%
% Verify that the wavelet family is created.
wavemngr('read')

%%
% Display the wavelet function.
wname = 'myw3';
clf ; wavefun(wname,'plot',7);

%%
% Now compute the scalogram of the |sumsin| signal using this new
% wavelet.
%%
load sumsin; x = sumsin;
scales = 1:1:128;
coefs = cwt(x,scales,wname); 
clf; colormap(jet(128));
wscalogram('image',coefs,'scales',scales,'ydata',x);

%%
% We can analyze the same signal using the Continuous Wavelet Transform 
% (CWT) with the new wavelet and display it as a 3-D representation.
clf; coefs = cwt(x,scales,wname,'3Dplot');

%%
% Delete the new wavelet. 
wavemngr('del',familyShortName);

%%
% Delete the files created in this example and return in the initial
% directory. 
delete('mywa.mat')
delete('mywb.mat')
delete('wavelets.asc')
delete('wavelets.inf')
delete('wavelets.prv')
cd(Current_DIR);        % Return to the initial directory.

%% Summary
% This example has shown you how to define new wavelet families and make
% them available for use with functions and GUIs that operate on wavelet
% families in Wavelet Toolbox.

displayEndOfDemoMessage(mfilename)