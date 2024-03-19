function [] = plotTriNoise(est_noise,timeAxis,data)
% This function plots an interactive view showing the DDWiener estimate of
% the noise distribution at each channel and in each trial. By clicking at
% any specific point in the color map you can visualize the ERP for that 
% specific trial (vs. other trials) on that specific channel. You can end 
% the interactive view by pressing any key on the keyboard (while having 
% the figure active).
%
% .........................................................................
% 20 March 2017 : Tuomas Mutanen, NBE, Aalto university  
% .........................................................................

% Plotting the color map showing the noise estimate in each trial and
% channel

figure;
subplot(1,2,1)
imagesc(log10(est_noise));
colormap(hot);
c = colorbar;
xlabel('Trials');
ylabel('Channels');
title('The DDWiener estimate of the noise distribution accross trials and channels')
title(c,'\muV');

dch = datacursormode(gcf);
set(dch,'enable','on','displaystyle','datatip','snaptodatavertex','on');
set(dch,'UpdateFcn',@UpdateCursor);

disp('Interactive trial-noise view started. To visualize any trial, click the corresponding pixel in the noise map.')
disp('To end the interactive trial-noise view, activate the corresponding figure and press any key.')

while 1 
    w = waitforbuttonpress;
    if w == 1
        disp('Interactive trial-noise view ended')
        break;
    end
    pause(0.1);
    pos=get(gcf,'userdata');
    
    subplot(1,2,2);
    plot(timeAxis,squeeze(data(pos(2),:,:)),'k')
    hold on;
    plot(timeAxis,data(pos(2),:,pos(1)),'r','LineWidth',2);
    xlim([min(timeAxis), max(timeAxis)])
    title(['Ch:',num2str(pos(2)),', Trial ',num2str(pos(1)),' vs. other trials.'])
    xlabel('Time [ms]')
    ylabel('V [\muV]');
    hold off;

end

function txt=UpdateCursor(empt,event_obj)

    %   A nested function that returns the coordinates of the clicked pixel 

pos = get(event_obj,'Position');
set(gcf,'userdata',pos);
txt={''};

