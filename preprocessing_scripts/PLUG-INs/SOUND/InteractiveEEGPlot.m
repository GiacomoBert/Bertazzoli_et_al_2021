function InteractiveEEGPlot(timeAxis,data1,data2,chanlocs,Y_Lims)
% This function creates an interactive plot of the ERPs. The time course of
% a specific channel can be visualized by clicking on its trace, producing
% a pop-up sub-axis view.
%
% .........................................................................
% 20 March 2017 : Tuomas Mutanen, NBE, Aalto university  
% .........................................................................

%   Making the EEG plot including all the channels

if data2
    if nargin < 5
        Y_Lims = [-max(max(abs(data2))), max(max(abs(data2)))];
    end
    
figure;
    
    for ii = 1:size(data1)
        subplot(4,5,ii)
        plot(timeAxis,data1(ii,:),'r',timeAxis,data2(ii,:),'k', 'ButtonDownFcn', {@newFigure1,ii,timeAxis,data1,data2,chanlocs(ii).labels})
        hold on;
        plot([0,0],Y_Lims,':');
        set(gca,'Visible','off');
        xlim([timeAxis(1), timeAxis(end)])
        ylim(Y_Lims)

        set(gca,'xtick',[])
        set(gca,'ytick',[])
        set(gca,'yticklabel',[])
        strmax = [chanlocs(ii).labels];
        text(0,double(Y_Lims(2)),strmax,'HorizontalAlignment','right');
        
    end
    
    subplot(4,5,20);
    plot(0,0);
    xlim([timeAxis(1), timeAxis(end)])
    ylim(Y_Lims)

else
    if nargin < 5
        Y_Lims = [-max(max(abs(data1))), max(max(abs(data1)))];
    end   

         
figure;
    for ii = 1:size(data1)
        subplot(4,5,ii)
        plot(timeAxis,data1(ii,:),'r', 'ButtonDownFcn', {@newFigure1,ii,timeAxis,data1,data2,chanlocs(ii).labels})
        hold on;
        plot([0,0],Y_Lims,':');
        set(gca,'Visible','off');
        xlim([timeAxis(1), timeAxis(end)])
        ylim(Y_Lims)

        set(gca,'xtick',[])
        set(gca,'ytick',[])
        set(gca,'yticklabel',[])
       
        strmax = [chanlocs(ii).labels];
        text(0,Y_Lims(2),strmax,'HorizontalAlignment','right');
        
    end
    
    subplot(4,5,20);
    plot(0,0);
    xlim([timeAxis(1), timeAxis(end)])
    ylim(Y_Lims)
    
end

    xlabel('Time [ms]')
    ylabel('V [\muV]')

end

function newFigure1(src,evt,ii,timeAxis,data1,data2,channame)

    %   A nested function that creates pop-up magnifications of the
    %   user-chosen channels. 

    if data2
    figure;
    p1 = plot(timeAxis,data1(ii,:),'r','LineWidth',1);
    hold on;
    p2 = plot(timeAxis,data2(ii,:),'k','LineWidth',1);
    title(channame)
    xlabel('Time [ms]')
    ylabel('V [\muV]')
    legend([p1 p2],['Original signal in Ch ',channame],['Signal in Ch ',channame, ' after SOUND'])
    
    else
    figure;
    plot(timeAxis,data1(ii,:),'r','LineWidth',1);
    title(channame)
    xlabel('Time [ms]')
    ylabel('V [\muV]')
    end
end   

