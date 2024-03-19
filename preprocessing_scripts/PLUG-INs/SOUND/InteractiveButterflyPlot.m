function InteractiveButterflyPlot(timeAxis,data1,data2,chanlocs,Y_Lims)

% This function creates an interactive Butterfly plot with the superimposed
% ERPs for all the channels. The original ERP signals are shown in red
% color and the cleaned ERP signals in black color. The topoplots of the 
% original and the cleaned ERP signals for a specific time point can be 
% visualized by clicking on that specific time point, producing a pop-up 
% view with the topoplots.
%
% .........................................................................
% 20 March 2017 : Tuomas Mutanen, NBE, Aalto university  
% .........................................................................


%   Making the butterfly plot

if data2
    if nargin < 5
        Y_Lims = [1.5*min(min(data2)), 1.5*max(max(data2))];
    end     
figure;

        plot(timeAxis,data1','r',timeAxis,data2','k', 'ButtonDownFcn', {@newFigure1,timeAxis,data1,data2,chanlocs})
        hold on;
        plot([0,0],Y_Lims,':','LineWidth',2);
        xlim([timeAxis(1), timeAxis(end)])
        ylim(Y_Lims)
        title('Butterflyplot');
        xlabel('Time [ms]')
        
        p1 = plot(timeAxis,data1(1,:),'r')
        p2 = plot(timeAxis,data2(1,:),'k')
        legend([p1 p2],'Original data','Data after SOUND')

else
    if nargin < 5
        Y_Lims = [1.5*min(min(data1)), 1.5*max(max(data1))];
    end   

         
figure;

        plot(timeAxis,data1','r', 'ButtonDownFcn', {@newFigure1,timeAxis,data1,data2,chanlocs})
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

end

function newFigure1(src,evt,timeAxis,data1,data2,chanlocs)
    
    %   A nested function that creates the pop-up-topoplot figure at the 
    %   wanted time points, when the user clicks the butterfly plots
    
    [~,PlotInd] = min(abs(timeAxis-evt.IntersectionPoint(1)));
    if data2
    figure;
    subplot(1,2,1)
    title(['Topography of the original data at time ',num2str(round(timeAxis(PlotInd))),' ms'])
    topoplot(data1(:,PlotInd), chanlocs); 
    subplot(1,2,2)
    title(['Topography of the cleaned data at time ',num2str(round(timeAxis(PlotInd))),' ms'])
    topoplot(data2(:,PlotInd), chanlocs); 
    else
    figure;
    title(['Topography of the original data at time ',num2str(round(timeAxis(PlotInd))),' ms'])
    topoplot(data1(:,PlotInd), chanlocs); 
    end
end   

