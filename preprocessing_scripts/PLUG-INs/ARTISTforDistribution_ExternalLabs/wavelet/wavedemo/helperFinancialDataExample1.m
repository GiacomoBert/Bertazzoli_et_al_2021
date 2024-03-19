function helperFinancialDataExample1(data,time,str)
%   This function helperFinancialDataExample1 is only in support of
%   FinancialDataExample. It may change in a future release.

plot(data)
set(gca,'xtick',1:40:260);
set(gca,'xticklabel',time(1:40:250))
set(gca,'xlim',[1 260])
xp = [146   260   260   146];
AX = gca;
yp = repelem(AX.YLim,2);
patch(xp,yp,[0.5 0.5 0.5],'facealpha',0.1);
title(str)
xlabel('Year')