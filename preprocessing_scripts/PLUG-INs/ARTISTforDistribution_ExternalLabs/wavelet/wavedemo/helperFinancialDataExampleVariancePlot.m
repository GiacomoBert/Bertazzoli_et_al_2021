function helperFinancialDataExampleVariancePlot(var1,var2,varargin)
%   This function helperFinancialDataExampleVariancePlot is only in support
%   of FinancialDataExample. It may change in a future release.

    if strncmpi(varargin{1},'table',1)
        var1 = table2array(var1);
        var1 = var1(1:end-1,2:end);
        var2 = table2array(var2);
        var2 = var2(1:end-1,2:end);
    else
        var1 = var1(1:end-1,:);
        var2 = var2(1:end-1,:);
    end
V1 = zeros(size(var1));
V2 = zeros(size(var2));

V1(:,1) = var1(:,2)-var1(:,1);
V1(:,3) = var1(:,3)-var1(:,2);
V2(:,1) = var2(:,2)-var2(:,1);
V2(:,3) = var2(:,3)-var2(:,2);
errorbar(1:size(var1,1),var1(:,2),V1(:,1),V1(:,3),'ko',...
    'markerfacecolor',[0 0 0]);
hold on;
errorbar(1.5:size(var1,1)+0.5,var2(:,2),V2(:,1),V2(:,3),'b^',...
    'markerfacecolor',[0 0 1]);
set(gca,'xtick',1.25:size(var1,1)+0.25);
set(gca,'xticklabel',varargin{2});
grid on;