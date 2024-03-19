function savefig(filename,textSize,titleSize,res,Title,type,inches)
% saveFig(filename,textFontSize,titleFontSize,resolution,Title,type,inches)
% type = 4 for png, 1 = eps


warning('OFF', 'all')
orient portrait
     
% SO YOU DONT LOSE AN AXIS!
subplot('Position', [.98 .98 .01 .01]),axis off

if ~isempty(textSize)
    all_axis = findall(gcf, 'Type', 'Axes');
    set(all_axis, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', ...
        'FontSize', textSize, 'XGrid', 'off', 'YGrid', 'off');
end
if ~isempty(titleSize)
    all_text = findall(gcf, 'Type', 'Text');
    set(all_text, 'FontSize', titleSize, 'FontName', 'Arial');
end

if ~isempty(Title)
    subplot('Position', [.1 .97 .8 .01])
    text(.25, 1, Title, 'FontSize', round(titleSize))
end; axis off; box off;
  
if isempty(inches)
    inches = [17 10];
end

% SET PAPER SIZE
set(gcf, 'papersize', [inches(1) inches(2)])
set(gcf, 'paperposition', [.25 .25 inches(1)-.5 inches(2)-.5])

all_axis = findall(gcf, 'Type', 'Axes'); set(all_axis, 'Linewidth', 3); % all axes

if type == 1
    theFilename = strcat(filename, '.eps');
    eval(['print -f1 -r', num2str(res), ' -loose -dpsc2 -painters ''', theFilename, ''';']);
elseif type == 2
    theFilename = strcat(filename, '.bmp');
    eval(['print -r', num2str(res), ' -loose -dbmp -noui ''', theFilename, ''';']);
elseif type == 3
    theFilename = strcat(filename, '.tif');
    print('-dtiff', ['-r' num2str(res)], theFilename)
elseif type == 4
    theFilename = strcat(filename, '.png');
    eval(['print -r', num2str(res), ' -loose -dpng -noui ''', theFilename, ''';']);
elseif type == 5
    theFilename = strcat(filename, '.jpg');
    eval(['print -r', num2str(res), ' -loose -djpeg -noui ''', theFilename, ''';']);
elseif type == 20
    theFilename = strcat(filename, '.fig');
    savefig(theFilename);
elseif type == 10
    theFilename = strcat(filename, '.pdf');
    eval(['print -r', num2str(res), ' -loose -dpdf -painters ''', theFilename, ''';']);
end

