run preamble;

caseFolders = dir('~/Downloads/experimentalSpray/OGI/');

caseFolders(1:3) = [];

testData = [];

for i = 1:height(caseFolders)
    testFiles = dir(['~/Downloads/experimentalSpray/OGI/', caseFolders(i).name]);
    testFiles(1:2) = [];
    
    for j = 1:height(testFiles)
        fileID = fopen(['~/Downloads/experimentalSpray/OGI/', caseFolders(i).name, '/' testFiles(j).name]);
        content = textscan(fileID, '%s', 'headerLines', 1, 'delimiter', '\n');
        fclose(fileID);
        
        for k = 1:height(content{1})
            testData = [testData; split(content{1}{k}, ',')']; %#ok<AGROW>
        end
        
    end
    
end

testData(:,3) = [];
testData = sortrows(testData, [1,3]);

testData = [testData, cell(height(testData),2)];


for i = 1:height(testData)
    
    if strcmp(testData{i,1}, '15/08/22')
        
        if strcmp(testData{i,3}(1:2), '14')
            testData{i,7} = '15.5';
        elseif strcmp(testData{i,3}(1:2), '15')
            testData{i,7} = '15.5';
        elseif strcmp(testData{i,3}(1:2), '16')
            testData{i,7} = '13.0';
        end
        
    elseif strcmp(testData{i,1}, '16/08/22')
        
        if strcmp(testData{i,3}(1:2), '09')
            testData{i,7} = '14.0';
        elseif strcmp(testData{i,3}(1:2), '11')
            testData{i,7} = '14.5';
        end
        
    elseif strcmp(testData{i,1}, '17/08/22')
        
        if strcmp(testData{i,3}(1:2), '10')
            testData{i,7} = '13.0';
        elseif strcmp(testData{i,3}(1:2), '11')
            testData{i,7} = '13.5';
        elseif strcmp(testData{i,3}(1:2), '13')
            testData{i,7} = '14.0';
        elseif strcmp(testData{i,3}(1:2), '14')
            testData{i,7} = '13.0';
        elseif strcmp(testData{i,3}(1:2), '15')
            testData{i,7} = '14.0';
        end
        
    elseif strcmp(testData{i,1}, '18/08/22')
        
        if strcmp(testData{i,3}(1:2), '12')
            testData{i,7} = '13.5';
        elseif strcmp(testData{i,3}(1:2), '13')
            testData{i,7} = '13.5';
        elseif strcmp(testData{i,3}(1:2), '14')
            testData{i,7} = '12.5';
        elseif strcmp(testData{i,3}(1:2), '16')
            testData{i,7} = '13.5';
        end
        
    elseif strcmp(testData{i,1}, '19/08/22')
        
        if strcmp(testData{i,3}(1:2), '10')
            testData{i,7} = '07.5';
        elseif strcmp(testData{i,3}(1:2), '11')
            testData{i,7} = '06.0';
        end
        
    end
    
end

T = str2num(cell2mat(testData(:,5))); %#ok<ST2NM>
dp = str2num(cell2mat(testData(:,7))); %#ok<ST2NM>
phi = zeros([height(testData),1]);


for i = 1:height(testData)
    phi(i) = 100 * (exp((17.625 * dp(i)) / (243.04 + dp(i))) / exp((17.625 * T(i)) / (243.04 + T(i))));
    testData{i,8} = num2str(phi(i));
end

index15 = find(strcmp(testData(:,1), '15/08/22'));
min15 = min(phi(index15));
mean15 = mean(phi(index15));
max15 = max(phi(index15));

index16 = find(strcmp(testData(:,1), '16/08/22'));
min16 = min(phi(index16));
mean16 = mean(phi(index16));
max16 = max(phi(index16));

index17 = find(strcmp(testData(:,1), '17/08/22'));
min17 = min(phi(index17));
mean17 = mean(phi(index17));
max17 = max(phi(index17));

index18 = find(strcmp(testData(:,1), '18/08/22'));
min18 = min(phi(index18));
mean18 = mean(phi(index18));
max18 = max(phi(index18));

index19 = find(strcmp(testData(:,1), '19/08/22'));
min19 = min(phi(index19));
mean19 = mean(phi(index19));
max19 = max(phi(index19));

minVals = [min15; min16; min17; min18; min19];
meanVals = [mean15; mean16; mean17; mean18; mean19];
maxVals = [max15; max16; max17; max18; max19];


% Initialise Figure
fig = fig + 1;
figName = 'Relative_Humidity_Variation';
set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
             'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
pause(0.5);
hold on;
set(gca, 'positionConstraint', 'outerPosition', 'plotBoxAspectRatio', [1, 0.75, 0.75], ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 22, 'layer', 'top');

% Plot
for i = 1:height(minVals)
    plot([i, i], [minVals(i), maxVals(i)], 'lineStyle', '-', 'lineWidth', 2, 'color', graphColours(1))
end

plot(minVals, 'lineStyle', 'none', 'lineWidth', 2, 'marker', '_', 'markerSize', 20, 'color', graphColours(1));
% plot(meanVals, 'lineStyle', 'none', 'lineWidth', 2, 'marker', '_', 'markerSize', 10, 'color', graphColours(1));
plot(maxVals, 'lineStyle', 'none', 'lineWidth', 2, 'marker', '_', 'markerSize', 20, 'color', graphColours(1));

% Format Figure
title('{-----}', 'interpreter', 'latex');
subtitle('{ }');
axis on;
box on;
grid off;
xlim([0; 6]);
ylim([0; 100]);
xticks(1:1:5);
yticks(20:20:80);
% xtickformat('%+.2g');
% ytickformat('%+.2g');
% ztickformat('%+.2g');
xlabel({'{Test Period}'; '{-----}'}, 'interpreter', 'latex');
ylabel({'{-----}'; '{$\varphi$ $(\%)$}'}, 'interpreter', 'latex');
tightInset = get(gca, 'TightInset');
set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                           (tightInset(2) + 0.00625), ...
                           (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                           (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
pause(0.5);
hold off;

% Save Figure
print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');