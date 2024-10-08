run preamble;


%%

meanSB = [0.1569; 0.1543; 0.1542; 0.1564];
meanST = [0.1472; 0.1449; 0.1449; 0.1454];
meanRSST = [0.1766; 0.1740; 0.1713; 0.1709];

run2runMeanVals = [meanSB, meanST, meanRSST];

minMeanVals = min(run2runMeanVals);
meanMeanVals = mean(run2runMeanVals);
maxMeanVals = max(run2runMeanVals);

%%%

rmsSB = [0.2091; 0.2052; 0.2053; 0.2064];
rmsST = [0.2010; 0.1971; 0.1971; 0.1968];
rmsRSST = [0.2263; 0.2214; 0.2165; 0.2171];

run2runRMSvals = [rmsSB, rmsST, rmsRSST];

minRMSvals = min(run2runRMSvals);
meanRMSvals = mean(run2runRMSvals);
maxRMSvals = max(run2runRMSvals);


%%

% Initialise Figure
fig = fig + 1;
figName = 'Run_2_Run_Variation_Mean';
set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
             'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
pause(0.5);
hold on;
set(gca, 'positionConstraint', 'outerPosition', 'plotBoxAspectRatio', [1, 0.75, 0.75], ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 22, 'layer', 'top');

% Plot
for i = 1:width(minMeanVals)
    plot([i, i], [minMeanVals(i), maxMeanVals(i)], 'lineStyle', '-', 'lineWidth', 2, 'color', graphColours(1))
end

plot(minMeanVals, 'lineStyle', 'none', 'lineWidth', 2, 'marker', '_', 'markerSize', 20, 'color', graphColours(1));
plot(maxMeanVals, 'lineStyle', 'none', 'lineWidth', 2, 'marker', '_', 'markerSize', 20, 'color', graphColours(1));

% Format Figure
title('{-----}', 'interpreter', 'latex');
subtitle('{ }');
axis on;
box on;
grid off;
xlim([0; 4]);
ylim([0.12; 0.2]);
tickData = (1:3);
xticks(tickData);
tickData = (0.14:0.02:0.18);
yticks(tickData);
xticklabels({'\textit{Config~A}', '\textit{Config~B}', '\textit{Config~C}'});
set(gca, 'tickLabelInterpreter', 'latex');
ytickformat('%.3f');
xlabel({'{-----}'}, 'interpreter', 'latex');
ylabel({'{-----}'; '{$\overline{\overline{\varrho_{_{n}}}}$}'}, 'interpreter', 'latex');
tightInset = get(gca, 'TightInset');
set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                           (tightInset(2) + 0.00625), ...
                           (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                           (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
pause(0.5);
hold off;

% Save Figure
print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');


%%

% Initialise Figure
fig = fig + 1;
figName = 'Run_2_Run_Variation_RMS';
set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
             'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
pause(0.5);
hold on;
set(gca, 'positionConstraint', 'outerPosition', 'plotBoxAspectRatio', [1, 0.75, 0.75], ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 22, 'layer', 'top');

% Plot
for i = 1:width(minRMSvals)
    plot([i, i], [minRMSvals(i), maxRMSvals(i)], 'lineStyle', '-', 'lineWidth', 2, 'color', graphColours(1))
end

plot(minRMSvals, 'lineStyle', 'none', 'lineWidth', 2, 'marker', '_', 'markerSize', 20, 'color', graphColours(1));
plot(maxRMSvals, 'lineStyle', 'none', 'lineWidth', 2, 'marker', '_', 'markerSize', 20, 'color', graphColours(1));

% Format Figure
title('{-----}', 'interpreter', 'latex');
subtitle('{ }');
axis on;
box on;
grid off;
xlim([0; 4]);
ylim([0.18; 0.24]);
tickData = (1:3);
xticks(tickData);
tickData = (0.195:0.015:0.225);
yticks(tickData);
xticklabels({'\textit{Config~A}', '\textit{Config~B}', '\textit{Config~C}'});
ytickformat('%.3f');
xlabel({'{-----}'}, 'interpreter', 'latex');
ylabel({'{-----}'; '{$\overline{\mathrm{RMS}(\varrho_{_{n}}'')}$}'}, 'interpreter', 'latex');
set(gca, 'tickLabelInterpreter', 'latex')
tightInset = get(gca, 'TightInset');
set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                           (tightInset(2) + 0.00625), ...
                           (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                           (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
pause(0.5);
hold off;

% Save Figure
print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');