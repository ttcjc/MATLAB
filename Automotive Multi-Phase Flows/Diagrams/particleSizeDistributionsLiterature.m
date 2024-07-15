%% Literature Particle Size Distribution Plots v1.0
% ----
% Lorem ipsum


%% Preamble

run preamble;

disp('================================================');
disp('Literature Particle Size Distribution Plots v1.0');
disp('================================================');

disp(' ');
disp(' ');


%% Changelog

% v1.0 - Initial Commit


%% Shearman

shearman = importdata('~/Data/Particle Size Data/Shearman.csv');

% Initialise Figure
fig = fig + 1;
figName = 'Shearman';
set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
             'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
pause(0.5);
hold on;
set(gca, 'positionConstraint', 'outerPosition', 'plotBoxAspectRatio', [1, 0.75, 0.75], ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 22, 'layer', 'top');

% Plot Distribution
plot(shearman(:,1), shearman(:,2), 'color', graphColours(1), 'lineStyle', '-', 'lineWidth', 2, ...
                                   'marker', 'o', 'markerSize', 7.5, 'markerFaceColor', graphColours(1));

% Format Figure
title('{-----}', 'interpreter', 'latex');
subtitle('{ }');
axis on;
box on;
grid off;
xlim([0; 600]);
ylim([0; 60]);
tickData = (120:120:480);
xticks(tickData);
tickData = ([]);
yticks(tickData);
xlabel({'{$D_{_{p}}$ ($\mu m$)}'; '{-----}'}, 'interpreter', 'latex');
ylabel({'{-----}'; '{Population}'}, 'interpreter', 'latex');
tightInset = get(gca, 'TightInset');
set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                           (tightInset(2) + 0.00625), ...
                           (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                           (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
pause(0.5);
hold off;

% Save Figure
print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');


%% Borg

borg25 = importdata('~/Data/Particle Size Data/Borg_25.csv');
borg25 = [(20:40:460)', (borg25.data * 100)];
borg50 = importdata('~/Data/Particle Size Data/Borg_50.csv');
borg50 = [(20:40:460)', (borg50.data * 100)];

% Initialise Figure
fig = fig + 1;
figName = 'Borg';
set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
             'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
pause(0.5);
hold on;
set(gca, 'positionConstraint', 'outerPosition', 'plotBoxAspectRatio', [1, 0.75, 0.75], ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 22, 'layer', 'top');

% Plot Distribution
plot(borg25(:,1), borg25(:,2), 'color', graphColours(1), 'lineStyle', '-', 'lineWidth', 2, ...
                               'marker', 'o', 'markerSize', 7.5, 'markerFaceColor', graphColours(1));
plot(borg50(:,1), borg50(:,2), 'color', graphColours(2), 'lineStyle', '-', 'lineWidth', 2, ...
                               'marker', 'o', 'markerSize', 7.5, 'markerFaceColor', graphColours(2));

% Format Figure
title('{-----}', 'interpreter', 'latex');
subtitle('{ }');
axis on;
box on;
grid off;
xlim([0; 600]);
ylim([0; 70]);
tickData = (120:120:480);
xticks(tickData);
tickData = ([]);
yticks(tickData);
xlabel({'{$D_{_{p}}$ ($\mu m$)}'; '{-----}'}, 'interpreter', 'latex');
ylabel({'{-----}'; '{Population}'}, 'interpreter', 'latex');
legend({'$25\,m$ Downstream', ...
        '$50\,m$ Downstream'}, 'location', 'northEast', ...
                               'orientation', 'vertical', 'interpreter', 'latex', ...
                               'fontSize', 16, 'box', 'off');
tightInset = get(gca, 'TightInset');
set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                           (tightInset(2) + 0.00625), ...
                           (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                           (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
pause(0.5);
hold off;

% Save Figure
print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');


%% Bouchet

bouchet40 = importdata('~/Data/Particle Size Data/Bouchet Data/40 kph.csv');
bouchet40(:,1) = bouchet40(:,1) * 1e3;
bouchet80 = importdata('~/Data/Particle Size Data/Bouchet Data/80 kph.csv');
bouchet80(:,1) = bouchet80(:,1) * 1e3;
bouchet140 = importdata('~/Data/Particle Size Data/Bouchet Data/140 kph.csv');
bouchet140(:,1) = bouchet140(:,1) * 1e3;

% Initialise Figure
fig = fig + 1;
figName = 'Bouchet';
set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
             'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
pause(0.5);
hold on;
set(gca, 'positionConstraint', 'outerPosition', 'plotBoxAspectRatio', [1, 0.75, 0.75], ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 22, 'layer', 'top');

% Plot Distribution
plot(bouchet40(:,1), bouchet40(:,2), 'color', graphColours(1), 'lineStyle', '-', 'lineWidth', 2);
plot(bouchet80(:,1), bouchet80(:,2), 'color', graphColours(2), 'lineStyle', '-', 'lineWidth', 2);
plot(bouchet140(:,1), bouchet140(:,2), 'color', graphColours(3), 'lineStyle', '-', 'lineWidth', 2);

% Format Figure
title('{-----}', 'interpreter', 'latex');
subtitle('{ }');
axis on;
box on;
grid off;
xlim([0; 600]);
ylim([0; 11]);
tickData = (120:120:480);
xticks(tickData);
tickData = ([]);
yticks(tickData);
xlabel({'{$D_{_{p}}$ ($\mu m$)}'; '{-----}'}, 'interpreter', 'latex');
ylabel({'{-----}'; '{Population}'}, 'interpreter', 'latex');
legend({'$40\,km {\cdot} hr^{-1}$', ...
        '$80\,km {\cdot} hr^{-1}$', ...
        '$140\,km {\cdot} hr^{-1}$'}, 'location', 'northEast', ...
                                      'orientation', 'vertical', 'interpreter', 'latex', ...
                                      'fontSize', 16, 'box', 'off');
tightInset = get(gca, 'TightInset');
set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                           (tightInset(2) + 0.00625), ...
                           (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                           (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
pause(0.5);
hold off;

% Save Figure
print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');


%% Strohbucker

strohbucker40 = importdata('~/Data/Particle Size Data/Strohbucker_40.csv');
strohbucker40 = [(strohbucker40(:,1) * 1e3), (strohbucker40(:,2) * 100)];
strohbucker80 = importdata('~/Data/Particle Size Data/Strohbucker_80.csv');
strohbucker80 = [(strohbucker80(:,1) * 1e3), (strohbucker80(:,2) * 100)];
strohbucker120 = importdata('~/Data/Particle Size Data/Strohbucker_120.csv');
strohbucker120 = [(strohbucker120(:,1) * 1e3), (strohbucker120(:,2) * 100)];

% Initialise Figure
fig = fig + 1;
figName = 'Strohbucker';
set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
             'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
pause(0.5);
hold on;
set(gca, 'positionConstraint', 'outerPosition', 'plotBoxAspectRatio', [1, 0.75, 0.75], ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 22, 'layer', 'top');

% Plot Distribution
plot(strohbucker40(:,1), strohbucker40(:,2), 'color', graphColours(1), 'lineStyle', '-', 'lineWidth', 2);
plot(strohbucker80(:,1), strohbucker80(:,2), 'color', graphColours(2), 'lineStyle', '-', 'lineWidth', 2);
plot(strohbucker120(:,1), strohbucker120(:,2), 'color', graphColours(3), 'lineStyle', '-', 'lineWidth', 2);

% Format Figure
title('{-----}', 'interpreter', 'latex');
subtitle('{ }');
axis on;
box on;
grid off;
xlim([0; 600]);
ylim([0; 22]);
tickData = (120:120:480);
xticks(tickData);
tickData = ([]);
yticks(tickData);
xlabel({'{$D_{_{p}}$ ($\mu m$)}'; '{-----}'}, 'interpreter', 'latex');
ylabel({'{-----}'; '{Population}'}, 'interpreter', 'latex');
legend({'$40\,km {\cdot} hr^{-1}$', ...
        '$80\,km {\cdot} hr^{-1}$', ...
        '$120\,km {\cdot} hr^{-1}$'}, 'location', 'northEast', ...
                                      'orientation', 'vertical', 'interpreter', 'latex', ...
                                      'fontSize', 16, 'box', 'off');
tightInset = get(gca, 'TightInset');
set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                           (tightInset(2) + 0.00625), ...
                           (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                           (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
pause(0.5);
hold off;

% Save Figure
print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');