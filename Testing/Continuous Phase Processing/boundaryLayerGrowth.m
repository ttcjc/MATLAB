%% Boundary Layer Development Plots v1.0
% ----
% Lorem Ipsum


%% Preamble

run preamble;

DataLoc = '/home/lunet/ttcjc/Data/Boundary Layer Measurements (Luckhurst)';

figSave = false; % Save .fig File(s)

disp('=====================================');
disp('Boundary Layer Development Plots v1.0');
disp('=====================================');

disp(' ');
disp(' ');


%% Changelog

% v1.0 - Initial Commit


%% Load & Normalise Data

U = 40;
L = 1.044;

content = importdata([DataLoc, '/boundaryLayerLocationA_U.xy']);
crickmoreCFD_A = [(content(:,1) / L), (sqrt(content(:,2).^2 + content(:,3).^2 + content(:,4).^2) / U)];
clear content;

content = importdata([DataLoc, '/boundaryLayerLocationB_U.xy']);
crickmoreCFD_B = [(content(:,1) / L), (sqrt(content(:,2).^2 + content(:,3).^2 + content(:,4).^2) / U)];
clear content;

content = importdata([DataLoc, '/Luckhurst_CFD_Location_A.csv']);
luckhurstCFD_A = [((content(:,2) / 1000) / L), content(:,1)];
clear content;

content = importdata([DataLoc, '/Luckhurst_CFD_Location_B.csv']);
luckhurstCFD_B = [((content(:,2) / 1000) / L), content(:,1)];
clear content;

content = importdata([DataLoc, '/Luckhurst_Exp_Location_A.csv']);
luckhurstExp_A = [((content(:,2) / 1000) / L), content(:,1)];
clear content;

content = importdata([DataLoc, '/Luckhurst_Exp_Location_B.csv']);
luckhurstExp_B = [((content(:,2) / 1000) / L), content(:,1)];


%% Plot Profile A

% Initialise Figure
fig = fig + 1;
figName = 'Boundary_Layer_Profile_A';
set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
                 'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
pause(0.5);
hold on;
set(gca, 'positionConstraint', 'outerPosition', 'plotBoxAspectRatio', [1, 0.75, 0.75], ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 22, 'layer', 'top');

% Plot Mean Profiles
plot(crickmoreCFD_A(:,2), crickmoreCFD_A(:,1), 'color', graphColours(1), 'lineWidth', 2);
plot(luckhurstCFD_A(:,2), luckhurstCFD_A(:,1), 'color', graphColours(2), 'lineWidth', 2);
plot(luckhurstExp_A(:,2), luckhurstExp_A(:,1), 'color', graphColours(3), 'lineWidth', 2, ...
                                                                         'lineStyle', 'none', ...
                                                                         'marker', 'o', ...
                                                                         'markerSize', 5);

% Format Figure
title('{-----}', 'interpreter', 'latex');
subtitle('{ }');
axis on;
box on;
grid off;
xlim([0; 1.1]);
ylim([0; 0.08]);
tickData = (0.22:0.22:0.88);
xticks(tickData);
tickData = (0.016:0.016:0.064);
yticks(tickData);
xtickformat('%.2f');
ytickformat('%.3f');
xlabel({'{$U_{_{n}}$}'; '{-----}'}, 'interpreter', 'latex');
ylabel({'{-----}'; '{$z_{_{\ell}}$}'}, 'interpreter', 'latex');
legend({'Numerical (Crickmore)', 'Numerical (Luckhurst)', 'Experimental (Luckhurst)'}, ...
       'location', 'northWest', 'orientation', 'vertical', 'interpreter', 'latex', ...
       'fontSize', 16, 'box', 'off')
tightInset = get(gca, 'TightInset');
set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                           (tightInset(2) + 0.00625), ...
                           (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                           (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
pause(0.5);
hold off;

% Save Figure
print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');


%% Plot Profile B

% Initialise Figure
fig = fig + 1;
figName = 'Boundary_Layer_Profile_B';
set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
                 'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
pause(0.5);
hold on;
set(gca, 'positionConstraint', 'outerPosition', 'plotBoxAspectRatio', [1, 0.75, 0.75], ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 22, 'layer', 'top');

% Plot Mean Profiles
plot(crickmoreCFD_B(:,2), crickmoreCFD_B(:,1), 'color', graphColours(1), 'lineWidth', 2);
plot(luckhurstCFD_B(:,2), luckhurstCFD_B(:,1), 'color', graphColours(2), 'lineWidth', 2);
plot(luckhurstExp_B(:,2), luckhurstExp_B(:,1), 'color', graphColours(3), 'lineWidth', 2, ...
                                                                         'lineStyle', 'none', ...
                                                                         'marker', 'o', ...
                                                                         'markerSize', 5);

% Format Figure
title('{-----}', 'interpreter', 'latex');
subtitle('{ }');
axis on;
box on;
grid off;
xlim([0; 1.1]);
ylim([0; 0.08]);
tickData = (0.22:0.22:0.88);
xticks(tickData);
tickData = (0.016:0.016:0.064);
yticks(tickData);
xtickformat('%.2f');
ytickformat('%.3f');
xlabel({'{$U_{_{n}}$}'; '{-----}'}, 'interpreter', 'latex');
ylabel({'{-----}'; '{$z_{\ell}$}'}, 'interpreter', 'latex');
legend({'Numerical (Crickmore)', 'Numerical (Luckhurst)', 'Experimental (Luckhurst)'}, ...
       'location', 'northWest', 'orientation', 'vertical', 'interpreter', 'latex', ...
       'fontSize', 16, 'box', 'off')
tightInset = get(gca, 'TightInset');
set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                           (tightInset(2) + 0.00625), ...
                           (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                           (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
pause(0.5);
hold off;

% Save Figure
print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');