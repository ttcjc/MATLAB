CLA;

run preamble;


%%

E_QS = 0.046304987335294;
E_FS = 0.004922468583662;

content = importdata('/mnt/Processing/Data/Numerical/ParaView/Windsor_Upstream_2023/Windsor_SB_wW_Upstream_SC/underbody_Velocity_Profile.csv');
QS.xPos = content.data(:,2) / 1.044;
QS.U = content.data(:,1) / (40 / ((1 - E_QS)^2));

content = importdata('/mnt/Processing/Data/Numerical/ParaView/Windsor_fullScale/Windsor_SB_fullScale_multiPhase_uncoupled/underbody_Velocity_Profile.csv');
FS.xPos = content.data(:,2) / 4.176;
FS.U = content.data(:,1) / (22.2222 / ((1 - E_QS)^2));

clear content;


%%

% Initialise Figure
fig = fig + 1;
figName = 'Underbody_Velocity_Profile';
set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
                 'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
pause(0.5);
hold on;
set(gca, 'positionConstraint', 'outerPosition', 'plotBoxAspectRatio', [1, 0.75, 0.75], ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 22, 'layer', 'top');

% Plot Mean Profiles
plot(QS.xPos, QS.U, 'color', graphColours(1), 'lineWidth', 2);
plot(FS.xPos, FS.U, 'color', graphColours(4), 'lineWidth', 2);
% xline(1 - (0.1645 / 1.044))
% xline(1 - ((0.1645 + 0.6375) / 1.044))

% Format Figure
title('{-----}', 'interpreter', 'latex');
subtitle('{ }');
axis on;
box on;
grid off;
xlim([-0.1; 1.1]);
ylim([0.76; 1.16]);
tickData = (0.14:0.24:0.86);
xticks(tickData);
tickData = (0.84:0.08:1.08);
yticks(tickData);
xtickformat('%.2f');
ytickformat('%.2f');
xlabel({'{$x_{_{\ell}}$}'; '{-----}'}, 'interpreter', 'latex');
ylabel({'{-----}'; '{$\overline{|u|}_{_{n}}$}'}, 'interpreter', 'latex');
legend({'Reduced-Scale', ...
        'Full-Scale'}, ...
       'location', 'northEast', 'orientation', 'vertical', 'interpreter', 'latex', ...
       'fontSize', 18, 'box', 'off');
tightInset = get(gca, 'TightInset');
set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                           (tightInset(2) + 0.00625), ...
                           (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                           (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
pause(0.5);
hold off;

print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');