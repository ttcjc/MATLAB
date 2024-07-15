run preamble;


%%

content = importdata('/mnt/Processing/Data/Numerical/ParaView/Windsor_Upstream_2023/Windsor_SB_wW_Upstream_SC/centreline_Pressure_Distribution.csv');
QS.positionData = content.data(:,[1,3]);
QS.p = content.data(:,4);

content = importdata('/mnt/Processing/Data/Numerical/ParaView/Windsor_fullScale/Windsor_SB_fullScale_multiPhase_uncoupled/centreline_Pressure_Distribution.csv');
FS.positionData = content.data(:,[1,3]);
FS.p = content.data(:,4);

clear content;

%%

QS.positionData(:,1) = (QS.positionData(:,1) + 0.56075) + 1.325;
QS.positionData = QS.positionData / 1.044;

FS.positionData(:,1) = FS.positionData(:,1) + 2.243;
FS.positionData = FS.positionData / 4.176;


%%

U = 40;
rho = 1.269;
pRef = 19.524 * rho;
QS.p = QS.p * rho;
QS.Cp = (QS.p - pRef) / (0.5 * rho * U^2);
            
U = 22.222222222222222;
rho = 1.269;
pRef = 0 * rho;
FS.p = FS.p * rho;
FS.Cp = (FS.p - pRef) / (0.5 * rho * U^2);

clear U rho pRef;


%%

Am = (0.289 * 0.389) + (2 * (0.046 * 0.055));
At = (2 * (0.9519083 + (3.283 * tan(atan(0.0262223 / 9.44)))) * 1.32);
QS.p = (QS.p + (2 * (Am / At))) / (1 + (2 * (Am / At)));
QS.Cp = (QS.Cp + (2 * (Am / At))) / (1 + (2 * (Am / At)));

Am = ((4 * 0.289) * (4 * 0.389)) + (2 * (((4 * 0.05) - 0.018) * (4 * 0.055)));
At = 14.336 * 26.624;
FS.p = (FS.p + (2 * (Am / At))) / (1 + (2 * (Am / At)));
FS.Cp = (FS.Cp + (2 * (Am / At))) / (1 + (2 * (Am / At)));

clear Am At;


%%

index = QS.positionData(:,2) < (0.1 / 1.044) & QS.positionData(:,1) > (1e-3 / 1.044) & QS.positionData(:,1) < (1.043 / 1.044);
QS.positionData = QS.positionData(index,:);
QS.p = QS.p(index);
QS.Cp = QS.Cp(index);

index = FS.positionData(:,2) < (0.1 / 1.044) & FS.positionData(:,1) > (1e-3 / 1.044) & FS.positionData(:,1) < (1.043 / 1.044);
FS.positionData = FS.positionData(index,:);
FS.p = FS.p(index);
FS.Cp = FS.Cp(index);

clear index;


%%

[QS.positionData, index] = sortrows(QS.positionData, 1);
QS.p = QS.p(index);
QS.Cp = QS.Cp(index);

[FS.positionData, index] = sortrows(FS.positionData, 1);
FS.p = FS.p(index);
FS.Cp = FS.Cp(index);

clear index;

%%

% Initialise Figure
fig = fig + 1;
figName = 'Cp';
set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
                 'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
pause(0.5);
hold on;
set(gca, 'positionConstraint', 'outerPosition', 'plotBoxAspectRatio', [1, 0.75, 0.75], ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 20, 'layer', 'top');

% Plot Mean Profiles
plot(QS.positionData(:,1), QS.Cp, 'color', graphColours(1), 'lineWidth', 2);
plot(FS.positionData(:,1), FS.Cp, 'color', graphColours(2), 'lineWidth', 2);

% Format Figure
title('{-----}', 'interpreter', 'latex');
subtitle('{ }');
axis on;
box on;
grid off;
xlim([-0.1; 1.1]);
ylim([-1.5; 1]);
tickData = (0.14:0.24:0.86);
xticks(tickData);
tickData = (-1:0.5:0.5);
yticks(tickData);
xtickformat('%.2f');
ytickformat('%.1f');
xlabel({'{$x_{_{\ell}}$}'; '{-----}'}, 'interpreter', 'latex');
ylabel({'{-----}'; '{$\overline{C_{_{p}}}$}'}, 'interpreter', 'latex');
legend({'Reduced-Scale', ...
        'Full-Scale'}, ...
       'location', 'northEast', 'orientation', 'vertical', 'interpreter', 'latex', ...
       'fontSize', 16, 'box', 'off');
tightInset = get(gca, 'TightInset');
set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                           (tightInset(2) + 0.00625), ...
                           (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                           (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
pause(0.5);
hold off;

% print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');