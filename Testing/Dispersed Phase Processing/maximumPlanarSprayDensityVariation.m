run preamble;

% QS_1L = load('/mnt/Processing/Data/PhD/Numerical/MATLAB/planarSprayMap/Windsor_Upstream_2023/Windsor_SB_wW_Upstream_SC/X_P0_20225/T12525_T40000_F400_D1_D147_cumulative.mat', 'mapData').mapData;
% QS_2L = load('/mnt/Processing/Data/PhD/Numerical/MATLAB/planarSprayMap/Windsor_Upstream_2023/Windsor_SB_wW_Upstream_SC/X_P1_24625/T12525_T40000_F400_D1_D147_cumulative.mat', 'mapData').mapData;

FS_Uncoupled_1L = load('/mnt/Processing/Data/PhD/Numerical/MATLAB/planarSprayMap/Windsor_fullScale/Windsor_SB_fullScale_multiPhase_uncoupled/X_P6_109/T1002_T3200_F50_D20_D400_cumulative.mat', 'mapData').mapData;
FS_Uncoupled_2L = load('/mnt/Processing/Data/PhD/Numerical/MATLAB/planarSprayMap/Windsor_fullScale/Windsor_SB_fullScale_multiPhase_uncoupled/X_P10_285/T1002_T3200_F50_D20_D400_cumulative.mat', 'mapData').mapData;
FS_Uncoupled_3L = load('/mnt/Processing/Data/PhD/Numerical/MATLAB/planarSprayMap/Windsor_fullScale/Windsor_SB_fullScale_multiPhase_uncoupled/X_P14_461/T1002_T3200_F50_D20_D400_cumulative.mat', 'mapData').mapData;
FS_Uncoupled_4L = load('/mnt/Processing/Data/PhD/Numerical/MATLAB/planarSprayMap/Windsor_fullScale/Windsor_SB_fullScale_multiPhase_uncoupled/X_P18_637/T1002_T3200_F50_D20_D400_cumulative.mat', 'mapData').mapData;

FS_Coupled_1L = load('/mnt/Processing/Data/PhD/Numerical/MATLAB/planarSprayMap/Windsor_fullScale/Windsor_SB_fullScale_multiPhase_coupled/X_P6_109/T1002_T3200_F50_D20_D400_cumulative.mat', 'mapData').mapData;
FS_Coupled_2L = load('/mnt/Processing/Data/PhD/Numerical/MATLAB/planarSprayMap/Windsor_fullScale/Windsor_SB_fullScale_multiPhase_coupled/X_P10_285/T1002_T3200_F50_D20_D400_cumulative.mat', 'mapData').mapData;
FS_Coupled_3L = load('/mnt/Processing/Data/PhD/Numerical/MATLAB/planarSprayMap/Windsor_fullScale/Windsor_SB_fullScale_multiPhase_coupled/X_P14_461/T1002_T3200_F50_D20_D400_cumulative.mat', 'mapData').mapData;
FS_Coupled_4L = load('/mnt/Processing/Data/PhD/Numerical/MATLAB/planarSprayMap/Windsor_fullScale/Windsor_SB_fullScale_multiPhase_coupled/X_P18_637/T1002_T3200_F50_D20_D400_cumulative.mat', 'mapData').mapData;


%%

% maxSpray_QS = [
%                max(full(QS_1L.areaDensity.mean));
%                max(full(QS_2L.areaDensity.mean));
%                NaN;
%                NaN
%               ];

maxSpray_FS_Uncoupled = [
                         max(full(FS_Uncoupled_1L.areaDensity.mean));
                         max(full(FS_Uncoupled_2L.areaDensity.mean));
                         max(full(FS_Uncoupled_3L.areaDensity.mean));
                         max(full(FS_Uncoupled_4L.areaDensity.mean))
                        ];

maxSpray_FS_Coupled = [
                       max(full(FS_Coupled_1L.areaDensity.mean));
                       max(full(FS_Coupled_2L.areaDensity.mean));
                       max(full(FS_Coupled_3L.areaDensity.mean));
                       max(full(FS_Coupled_4L.areaDensity.mean))
                      ];


%%

% Initialise Figure
fig = fig + 1;
figName = 'Maximum_Per_Plane_Area_Density';
set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
             'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
pause(0.5);
hold on;
set(gca, 'positionConstraint', 'outerPosition', 'plotBoxAspectRatio', [1, 0.75, 0.75], ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 22, 'layer', 'top');


% Plot Mass
% plot(maxSpray_QS, 'color', graphColours(1), 'lineStyle', '-', 'lineWidth', 2, ...
%                   'marker', 'o', 'markerSize', 10, 'markerFaceColor', graphColours(1));
plot(maxSpray_FS_Uncoupled, 'color', graphColours(4), 'lineStyle', '-', 'lineWidth', 2, ...
                            'marker', 'o', 'markerSize', 10, 'markerFaceColor', graphColours(4));
plot(maxSpray_FS_Coupled, 'color', graphColours(5), 'lineStyle', '-', 'lineWidth', 2, ...
                          'marker', 'o', 'markerSize', 10, 'markerFaceColor', graphColours(5));

% Format Figure
title('{-----}', 'interpreter', 'latex');
subtitle('{ }');
axis on;
box on;
grid off;
xlim([0.5; 4.5]);
ylim([0; 1]);
tickData = [1; 2; 3; 4];
xticks(tickData);
tickData = (200e-3:200e-3:800e-3);
yticks(tickData);
currentAxisX = get(gca, 'XAxis'); currentAxisX.TickLabelInterpreter = 'latex'; clear currentAxisX;
set(gca, 'XTickLabel', {'$1.0\,\ell$', '$2.0\,\ell$', '$3.0\,\ell$', '$4.0\,\ell$'});
% ytickformat('%.1f');
currentFig = gca; currentFig.YAxis.Exponent = -3; clear currentAxis;
xlabel({'{Measurement Plane}'; '{-----}'}, 'interpreter', 'latex');
ylabel({'{-----}'; '{$\max{(\overline{\varrho_{_{A}}})}$ ($kg {\cdot} m^{-2}$)}'}, 'interpreter', 'latex');
legend({'Uncoupled', ...
        'Coupled'}, ...
       'location', 'northEast', 'orientation', 'vertical', 'interpreter', 'latex', ...
       'fontSize', 18, 'box', 'off');
tightInset = get(gca, 'TightInset');
set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                           (tightInset(2) + 0.00625), ...
                           (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                           (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
pause(0.5);
hold off;

% Save Figure
print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');