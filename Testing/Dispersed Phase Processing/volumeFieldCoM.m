run preamble;

volumeData.uncoupled = load('/mnt/Processing/Data/Numerical/MATLAB/volumeField/Windsor_fullScale/Windsor_SB_fullScale_multiPhase_uncoupled/farWake/T1002_T3200_F50_D20_D400.mat', 'volumeData').volumeData;
volumeData.coupled = load('/mnt/Processing/Data/Numerical/MATLAB/volumeField/Windsor_fullScale/Windsor_SB_fullScale_multiPhase_coupled/farWake/T1002_T3200_F50_D20_D400.mat', 'volumeData').volumeData;
volumeData.halfTread = load('/mnt/Processing/Data/Numerical/MATLAB/volumeField/Windsor_fullScale/Windsor_SB_fullScale_multiPhase_halfTread/farWake/T1002_T3200_F50_D20_D400.mat', 'volumeData').volumeData;
volumeData.deg20 = load('/mnt/Processing/Data/Numerical/MATLAB/volumeField/Windsor_fullScale/Windsor_SB_fullScale_multiPhase_20deg/farWake/T1002_T3200_F50_D20_D400.mat', 'volumeData').volumeData;

cases = fieldnames(volumeData);

[geometry, xDims, yDims, zDims, spacePrecision, normLength] = selectGeometry(geoLoc);


%%

for i = 1:height(cases)
    xPos = sort(unique(volumeData.(cases{i}).positionGrid(:,1))); xPos = xPos(xPos > xDims(2));

    CoM.(cases{i}) = zeros([height(xPos),3]);

    for j = 1:height(xPos)
        index = find(volumeData.(cases{i}).positionGrid(:,1) == xPos(j));

        CoM.(cases{i})(j,1) = xPos(j);
        CoM.(cases{i})(j,2) = full(sum(volumeData.(cases{i}).density.mean(index) .* volumeData.(cases{i}).positionGrid(index,2)) / ...
                                   sum(volumeData.(cases{i}).density.mean(index)));
        CoM.(cases{i})(j,3) = full(sum(volumeData.(cases{i}).density.mean(index) .* volumeData.(cases{i}).positionGrid(index,3)) / ...
                                   sum(volumeData.(cases{i}).density.mean(index)));
    end
    clear j;
    
end
clear i;


%%

parts = fieldnames(geometry);
for i = 1:height(parts)
    geometry.(parts{i}).vertices = geometry.(parts{i}).vertices / normLength;
end
clear i parts;

xDims = xDims / normLength;
yDims = yDims / normLength;
zDims = zDims / normLength;

for i = 1:height(cases)
    CoM.(cases{i}) = CoM.(cases{i}) / normLength;
end
clear i;
     

%%

xLimsPlot = [0.3; 1.2];
yLimsPlot = [-0.3; 0.3];
zLimsPlot = [0; 0.5];

% xLimsPlot = [-0.637116858237548; 4.562883141762452];
% yLimsPlot = [-0.7; 0.7];
% zLimsPlot = [0; 0.7];

% Initialise Figure
fig = fig + 1;
figName = 'CoM_Evolution';
set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
                 'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
pause(0.5);
hold on;
set(gca, 'positionConstraint', 'outerPosition', 'dataAspectRatio', [1, 1, 1], ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 22, 'layer', 'top');
lighting gouraud;

% Plot Geometry
parts = fieldnames(geometry);
for i = 1:height(parts)
    patch('faces', geometry.(parts{i}).faces, ...
          'vertices', geometry.(parts{i}).vertices, ...
          'faceColor', [0.5, 0.5, 0.5], ...
          'edgeColor', [0.5, 0.5, 0.5], ...
          'lineStyle', 'none');
end
clear i parts;

% Plot CoM
for i = 1:height(cases)
    plot3(CoM.(cases{i})(:,1), CoM.(cases{i})(:,2), CoM.(cases{i})(:,3), 'color', graphColours(3+i), 'lineWidth', 2);
end
clear i;

% Format Figure
title('{-----}', 'interpreter', 'latex');
subtitle('{ }');
lightangle(0, 45);
axis on;
box on;
grid off
view([0,0]);
xlim([xLimsPlot(1), xLimsPlot(2)]);
ylim([yLimsPlot(1), yLimsPlot(2)]);
zlim([zLimsPlot(1), zLimsPlot(2)]);
tickData = xLimsPlot(1):(diff(xLimsPlot) / 5):xLimsPlot(2);
xticks(tickData(2:5));
tickData = [];
yticks(tickData);
tickData = zLimsPlot(1):(diff(zLimsPlot) / 5):zLimsPlot(2);
zticks(tickData(2:5));
xtickformat('%+.2g');
ytickformat('%+.2g');
ztickformat('%+.2g');
xlabel({'{$x_{_{\ell}}$}'; '{-----}'}, 'interpreter', 'latex');
zlabel({'{-----}'; '{$z_{_{\ell}}$}'}, 'interpreter', 'latex');
tightInset = get(gca, 'TightInset');
set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                           (tightInset(2) + 0.00625), ...
                           (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                           (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
legProps = legend({'', '', '', '', '', ...
                   'Uncoupled', ...
                   'Coupled', ...
                   'Reduced Mass', ...
                   'Increased Angle'}, ...
                   'location', 'northEast', 'orientation', 'vertical', 'interpreter', 'latex', ...
                   'fontSize', 18, 'box', 'off');
legProps.Position(2) = legProps.Position(2) - 0.26;

pause(0.5);
hold off;

% Save Figure
print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');