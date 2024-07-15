run preamble;

load('/mnt/Processing/Data/Numerical/MATLAB/probeData/Windsor_fullScale/Windsor_SB_fullScale_multiPhase_uncoupled/uProbes/wake/T1002_T3200_F50.mat');
% load('/mnt/Processing/Data/Numerical/MATLAB/probeData/Windsor_Upstream_2023/Windsor_SB_wW_Upstream_SC/uProbes/wake/T12525_T15000_F400.mat');

[geometry, xDims, yDims, zDims, spacePrecision, normLength] = selectGeometry(geoLoc);


%%

parts = fieldnames(geometry);
for i = 1:height(parts)
    geometry.(parts{i}).vertices = geometry.(parts{i}).vertices / normLength;
end
clear i parts;

xDims = xDims / normLength;
yDims = yDims / normLength;
zDims = zDims / normLength;

probeData.positionGrid = probeData.positionGrid / normLength;


%%

xLimsPlot = [(xDims(1) - 0.25); (max(probeData.positionGrid(:,1)) + 0.25)];
yLimsPlot = [(min(probeData.positionGrid(:,2)) - 0.25); (max(probeData.positionGrid(:,2)) + 0.25)];
zLimsPlot = [0; (max(probeData.positionGrid(:,3)) + 0.25)];


%%
figName = 'Continuous_Phase_Sampling_Probes';

fig = fig + 1;
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

% Plot Mesh
scatter3(probeData.positionGrid(:,1), probeData.positionGrid(:,2), probeData.positionGrid(:,3), 1, ...
         'markerFaceColor', graphColours(2), 'markerEdgeColor', 'none');

title('{-----}', 'interpreter', 'latex');
subtitle('{ }');
lightangle(0, 45);
axis on;
box on;
grid off;
view([0, 0]);
xlim([xLimsPlot(1), xLimsPlot(2)]);
ylim([yLimsPlot(1), yLimsPlot(2)]);
zlim([zLimsPlot(1), zLimsPlot(2)]);
tickData = [];
xticks(tickData);
tickData = [];
yticks(tickData);
tickData = [];
zticks(tickData);
tightInset = get(gca, 'TightInset');
set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                           (tightInset(2) + 0.00625), ...
                           (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                           (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
pause(0.5);
hold off;

viewAngle = [0, 0; 0, 90];
for i = 1:height(viewAngle)
    view(viewAngle(i,:));
    
    print(gcf, [userpath, '/Output/Figures/', figName, '_', num2str(i), '.png'], '-dpng', '-r300');
    
    pause(0.5);
end
clear i viewAngle;
    


