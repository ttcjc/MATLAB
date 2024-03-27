% run preamble;

refValue = 8.996259860381801e-05;

figSave = false; % Save .fig File(s);

%%%

% Select Relevant Geometry and Define Bounding Box
[geometry, xDims, yDims, zDims, spacePrecision, normLength] = selectGeometry(geoLoc);

parts = fieldnames(geometry);
for i = 1:height(parts)
    geometry.(parts{i}).vertices = geometry.(parts{i}).vertices / normLength;
end
clear i parts;

xDims = xDims / normLength;
yDims = yDims / normLength;
zDims = zDims / normLength;


%%

load('/mnt/Processing/Data/Numerical/MATLAB/planarSprayMap/Windsor_Upstream_2023/Windsor_SB_wW_Upstream_SC/X_P1_24625/T12525_T40000_F400_D1_D147_cumulative.mat');
% load('/mnt/Processing/Data/Numerical/MATLAB/planarSprayMap/Windsor_fullScale/Windsor_SB_fullScale_multiPhase_coupled/X_P18_637/T1002_T3200_F50_D20_D400_cumulative.mat');



%%

xLimsData = mapData.positionGrid(1,1);
yLimsData = [min(mapData.positionGrid(:,2)); max(mapData.positionGrid(:,2))];
zLimsData = [min(mapData.positionGrid(:,3)); max(mapData.positionGrid(:,3))];

xLimsPlot = [0.3; 4.6257662];
yLimsPlot = [-0.5; 0.5];
zLimsPlot = [0; 0.5];

spatialRes = 0.5e-3 / normLength;
    
positionData = mapData.positionGrid / normLength;

figTitle = '{ }'; % Leave Blank ('{ }') for Formatting Purposes

figName = 'Convergence_Probes';

contourlines = [0.02; 0.02];


%%

fig = fig + 1;
set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
                 'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
pause(0.5);
hold on;
set(gca, 'positionConstraint', 'outerPosition', 'dataAspectRatio', [1, 1, 1], ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 22, 'layer', 'top');
lighting gouraud;

% Plot Geometry
if ~isempty(geometry)
    parts = fieldnames(geometry);

    for i = 1:height(parts)
        patch('faces', geometry.(parts{i}).faces, ...
              'vertices', geometry.(parts{i}).vertices, ...
              'faceColor', [0.5, 0.5, 0.5], ...
              'edgeColor', [0.5, 0.5, 0.5], ...
              'lineStyle', 'none');
    end
    clear i;

end

% Plot Contour Lines
scalarData = full(mapData.areaDensity.mean) / refValue;

% Reshape Data for Improved Interpolation Performance
gridShape = [height(unique(positionData(:,2))), ...
             height(unique(positionData(:,3)))];

y = reshape(positionData(:,2), gridShape);
z = reshape(positionData(:,3), gridShape);

scalar = reshape(scalarData, gridShape);

% Perform Interpolation
interp = griddedInterpolant(y, z, scalar, 'linear', 'none');

% Generate 3D Surface
cellSizeX = spatialRes;
cellSizeY = (yLimsData(2) - yLimsData(1)) / round(((yLimsData(2) - yLimsData(1)) / spatialRes));
cellSizeZ = (zLimsData(2) - zLimsData(1)) / round(((zLimsData(2) - zLimsData(1)) / spatialRes));

[x, y, z] = ndgrid((xLimsData - cellSizeX):cellSizeX:(xLimsData + cellSizeX), ...
                   yLimsData(1):cellSizeY:yLimsData(2), ...
                   zLimsData(1):cellSizeZ:zLimsData(2));

% Map Data on to 3D Surface
scalar = zeros(size(x));
scalar(2,:,:) = interp(y(2,:,:), z(2,:,:));
x = permute(x, [2,1,3]);
y = permute(y, [2,1,3]);
z = permute(z, [2,1,3]);
scalar = permute(scalar, [2,1,3]);

contours = contourslice(x, y, z, scalar, xLimsData, [], [], contourlines);
set(contours, 'edgeColor', graphColours(7), 'lineStyle', '-', 'lineWidth', 2);

% Plot Sample Points for Convergence Testing
for i = 1:height(index)
    plot3(positionData(index(i),1), positionData(index(i),2), positionData(index(i),3), ...
          'lineStyle', 'none', 'marker', 'o', 'markerSize', 7.5, 'markerFaceColor', graphColours(i), ...
          'markerEdgeColor', graphColours(i));
end
clear i;

title('{-----}', 'interpreter', 'latex');
subtitle(figTitle);
lightangle(90, 45);
axis on;
box on;
grid off;
view([90, 0]);
xlim([xLimsPlot(1), xLimsPlot(2)]);
ylim([yLimsPlot(1), yLimsPlot(2)]);
zlim([zLimsPlot(1), zLimsPlot(2)]);
tickData = [];
xticks(tickData);
tickData = ([]);
yticks(tickData);
tickData = ([]);
zticks(tickData);
tightInset = get(gca, 'TightInset');
set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                           (tightInset(2) + 0.00625), ...
                           (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                           (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
pause(0.5);
hold off;

% Save Figure
print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');