clear variables;
close all;
clc;
evalc('delete(gcp(''nocreate''));');

fig = 0; % Initialise Figure Tracking
figHold = 0; % Enable Overwriting of Figures

figSave = false; % Save .fig File(s);

%%%

% Select Relevant Geometry and Define Bounding Box
[geometry, xDims, yDims, zDims, spacePrecision, normDims, normLength] = selectGeometry(true);


%%

% caseA = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/SB_1.0L_120s_15Hz_01/T0067_T120000_F15_Norm.mat';
% caseB = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/SB_1.0L_120s_15Hz_02/T0067_T120000_F15_Norm.mat';
% caseC = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/SB_1.0L_120s_15Hz_03/T0067_T120000_F15_Norm.mat';

% caseA = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/ST_1.0L_120s_15Hz_01/T0067_T120000_F15_Norm.mat';
% caseB = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/ST_1.0L_120s_15Hz_02/T0067_T120000_F15_Norm.mat';
% caseC = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/ST_1.0L_120s_15Hz_03/T0067_T120000_F15_Norm.mat';

caseA = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/RSST_1.0L_120s_15Hz_01/T0067_T120000_F15_Norm.mat';
caseB = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/RSST_1.0L_120s_15Hz_02/T0067_T120000_F15_Norm.mat';
caseC = '/mnt/Processing/Data/Experimental/MATLAB/planarSprayMap/RSST_1.0L_120s_15Hz_03/T0067_T120000_F15_Norm.mat';

% Load Data
R2R.mapDataA = load(caseA, 'mapData').mapData;
R2R.mapDataB = load(caseB, 'mapData').mapData;
R2R.mapDataC = load(caseC, 'mapData').mapData;

sprayMaps = fieldnames(R2R);


%%

xLimsData = R2R.mapDataA.positionGrid(1,1);
yLimsData = [min(R2R.mapDataA.positionGrid(:,2)); max(R2R.mapDataA.positionGrid(:,2))];
zLimsData = [min(R2R.mapDataA.positionGrid(:,3)); max(R2R.mapDataA.positionGrid(:,3))];

xLimsPlot = [0.31875; 2.73575];
yLimsPlot = [-0.522; 0.522];
zLimsPlot = [0; 0.522];

xLimsPlot = round((xLimsPlot / normLength), spacePrecision);
yLimsPlot = round((yLimsPlot / normLength), spacePrecision);
zLimsPlot = round((zLimsPlot / normLength), spacePrecision);

spatialRes = 0.5e-3;
    
positionData = R2R.mapDataA.positionGrid;

figTitle = '{ }'; % Leave Blank ('{ }') for Formatting Purposes

figName = 'Spray_Boundary_Variation_RSST_1.0L';

contourlines = [0.02; 0.02];


%%

fig = fig + 1;
set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
                 'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
pause(0.5);
hold on;
set(gca, 'positionConstraint', 'outerPosition', 'dataAspectRatio', [1, 1, 1], ...
         'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 20, 'layer', 'top');
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
for i = 1:height(sprayMaps)
    sprayMap = sprayMaps{i};
    
    scalarData = R2R.(sprayMap).density.mean;
    
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
    set(contours, 'edgeColor', graphColours(i), 'lineStyle', '-', 'lineWidth', 2);
end

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
tickData = (yLimsPlot(1):((yLimsPlot(2) - yLimsPlot(1)) / 5):yLimsPlot(2));
yticks(tickData(2:(end-1)));
tickData = (zLimsPlot(1):((zLimsPlot(2) - zLimsPlot(1)) / 5):zLimsPlot(2));
zticks(tickData(2:(end-1)));
xtickformat('%+.2f');
ytickformat('%+.2f');
ztickformat('%+.2f');

if normDims
    ylabel({'{$y_{\ell}$}'; '{-----}'}, 'interpreter', 'latex');
    zlabel({'{-----}'; '{$z_{\ell}$}'}, 'interpreter', 'latex');
else
    ylabel({'{$y$ ($m$)}'; '{-----}'}, 'interpreter', 'latex');
    zlabel({'{-----}'; '{$z$ ($m$)}'}, 'interpreter', 'latex');
end

tightInset = get(gca, 'TightInset');
set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                           (tightInset(2) + 0.00625), ...
                           (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                           (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
pause(0.5);
hold off;

% Save Figure
print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');