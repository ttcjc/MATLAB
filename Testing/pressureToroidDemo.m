%% Preamble

run preamble;


%% Load Data

content = importdata('/mnt/Processing/Data/Numerical/ParaView/Windsor_fullScale/Windsor_SB_fullScale_multiPhase_uncoupled/volumeData/p.csv');

rawData.positionGrid = content.data(:,[1,2,3]);
rawData.p = content.data(:,4);

clear content;

[~, index] = unique(rawData.positionGrid, 'rows');

rawData.positionGrid = rawData.positionGrid(index,:);
rawData.p = rawData.p(index);


%% Select Relevant Geometry and Define Bounding Box

[geometry, xDims, yDims, zDims, spacePrecision, normLength] = selectGeometry(geoLoc);


%% Adjust Data Origin

% rawData.positionGrid(:,1) = rawData.positionGrid(:,1) + 1.325;


%% Limit Volume Boundaries

xLimsData = [xDims(2); (xDims(2) + (normLength / 2))];
yLimsData = [-(normLength / 4); (normLength / 4)];
zLimsData = [0; (normLength / 2)];

index = ((rawData.positionGrid(:,1) > xLimsData(1) & rawData.positionGrid(:,1) < xLimsData(2)) & ...
         (rawData.positionGrid(:,2) > yLimsData(1) & rawData.positionGrid(:,2) < yLimsData(2)) & ...
         (rawData.positionGrid(:,3) > zLimsData(1) & rawData.positionGrid(:,3) < zLimsData(2)));

rawData.positionGrid = rawData.positionGrid(index,:);
rawData.p = rawData.p(index);


%% Sort Position Grid for 'ndgrid' Compatibility

[rawData.positionGrid, index] = sortrows(rawData.positionGrid, [3, 2, 1]);
rawData.p = rawData.p(index);


%% Interpolate Volume Data Onto Uniform Grid

cellSize.target = 16e-3;

nPx = (diff(xLimsData) / cellSize.target) + 1;
nPy = (diff(yLimsData) / cellSize.target) + 1;
nPz = (diff(zLimsData) / cellSize.target) + 1;

sizeX = diff(linspace(xLimsData(1), xLimsData(2), nPx));
sizeY = diff(linspace(yLimsData(1), yLimsData(2), nPy));
sizeZ = diff(linspace(zLimsData(1), zLimsData(2), nPz));

cellSize.x = sizeX(1); clear sizeX;
cellSize.y = sizeY(1); clear sizeY;
cellSize.z = sizeZ(1); clear sizeZ;
cellSize.volume = cellSize.x * cellSize.y * cellSize.z;

[x, y, z] = ndgrid(linspace(xLimsData(1), xLimsData(2), nPx), ...
                   linspace(yLimsData(1), yLimsData(2), nPy), ...
                   linspace(zLimsData(1), zLimsData(2), nPz));

volumeData.positionGrid = [x(:), y(:), z(:)]; clear x y z;

pInterp = scatteredInterpolant(rawData.positionGrid(:,1), ...
                               rawData.positionGrid(:,2), ...
                               rawData.positionGrid(:,3), ...
                               rawData.p, ...
                               'linear', 'none');

volumeData.p = pInterp(volumeData.positionGrid(:,1), ...
                       volumeData.positionGrid(:,2), ...
                       volumeData.positionGrid(:,3));


%% Calculate Cp

volumeData.Cp = (volumeData.p * 1.269) / (0.5 * 1.269 * 22.22^2);


%% Normalise Coordinate System

parts = fieldnames(geometry);
for i = 1:height(parts)
    geometry.(parts{i}).vertices = geometry.(parts{i}).vertices / normLength;
end
clear i parts;

xDims = xDims / normLength;
yDims = yDims / normLength;
zDims = zDims / normLength;

xLimsData = xLimsData / normLength;
yLimsData = yLimsData / normLength;
zLimsData = zLimsData / normLength;

cellSize.target = cellSize.target / normLength;
cellSize.x = cellSize.x / normLength;
cellSize.y = cellSize.y / normLength;
cellSize.z = cellSize.z / normLength;
cellSize.volume = cellSize.volume / (normLength^3);

volumeData.positionGrid = volumeData.positionGrid / normLength;


%% Present Toroid

gridShape = [height(unique(volumeData.positionGrid(:,1))), ...
             height(unique(volumeData.positionGrid(:,2))), ...
             height(unique(volumeData.positionGrid(:,3)))];

spatialRes = cellSize.target / 2;
xOrig = reshape(volumeData.positionGrid(:,1), gridShape);
yOrig = reshape(volumeData.positionGrid(:,2), gridShape);
zOrig = reshape(volumeData.positionGrid(:,3), gridShape);
POD = false;
fieldData = reshape(volumeData.Cp, gridShape);
nSurfaces = 1;
surfaceNo = 1;
isoValue = -0.343;
cMap = graphColours(1);
figTitle = '{ }'; % Leave Blank ('{ }') for Formatting Purposes
viewAngle = [60, 5];
multiView = false;
xLimsPlot = [(xDims(1) - 0.15); (xDims(2) + 0.3)];
yLimsPlot = [(yDims(1) - 0.15); (yDims(2) + 0.15)];
zLimsPlot = [-0.15; (zDims(2) + 0.15)];
figSave = false;

for i = 1:height(isoValue)
    figName = ['FS_Cp_', num2str(isoValue(i))];

    [fig, surfaceNo] = plotVolumeField(xLimsData, yLimsData, zLimsData, spatialRes, ...
                                       xOrig, yOrig, zOrig, POD, fieldData, nSurfaces, surfaceNo, ...
                                       fig, figName, geometry, isoValue(i), cMap, figTitle, viewAngle, ...
                                       multiView, xLimsPlot, yLimsPlot, zLimsPlot, figSave);
end