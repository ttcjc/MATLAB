%% Preamble

run preamble;

cellSize.target = 8e-3;
% cellSize.target = 32e-3;

figSave = false;


%%

content = importdata('/mnt/Processing/Data/Numerical/ParaView/Windsor_Upstream_2023/Windsor_SB_wW_Upstream_SC/volumeData/localFlowField.csv');
% content = importdata('/mnt/Processing/Data/Numerical/ParaView/Windsor_fullScale/Windsor_SB_fullScale_multiPhase_uncoupled/volumeData/localFlowField.csv');

rawData.positionGrid = content.data(:,[1,2,3]); rawData.positionGrid(:,1) = rawData.positionGrid(:,1) + 1.325;
rawData.u.mean = content.data(:,4);
rawData.v.mean = content.data(:,5);
rawData.w.mean = content.data(:,6);

clear content;

[~, index] = unique(rawData.positionGrid, 'rows');

rawData.positionGrid = rawData.positionGrid(index,:);
rawData.u.mean = rawData.u.mean(index);
rawData.v.mean = rawData.v.mean(index);
rawData.w.mean = rawData.w.mean(index);

[rawData.positionGrid, index] = sortrows(rawData.positionGrid, [3, 2, 1]);
rawData.u.mean = rawData.u.mean(index);
rawData.v.mean = rawData.v.mean(index);
rawData.w.mean = rawData.w.mean(index);


%% Select Relevant Geometry and Define Bounding Box

[geometry, xDims, yDims, zDims, spacePrecision, normLength] = selectGeometry(geoLoc);


%% Interpolate Volume Data Onto Uniform Grid

xLimsData = [min(rawData.positionGrid(:,1)); max(rawData.positionGrid(:,1))];
yLimsData = [min(rawData.positionGrid(:,2)); max(rawData.positionGrid(:,2))];
zLimsData = [min(rawData.positionGrid(:,3)); max(rawData.positionGrid(:,3))];

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

nCells = height(volumeData.positionGrid);

uInterp = scatteredInterpolant(rawData.positionGrid(:,1), ...
                               rawData.positionGrid(:,2), ...
                               rawData.positionGrid(:,3), ...
                               rawData.u.mean, ...
                               'linear', 'none');

vInterp = scatteredInterpolant(rawData.positionGrid(:,1), ...
                               rawData.positionGrid(:,2), ...
                               rawData.positionGrid(:,3), ...
                               rawData.v.mean, ...
                               'linear', 'none');

wInterp = scatteredInterpolant(rawData.positionGrid(:,1), ...
                               rawData.positionGrid(:,2), ...
                               rawData.positionGrid(:,3), ...
                               rawData.w.mean, ...
                               'linear', 'none');

volumeData.u.mean = uInterp(volumeData.positionGrid(:,1), ...
                            volumeData.positionGrid(:,2), ...
                            volumeData.positionGrid(:,3));
volumeData.u.mean = single(volumeData.u.mean);

volumeData.v.mean = vInterp(volumeData.positionGrid(:,1), ...
                            volumeData.positionGrid(:,2), ...
                            volumeData.positionGrid(:,3));
volumeData.v.mean = single(volumeData.v.mean);

volumeData.w.mean = wInterp(volumeData.positionGrid(:,1), ...
                            volumeData.positionGrid(:,2), ...
                            volumeData.positionGrid(:,3));
volumeData.uwmean = single(volumeData.w.mean);


%% Calculate Time-Averaged Lambda2

gridShape = [height(unique(volumeData.positionGrid(:,1))), ...
             height(unique(volumeData.positionGrid(:,2))), ...
             height(unique(volumeData.positionGrid(:,3)))];

volumeData.lambda2.mean = zeros([nCells,1], 'single');

x = unique(volumeData.positionGrid(:,1));
y = unique(volumeData.positionGrid(:,2));
z = unique(volumeData.positionGrid(:,3));

% Calculate grad(U)
u = reshape(volumeData.u.mean, gridShape); u = permute(u, [2, 1, 3]);
v = reshape(volumeData.v.mean, gridShape); v = permute(v, [2, 1, 3]);
w = reshape(volumeData.w.mean, gridShape); w = permute(w, [2, 1, 3]);

[dudx, dudy, dudz] = gradient(u, x, y, z);
[dvdx, dvdy, dvdz] = gradient(v, x, y, z);
[dwdx, dwdy, dwdz] = gradient(w, x, y, z);

dudx = permute(dudx, [2, 1, 3]); dudx = dudx(:);
dudy = permute(dudy, [2, 1, 3]); dudy = dudy(:);
dudz = permute(dudz, [2, 1, 3]); dudz = dudz(:);

dvdx = permute(dvdx, [2, 1, 3]); dvdx = dvdx(:);
dvdy = permute(dvdy, [2, 1, 3]); dvdy = dvdy(:);
dvdz = permute(dvdz, [2, 1, 3]); dvdz = dvdz(:);

dwdx = permute(dwdx, [2, 1, 3]); dwdx = dwdx(:);
dwdy = permute(dwdy, [2, 1, 3]); dwdy = dwdy(:);
dwdz = permute(dwdz, [2, 1, 3]); dwdz = dwdz(:);

% Calculate Lambda2    
for i = 1:nCells
    J = [dudx(i), dudy(i), dudz(i); dvdx(i), dvdy(i), dvdz(i); dwdx(i), dwdy(i), dwdz(i)];
    J(isnan(J)) = 0;

    S = 0.5 * (J + J');
    W = 0.5 * (J - J');

    lambda = eig(S.^2 + W.^2); lambda = sort(lambda);

    volumeData.lambda2.mean(i) = lambda(2);
end
clear i;


%% Select Presentation Options

valid = false;
while ~valid
    disp(' ');
    
    selection = input('Plot Time-Averaged Lambda2? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        plotMean = false;
        
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        plotMean = true;
        
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end

if plotMean
    
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
end


%% Present Volume Fields

if plotMean
    spatialRes = cellSize.target / 2;
    xOrig = reshape(volumeData.positionGrid(:,1), gridShape);
    yOrig = reshape(volumeData.positionGrid(:,2), gridShape);
    zOrig = reshape(volumeData.positionGrid(:,3), gridShape);
    POD = false;
    fieldData = reshape((full(volumeData.lambda2.mean) / -prctile(volumeData.lambda2.mean(volumeData.lambda2.mean < 0), 1)), gridShape);
    isoValue = -0.01;
    nSurfaces = 1;
    surfaceNo = 1;
    cMap = graphColours(1);
    figTitle = '{ }'; % Leave Blank ('{ }') for Formatting Purposes
    viewAngle = [30, 30];
    multiView = true;
    xLimsPlot = [-0.637116858237548; 1.562883141762452];
    yLimsPlot = [-0.5; 0.5];
    zLimsPlot = [0; 0.5];
    
    for i = 1:height(isoValue)
        figName = ['Local_Lambda2_N', num2str(abs(isoValue(i)))];
        
        [fig, surfaceNo] = plotVolumeField(xLimsData, yLimsData, zLimsData, spatialRes, ...
                                           xOrig, yOrig, zOrig, POD, fieldData, nSurfaces, surfaceNo, ...
                                           fig, figName, geometry, isoValue(i), cMap, figTitle, viewAngle, ...
                                           multiView, xLimsPlot, yLimsPlot, zLimsPlot, figSave);
    end
    clear i;
    
    disp(' ');
end