%% Volumetric Velocity Field Processing v1.0
% ----
% Load, Process and Present Volumetric Flow Field Data Acquired Using OpenFOAM v7


%% Preamble

run preamble;

%#ok<*UNRCH>

normDims = true; % Normalise Spatial Dimensions

figSave = false; % Save .fig File(s)

disp('=========================================');
disp('Volumetric Velocity Field Processing v1.0');
disp('=========================================');

disp(' ');
disp(' ');


%% Changelog

% v1.0 - Origial Commit


%% Initialise Case

maxNumCompThreads(nProc);

[caseFolder, campaignID, caseID, timeDirs, deltaT, ...
 timePrecision, dataID, probeData, sampleInt] = initialiseProbeData(saveLoc, maxNumCompThreads, ...
                                                                    'uProbes', 'wake');

disp(' ');
disp(' ');


%% Select Relevant Geometry and Define Bounding Box

[geometry, xDims, yDims, zDims, spacePrecision, normLength] = selectGeometry(geoLoc);

disp(' ');
disp(' ');

%% Select Region of Interest

disp('Region of Interest');
disp('-------------------');

disp(' ');

disp('Possible Regions of Interest:');
disp('    A: Near Wake');
disp('    B: Mid Wake');
disp('    C: Far Wake (Full-Scale Only)');

valid = false;
while ~valid
    disp(' ');
    
    selection = input('Select Region of Interest [A/B/C]: ', 's');

    if selection == 'a' | selection == 'A' %#ok<OR2>
        format = 'A';
        
        valid = true;
    elseif selection == 'b' | selection == 'B' %#ok<OR2>
        format = 'B';
        
        valid = true;
    elseif selection == 'c' | selection == 'C' %#ok<OR2>
        format = 'C';
        
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;

if strcmp(format, 'C') && ~strcmp(campaignID, 'Windsor_fullScale')
    disp(' ');
    
    disp('WARNING: Far-Wake Data Is Unavailable for This Case');
    disp('         Performing Analysis on Mid-Wake Data Instead');
    
    format = 'B';
end

disp(' ');
disp(' ');


%% Generate Volume Field

disp('Volume Field Generation');
disp('------------------------');

disp(' ');

disp('***********');
disp('  RUNNING ');

tic;

evalc('parpool(''threads'');');

%%%%

disp(' ');

disp('    Initialising...');

nTimes = height(probeData.time);

% Identify Empty Time Instances Not Related to Case Initialisation
emptyTimes = cellfun(@isempty, probeData.u.inst);
firstValidTime = find(emptyTimes == false, 1, 'first');

suspectTimes = find(emptyTimes(firstValidTime:end) == true) + (firstValidTime - 1);

if ~isempty(suspectTimes)
    
    for i = 1:height(suspectTimes)
        disp(['        WARNING: Time ''', num2str(probeData.time(suspectTimes(i))), ''' Is Unexpectedly Empty']);
    end
    clear i;
    
end

% Adjust Data Origin
if strcmp(campaignID, 'Windsor_Upstream_2023')
    probeData.positionGrid(:,1) = probeData.positionGrid(:,1) + 1.325;
end

% Remove Duplicate Entries
[probeData.positionGrid, index] = unique(probeData.positionGrid, 'rows', 'stable');

for i = 1:nTimes
    probeData.u.inst{i} = probeData.u.inst{i}(index);
    probeData.v.inst{i} = probeData.v.inst{i}(index);
    probeData.w.inst{i} = probeData.w.inst{i}(index);
end
clear i;

% Specify Region Boundaries
switch format
    
    case 'A' % 1 L
        xLimsData = [0.462931034482759; 1.462883141762452] * normLength;
        yLimsData = [-0.4; 0.4] * normLength;
        zLimsData = [0; 0.4] * normLength;
        
    case 'B' % 2 L
        xLimsData = [0.462931034482759; 2.462883141762452] * normLength;
        yLimsData = [-0.5; 0.5] * normLength;
        zLimsData = [0; 0.5] * normLength;
        
    case 'C' % 4 L
        xLimsData = [0.462931034482759; 4.462883141762452] * normLength;
        yLimsData = [-0.6; 0.6] * normLength;
        zLimsData = [0; 0.6] * normLength;

end

if strcmp(campaignID, 'Windsor_fullScale')
    zLimsData(1) = 0.015325670498084 * normLength;
elseif strcmp(campaignID, 'Windsor_Upstream_2023')
    zLimsData(1) = 0.007662835249042 * normLength;
end

% Remove Excess Data
index = ((probeData.positionGrid(:,1) > xLimsData(1) & probeData.positionGrid(:,1) < xLimsData(2)) & ...
         (probeData.positionGrid(:,2) > yLimsData(1) & probeData.positionGrid(:,2) < yLimsData(2)) & ...
         (probeData.positionGrid(:,3) > zLimsData(1) & probeData.positionGrid(:,3) < zLimsData(2)));

probeData.positionGrid = probeData.positionGrid(index,:);

for i = 1:nTimes
    probeData.u.inst{i} = probeData.u.inst{i}(index);
    probeData.v.inst{i} = probeData.v.inst{i}(index);
    probeData.w.inst{i} = probeData.w.inst{i}(index);
end
clear i;

% Sort Position Grid for 'ndgrid' Compatibility
[probeData.positionGrid, index] = sortrows(probeData.positionGrid, [3, 2, 1]);

for i = 1:nTimes
    probeData.u.inst{i} = probeData.u.inst{i}(index);
    probeData.v.inst{i} = probeData.v.inst{i}(index);
    probeData.w.inst{i} = probeData.w.inst{i}(index);
end
clear i;

disp(' ');

% Generate Instantaneous Volume Field
disp('    Generating Presentation Grid...');

% Set Target Spatial Resolution
if strcmp(campaignID, 'Windsor_fullScale')
    cellSize.target = 64e-3;
elseif strcmp(campaignID, 'Windsor_Upstream_2023')
    cellSize.target = 8e-3;
else
    cellSize.target = 8e-3;
end

% Adjust Uniform Cell Size to Fit Region of Interest
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

% Generate Grid
[x, y, z] = ndgrid(linspace(xLimsData(1), xLimsData(2), nPx), ...
                   linspace(yLimsData(1), yLimsData(2), nPy), ...
                   linspace(zLimsData(1), zLimsData(2), nPz));

uData.positionGrid = [x(:), y(:), z(:)];

nCells = height(uData.positionGrid);

disp(' ');

% Interpolate Volume Data Onto Uniform Grid
disp('    Interpolating Data Onto Uniform Grid...');

uData.time = probeData.time; probeData.time = -1;

gridShape = [height(unique(probeData.positionGrid(:,1))), ...
             height(unique(probeData.positionGrid(:,2))), ...
             height(unique(probeData.positionGrid(:,3)))];

% Initialise Progress Bar
wB = waitbar(0, 'Interpolating Data Onto Uniform Grid', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';
dQ = parallel.pool.DataQueue;
afterEach(dQ, @parforWaitBar);
parforWaitBar(wB, nTimes);

% Perform Interpolation
u = cell(nTimes,1); u(:) = {zeros([nCells, 1], 'single')};
v = u;
w = u;

xOrig = double(reshape(probeData.positionGrid(:,1), gridShape));
yOrig = double(reshape(probeData.positionGrid(:,2), gridShape));
zOrig = double(reshape(probeData.positionGrid(:,3), gridShape));
uOrig = probeData.u.inst; probeData.u.inst = -1;
vOrig = probeData.v.inst; probeData.v.inst = -1;
wOrig = probeData.w.inst; probeData.w.inst = -1;
parfor i = 1:nTimes    
    uInterp = griddedInterpolant(xOrig, yOrig, zOrig, double(reshape(uOrig{i}, gridShape)), 'linear', 'none');
    vInterp = griddedInterpolant(xOrig, yOrig, zOrig, double(reshape(vOrig{i}, gridShape)), 'linear', 'none');
    wInterp = griddedInterpolant(xOrig, yOrig, zOrig, double(reshape(wOrig{i}, gridShape)), 'linear', 'none');
    
    u{i} = uInterp(x, y, z); u{i} = single(u{i}(:));
    v{i} = vInterp(x, y, z); v{i} = single(v{i}(:));
    w{i} = wInterp(x, y, z); w{i} = single(w{i}(:));
    
    % Remove Unnecessary Data
    uOrig{i} = -1;
    vOrig{i} = -1;
    wOrig{i} = -1;
    
    % Update Waitbar
    send(dQ, []);
end
clear xOrig yOrig zOrig uOrig vOrig wOrig;

delete(wB);

clear probeData x y z;

uData.u.inst = u; clear u;
uData.v.inst = v; clear v;
uData.w.inst = w; clear w;

disp(' ');

% Calculate Instantaneous Lambda2
disp('    Calculating Instantaneous Lambda2...');

gridShape = [height(unique(uData.positionGrid(:,1))), ...
             height(unique(uData.positionGrid(:,2))), ...
             height(unique(uData.positionGrid(:,3)))];

lambda2 = cell(nTimes,1); lambda2(:) = {zeros([nCells, 1], 'single')};

% Initialise Progress Bar
wB = waitbar(0, 'Calculating Instantaneous Lambda2 Field', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';
dQ = parallel.pool.DataQueue;
afterEach(dQ, @parforWaitBar);
parforWaitBar(wB, nTimes);

% Perform Calculation
uInst = uData.u.inst;
vInst = uData.v.inst;
wInst = uData.w.inst;
x = unique(uData.positionGrid(:,1));
y = unique(uData.positionGrid(:,2));
z = unique(uData.positionGrid(:,3));
parfor i = 1:nTimes

    % Calculate grad(U)
    u = reshape(uInst{i}, gridShape); u = permute(u, [2, 1, 3]);
    v = reshape(vInst{i}, gridShape); v = permute(v, [2, 1, 3]);
    w = reshape(wInst{i}, gridShape); w = permute(w, [2, 1, 3]);

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
    for j = 1:nCells
        J = [dudx(j), dudy(j), dudz(j); dvdx(j), dvdy(j), dvdz(j); dwdx(j), dwdy(j), dwdz(j)];
        J(isnan(J)) = 0;

        S = 0.5 * (J + J');
        W = 0.5 * (J - J');

        lambda = eig(S.^2 + W.^2); lambda = sort(lambda);

        lambda2{i}(j) = lambda(2);
    end

    % Update Waitbar
    send(dQ, []);
end
clear uInst vInst wInst x y z;

delete(wB);

uData.lambda2.inst = lambda2; clear lambda2;

disp(' ');

% Calculate Instantaneous Field Variables
disp('    Calculating Time-Averaged Velocity Field...');

% Initialise Progress Bar
wB = waitbar(0, 'Calculating Time-Averaged Velocity Field', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';

% Perform Calculation
uData.u.mean = zeros([nCells,1], 'single');
uData.v.mean = uData.u.mean;
uData.w.mean = uData.u.mean;

for i = 1:nTimes
    uData.u.mean = uData.u.mean + uData.u.inst{i};
    uData.v.mean = uData.v.mean + uData.v.inst{i};
    uData.w.mean = uData.w.mean + uData.w.inst{i};

    % Update Waitbar
    waitbar((i / nTimes), wB);
end
clear i;

delete(wB);

uData.u.mean = uData.u.mean / nTimes;
uData.v.mean = uData.v.mean / nTimes;
uData.w.mean = uData.w.mean / nTimes;

disp(' ');

% Calculate Time-Averaged Lambda2
disp('    Calculating Time-Averaged Lambda2...');

uData.lambda2.mean = zeros([nCells,1], 'single');

x = unique(uData.positionGrid(:,1));
y = unique(uData.positionGrid(:,2));
z = unique(uData.positionGrid(:,3));

% Calculate grad(U)
u = reshape(uData.u.mean, gridShape); u = permute(u, [2, 1, 3]);
v = reshape(uData.v.mean, gridShape); v = permute(v, [2, 1, 3]);
w = reshape(uData.w.mean, gridShape); w = permute(w, [2, 1, 3]);

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

    uData.lambda2.mean(i) = lambda(2);
end
clear i;

%%%%

evalc('delete(gcp(''nocreate''));');

executionTime = toc;

disp(' ');

disp(['    Run Time: ', num2str(executionTime), 's']);

disp(' ');

disp('  SUCCESS  ');
disp('***********');

disp(' ');
disp(' ');


%% Select Presentation Options

disp('Presentation Options');
disp('---------------------');

valid = false;
while ~valid
    disp(' ');
    
    selection = input('Plot Time-Averaged Velocity Deficit? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        plotMeanU = false;
        
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        plotMeanU = true;
        
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end

valid = false;
while ~valid
    disp(' ');
    
    selection = input('Plot Time-Averaged Lambda2? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        plotMeanL2 = false;
        
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        plotMeanL2 = true;
        
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
    
valid = false;
while ~valid
    disp(' ');

    selection = input('Plot Instantaneous Lambda2? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        plotInstL2 = false;

        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        plotInstL2 = true;

        startFrame = inputFrames(nTimes, 'Start');

        if startFrame == -1
            continue;
        end

        endFrame = inputFrames(nTimes, 'End');

        if endFrame == -1
            continue;
        elseif endFrame < startFrame
            disp('        WARNING: Invalid Time Format (''endFrame'' Precedes ''startFrame'')');

            continue;
        end

        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;

if plotMeanU || plotMeanL2 || plotInstL2
    
    % Normalise Coordinate System
    if normDims
        disp(' ');

        disp('Normalising Spatial Dimensions...');

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

        uData.positionGrid = uData.positionGrid / normLength;
    end
    
    disp(' ');
    
    % Normalise Velocity
    disp('Normalising Velocity...');
    
    if strcmp(campaignID, 'Windsor_fullScale')
        U = 22.2222; % m/s
    elseif strcmp(campaignID, 'Windsor_Upstream_2023')
        U = 40; % m/s
    else
        U = max(uData.u.mean);
    end
    
    uData.u.mean = uData.u.mean / U;
    uData.v.mean = uData.v.mean / U;
    uData.w.mean = uData.w.mean / U;
    
    for i = 1:nTimes
        uData.u.inst{i} = uData.u.inst{i} / U;
        uData.v.inst{i} = uData.v.inst{i} / U;
        uData.w.inst{i} = uData.w.inst{i} / U;
    end
    clear i;
    
end

disp(' ');
disp(' ');


%% Present Volume Fields

disp('Volume Field Presentation');
disp('--------------------------');

disp(' ');

if plotMeanU || plotMeanL2 || plotInstL2
             
    spatialRes = cellSize.target / 2;
    xOrig = reshape(uData.positionGrid(:,1), gridShape);
    yOrig = reshape(uData.positionGrid(:,2), gridShape);
    zOrig = reshape(uData.positionGrid(:,3), gridShape);
    POD = false;
    nSurfaces = 1;
    surfaceNo = 1;
    
    if strcmp(campaignID, 'Windsor_fullScale')
        
        if strcmp(caseID, 'Windsor_SB_fullScale_multiPhase_uncoupled')
            cMap = graphColours(4);
        elseif strcmp(caseID, 'Windsor_SB_fullScale_multiPhase_coupled')
            cMap = graphColours(5);
        end
        
    elseif strcmp(campaignID, 'Windsor_Upstream_2023')
        
        if strcmp(caseID, 'Windsor_SB_wW_Upstream_SC')
            cMap = graphColours(1);
        elseif strcmp(caseID, 'Windsor_ST_wW_Upstream_SC')
            cMap = graphColours(2);
        elseif strcmp(caseID, 'Windsor_RSST_wW_Upstream_SC')
            cMap = graphColours(3);
        end
        
    end
    
    if ~exist('cMap', 'var')
        cMap = graphColours(1);
    end
    
    viewAngle = [30, 30];
    
    switch format

        case 'A' % 1 L
            xLimsPlot = [0.3; 1.562883141762452];
            yLimsPlot = [-0.5; 0.5];
            zLimsPlot = [0; 0.5];

        case 'B' % 2 L
            xLimsPlot = [-0.637116858237548; 2.562883141762452];
            yLimsPlot = [-0.6; 0.6];
            zLimsPlot = [0; 0.6];

        case 'C' % 4 L
            xLimsPlot = [-0.637116858237548; 4.562883141762452];
            yLimsPlot = [-0.7; 0.7];
            zLimsPlot = [0; 0.7];

    end
    
    if ~normDims
        xLimsPlot = xLimsPlot * normLength;
        yLimsPlot = yLimsPlot * normLength;
        zLimsPlot = zLimsPlot * normLength;
    end
    
end

if plotMeanU
    disp('Presenting Time-Averaged Velocity Deficit...');
    
    % Calculate Velocity Magnitude
    fieldData = reshape(sqrt(uData.u.mean.^2 + uData.v.mean.^2 + uData.w.mean.^2), gridShape);
    
    % Calculate Velocity Deficit
    fieldData = fieldData - 1;
    isoValue = -0.05;
    figTitle = '{ }'; % Leave Blank ('{ }') for Formatting Purposes
    multiView = true;
             
    for i = 1:height(isoValue)
        
        switch format

            case 'A'
                figName = ['NW_Average_Velocity_Deficit_', num2str(isoValue(i)), '_', caseID];

            case 'B'
                figName = ['MW_Average_Velocity_Deficit_', num2str(isoValue(i)), '_', caseID];

            case 'C'
                figName = ['FW_Average_Velocity_Deficit_', num2str(isoValue(i)), '_', caseID];

        end
        
        [fig, surfaceNo] = plotVolumeField(xLimsData, yLimsData, zLimsData, spatialRes, ...
                                           xOrig, yOrig, zOrig, POD, fieldData, nSurfaces, surfaceNo, ...
                                           fig, figName, geometry, isoValue(i), cMap, figTitle, viewAngle, ...
                                           multiView, xLimsPlot, yLimsPlot, zLimsPlot, figSave);
    end
    clear i;
                       
    disp(' ');
end

if plotMeanL2
    disp('Presenting Time-Averaged Lambda2...');
    
    fieldData = reshape((full(uData.lambda2.mean) / -prctile(uData.lambda2.mean(uData.lambda2.mean < 0), 1)), gridShape);
    isoValue = -0.2;
    figTitle = '{ }'; % Leave Blank ('{ }') for Formatting Purposes
    multiView = true;
             
    for i = 1:height(isoValue)
        
        switch format

            case 'A'
                figName = ['NW_Average_Lambda2_', num2str(isoValue(i)), '_', caseID];

            case 'B'
                figName = ['MW_Average_Lambda2_', num2str(isoValue(i)), '_', caseID];

            case 'C'
                figName = ['FW_Average_Lambda2_', num2str(isoValue(i)), '_', caseID];

        end
        
        [fig, surfaceNo] = plotVolumeField(xLimsData, yLimsData, zLimsData, spatialRes, ...
                                           xOrig, yOrig, zOrig, POD, fieldData, nSurfaces, surfaceNo, ...
                                           fig, figName, geometry, isoValue(i), cMap, figTitle, viewAngle, ...
                                           multiView, xLimsPlot, yLimsPlot, zLimsPlot, figSave);
    end
    clear i;
                       
    disp(' ');
end

if plotInstL2
    disp('Presenting Instantaneous Lambda2...');
    
    isoValue = -5;
    multiView = false;
    
    for i = 1:height(isoValue)
        figHold = fig;
    
        for j = startFrame:endFrame

            if j ~= startFrame
                clf(fig);
                fig = figHold;
            end
            
            fieldData = reshape((full(uData.lambda2.inst{j}) / -prctile(uData.lambda2.mean(uData.lambda2.mean < 0), 1)), gridShape);
            figTime = num2str(uData.time(j), ['%.', num2str(timePrecision), 'f']);
            
            switch format

                case 'A'
                    figName = ['NW_Inst_Lambda2_', num2str(isoValue(i)), '_T'...
                               erase(figTime, '.'), '_', caseID];

                case 'B'
                    figName = ['MW_Inst_Lambda2_', num2str(isoValue(i)), '_T'...
                               erase(figTime, '.'), '_', caseID];

                case 'C'
                    figName = ['FW_Inst_Lambda2_', num2str(isoValue(i)), '_T'...
                               erase(figTime, '.'), '_', caseID];

            end
            
            figTitle = ['{', figTime, ' \it{s}}'];
            
            [fig, surfaceNo] = plotVolumeField(xLimsData, yLimsData, zLimsData, spatialRes, ...
                                               xOrig, yOrig, zOrig, POD, fieldData, nSurfaces, surfaceNo, ...
                                               fig, figName, geometry, isoValue(i), cMap, figTitle, viewAngle, ...
                                               multiView, xLimsPlot, yLimsPlot, zLimsPlot, figSave);

        end
        clear j;
        
    end
    clear i;
    
    disp(' ');
end

if ~plotMeanU && ~plotMeanL2 && ~plotInstL2
    disp('Skipping Volume Field Presentation...');

    disp(' ');
end


%% Local Functions

function frameNo = inputFrames(Nt, type)

    frameNo = str2double(input(['    Input Desired ', type, ' Frame [1-', num2str(Nt), ']: '], 's'));
    
    if isnan(frameNo) || frameNo < 1 || frameNo > Nt
        disp('        WARNING: Invalid Entry');
        
        frameNo = -1;
    end

end
