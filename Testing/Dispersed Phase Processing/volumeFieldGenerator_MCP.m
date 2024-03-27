%% Lagrangian Volume Field Generator v4.2
% ----
% Load, Process and Present Volumetric Lagrangian Data Acquired Using OpenFOAM v7


%% Preamble

run preamble;

%#ok<*UNRCH>

cloudName = 'kinematicCloud'; % OpenFOAM Cloud Name

normDims = true; % Normalise Spatial Dimensions

normDensity = false; % Normalise Area Density in Plots

refValue = 1;

figSave = false; % Save .fig File(s)

disp('===========================');
disp('Volume Field Generator v4.2');
disp('===========================');

disp(' ');
disp(' ');


%% Changelog

% v1.0 - Initial Commit
% v1.1 - Updated calls to 'globalPos' to 'positionCartesian'
% v1.2 - Updated to Support Changes to 'timeDirectories.m'
% v2.0 - Rewritten to Follow Recent Lagrangian Processing Structure Changes
% v2.1 - Added Time-Averaging Functionality
% v3.0 - Rewrite, Accommodating New OpenFOAM Data Formats
% v3.1 - Added Support for Arithmetic and Sauter Mean Particle Diameters
% v3.2 - Changed to Thread-Based Parallelization to Reduce Memory Requirements
% v3.3 - Added Support for Full-Scale Windsor Model Simulations
% v3.4 - Changed Primary Output From Mass to Density
% v4.0 - Rewrite, Making Use of Sparse Arrays to Reduce Memory Requirements
% v4.1 - Minor Update to Shift Preamble Into Separate Script
% v4.2 - Update To Correct Inconsistent Normalisation Throughout Repository


%% Initialise Case

[caseFolder, campaignID, caseID, timeDirs, deltaT, timePrecision, geometry, ...
 xDims, yDims, zDims, spacePrecision, normLength] = initialiseCaseData(geoLoc);

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
    
    selection = input('Select Region of Interest [A/B]: ', 's');

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

disp(' ');
disp(' ');


%% Initialise Lagrangian Data

maxNumCompThreads(nProc);

[dataID, LagProps, ~, ~, ...
 LagData, sampleInt, ~] = initialiseLagData(saveLoc, caseFolder, campaignID, ...
                                            caseID, cloudName, false, false, ...
                                            true, timeDirs, deltaT, ...
                                            timePrecision, maxNumCompThreads);

if strcmp(format, 'C') && ~strcmp(campaignID, 'Windsor_fullScale')
    disp(' ');
    
    disp('WARNING: Far-Wake Data Is Unavailable for This Case');
    disp('         Performing Analysis on Mid-Wake Data Instead');
    
    format = 'B';
end
    

disp(' ');
disp(' ');


%% Select Field Options

disp('Field Options');
disp('--------------');

if strcmp(campaignID, 'Windsor_fullScale')
    dLimsDefault = [20; 400]; % um
else
    dLimsDefault = [1; 147];
end

valid = false;
while ~valid
    disp(' ');
    
    selection = input('Filter Particle Diameters? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        dLims = dLimsDefault;

        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        dLims = zeros([2,1]);
        
        dLims(1) = inputD('Min');

        if dLims(1) == -1
            continue;
        end

        dLims(2) = inputD('Max');

        if dLims(2) == -1
            continue;
        end
        
        dLims = sort(dLims);
        dLims(1) = floor(dLims(1));
        dLims(2) = ceil(dLims(2));
        
        if dLims(2) < dLimsDefault(1) || dLims(1) > dLimsDefault(2)
            disp('        WARNING: No Lagrangian Data in Diameter Range');
            continue;
        end
        
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid dLimsDefault;

% Generate Dataset ID
dataID = [dataID, '_D', num2str(dLims(1)), '_D', num2str(dLims(2))];

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

nTimes = height(LagData.time);

% Identify Empty Time Instances Not Related to Case Initialisation
emptyTimes = cellfun(@isempty, LagData.origId);
firstValidTime = find(emptyTimes == false, 1, 'first');

suspectTimes = find(emptyTimes(firstValidTime:end) == true) + (firstValidTime - 1);

if ~isempty(suspectTimes)
    
    for i = 1:height(suspectTimes)
        disp(['        WARNING: Time ''', num2str(LagData.time(suspectTimes(i))), ''' Is Unexpectedly Empty']);
    end
    clear i;
    
end

% Adjust Data Origin
if strcmp(campaignID, 'Windsor_Upstream_2023')
    
    for i = 1:nTimes
        
        if ~isempty(LagData.positionCartesian{i})
            LagData.positionCartesian{i}(:,1) = LagData.positionCartesian{i}(:,1) + 1.325;
        end
        
    end
    clear i;
    
end

% Specify Region Boundaries
switch format
    
    case 'A' % 1 L
        xLimsData = [-0.537116858237548; 1.462883141762452] * normLength;
        yLimsData = [-0.4; 0.4] * normLength;
        zLimsData = [0; 0.4] * normLength;
        
    case 'B' % 2 L
        xLimsData = [-0.537116858237548; 2.462883141762452] * normLength;
        yLimsData = [-0.5; 0.5] * normLength;
        zLimsData = [0; 0.5] * normLength;
        
    case 'C' % 4 L
        xLimsData = [-0.537116858237548; 4.462883141762452] * normLength;
        yLimsData = [-0.6; 0.6] * normLength;
        zLimsData = [0; 0.6] * normLength;

end

disp(' ');

% Collate Particles of Interest
disp('    Collating Particles of Interest...');

% Initialise Progress Bar
wB = waitbar(0, 'Collating Particles of Interest', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';
dQ = parallel.pool.DataQueue;
afterEach(dQ, @parforWaitBar);
parforWaitBar(wB, nTimes);

% Perform Collation
index = cell(nTimes,1);

d = LagData.d;
positionCartesian = LagData.positionCartesian;
parfor i = 1:nTimes
    
    if ~isempty(positionCartesian{i})
        index{i} = find(((d{i} * 1e6) >= dLims(1)) & ...
                        ((d{i} * 1e6) <= dLims(2)) & ...
                        (positionCartesian{i}(:,1) >= xLimsData(1)) & ...
                        (positionCartesian{i}(:,1) <= xLimsData(2)) & ...
                        (positionCartesian{i}(:,2) >= yLimsData(1)) & ...
                        (positionCartesian{i}(:,2) <= yLimsData(2)) & ...
                        (positionCartesian{i}(:,3) >= zLimsData(1)) & ...
                        (positionCartesian{i}(:,3) <= zLimsData(2))); %#ok<PFBNS>
    end
    
    % Remove Unnecessary Data
    d{i} = -1;
    positionCartesian{i} = -1;
    
    % Update Waitbar
    send(dQ, []);
end
clear d positionCartesian;
        
delete(wB);

% Remove Unnecessary Data
disp('        Removing Unnecessary Data...');

LagFields = fieldnames(LagData);
reqFields = {'time'; 'd'; 'nParticle'; 'positionCartesian'; 'U'; 'Uslip'};

LagData = rmfield(LagData, LagFields(~ismember(LagFields, reqFields)));
LagProps = LagProps(ismember(LagProps, reqFields));

for i = 1:nTimes
    
    for j = 1:height(LagProps)
        LagData.(LagProps{j}){i} = LagData.(LagProps{j}){i}(index{i},:);
    end
    clear j;

end
clear i;

disp(' ');

% Generate Instantaneous Volume Field
disp('    Generating Presentation Grid...');

% Set Target Spatial Resolution
if strcmp(campaignID, 'Windsor_fullScale')
    cellSize.target = 32e-3;
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

volumeData.positionGrid = [x(:), y(:), z(:)]; clear x y z;

nCells = height(volumeData.positionGrid);

% Assign Particles to Mesh Nodes
disp('        Assigning Particles to Grid Cells');

% Initialise Progress Bar
wB = waitbar(0, 'Assigning Particles to Grid Cells', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';
dQ = parallel.pool.DataQueue;
afterEach(dQ, @parforWaitBar);
parforWaitBar(wB, nTimes);

% Perform Assignment
totalParticles = cellfun(@height, LagData.positionCartesian);
index = cell(nTimes,1); % Array Position of Nearest Mesh Node

xVals = unique(volumeData.positionGrid(:,1));
yVals = unique(volumeData.positionGrid(:,2));
zVals = unique(volumeData.positionGrid(:,3));
positionCartesian = LagData.positionCartesian;
positionGrid = volumeData.positionGrid;
parfor i = 1:nTimes
    
    if totalParticles(i) > 0
        index3D = zeros([totalParticles(i),3], 'uint32');
        positions = zeros([totalParticles(i),3]);
        
        for j = 1:totalParticles(i)
            [~, index3D(j,1)] = min(abs(positionCartesian{i}(j,1) - xVals));
            [~, index3D(j,2)] = min(abs(positionCartesian{i}(j,2) - yVals));
            [~, index3D(j,3)] = min(abs(positionCartesian{i}(j,3) - zVals));
            
            positions(j,:) = [xVals(index3D(j,1)), yVals(index3D(j,2)), zVals(index3D(j,3))];
        end
        
        [~, index{i}] = ismember(positions, positionGrid, 'rows');
    end
    
    % Remove Unnecessary Data
    positionCartesian{i} = -1;
    
    % Update Waitbar
    send(dQ, []);
end
clear xVals yVals zVals positionCartesian positionGrid;

delete(wB);

disp(' ');

% Generate Instantaneous Volume Field
disp('    Generating Instantaneous Volume Field...');

volumeData.time = LagData.time; LagData.time = -1;

% Calculate Instantaneous Field Variables
disp('        Calculating Instantaneous Field Variables');

% Initialise Progress Bar
wB = waitbar(0, 'Calculating Instantaneous Field Variables', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';
dQ = parallel.pool.DataQueue;
afterEach(dQ, @parforWaitBar);
parforWaitBar(wB, nTimes);

% Perform Calculation
nParticles = cell(nTimes,1); nParticles(:) = {zeros([nCells,1])}; % Number of Particles in Cell
density = nParticles; % Spray Density in Cell
d10 = nParticles; % Mean Particle Diameter in Cell
Up = cell(nTimes,1); Up(:) = {zeros([nCells,3])}; % Mean Particle Velocity in Cell
Us = cell(nTimes,1); Us(:) = {zeros([nCells,3])}; % Mean Slip Velocity in Cell
Uf = Up; % Mean Fluid Velocity in Cell

nParticle = cellfun(@double, LagData.nParticle, 'uniformOutput', false); LagData.nParticle = -1;
d = cellfun(@double, LagData.d, 'uniformOutput', false); LagData.d = -1;
U = cellfun(@double, LagData.U, 'uniformOutput', false); LagData.U = -1;
Uslip = cellfun(@double, LagData.Uslip, 'uniformOutput', false); LagData.Uslip = -1;
cellVolume = cellSize.volume;
parfor i = 1:nTimes
    
    if totalParticles(i) > 0

        for j = 1:totalParticles(i)
            nParticles{i}(index{i}(j)) = nParticles{i}(index{i}(j)) + ...
                                         nParticle{i}(j);

            density{i}(index{i}(j)) = density{i}(index{i}(j)) + ...
                                      (nParticle{i}(j) * ((1 / 12) * tau * (d{i}(j)^3)));

            d10{i}(index{i}(j)) = d10{i}(index{i}(j)) + ...
                                  (nParticle{i}(j) * d{i}(j));

            Up{i}(index{i}(j),:) = Up{i}(index{i}(j),:) + ...
                                   (nParticle{i}(j) * U{i}(j,:));
            
            Us{i}(index{i}(j),:) = Us{i}(index{i}(j),:) + ...
                                   (nParticle{i}(j) * Uslip{i}(j,:));
            
            Uf{i}(index{i}(j),:) = Uf{i}(index{i}(j),:) + ...
                                   (nParticle{i}(j) * (U{i}(j,:) - Uslip{i}(j,:)));
        end
        
        % Calculate Derived Variables
        density{i} = (1000 * density{i}) / cellVolume;
        d10{i} = d10{i} ./ nParticles{i};
        Up{i} = Up{i} ./ nParticles{i};
        Us{i} = Us{i} ./ nParticles{i};
        Uf{i} = Uf{i} ./ nParticles{i};
        
        % Set Empty Cells Back to Zero
        indexNaN = isnan(d10{i});
        
        d10{i}(indexNaN) = 0;
        Up{i}(indexNaN,:) = 0;
        Us{i}(indexNaN,:) = 0;
        Uf{i}(indexNaN,:) = 0;
    end
    
    % Make Arrays Sparse
    nParticles{i} = sparse(nParticles{i});
    density{i} = sparse(density{i});
    d10{i} = sparse(d10{i});
    Up{i} = sparse(Up{i});
    Us{i} = sparse(Us{i});
    Uf{i} = sparse(Uf{i});
    
    % Remove Unnecessary Data
    index{i} = -1;
    nParticle{i} = -1;
    d{i} = -1;
    U{i} = -1;
    Uslip{i} = -1;
    
    % Update Waitbar
    send(dQ, []);
end
clear nParticle d U Uslip cellVolume;

delete(wB);

clear LagData index;

volumeData.nParticles.inst = nParticles; clear nParticles;
volumeData.density.inst = density; clear density;
volumeData.d10.inst = d10; clear d10;
volumeData.Up.inst = Up; clear Up;
volumeData.Us.inst = Us; clear Us;
volumeData.Uf.inst = Uf; clear Uf;

% Calculate Instantaneous MCP
disp('        Calculating Instantaneous Momentum Coupling Parameter');

% Initialise Progress Bar
wB = waitbar(0, 'Calculating Instantaneous Momentum Coupling Parameter', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';
dQ = parallel.pool.DataQueue;
afterEach(dQ, @parforWaitBar);
parforWaitBar(wB, nTimes);

% Perform Calculation
Pi = cell(nTimes,1); nParticles(:) = {zeros([nCells,1])}; % Mean MCP in Cell

density = volumeData.density.inst;
Up = volumeData.Up.inst;
Uf = volumeData.Uf.inst;
d10 = volumeData.d10.inst;
cellVolume = cellSize.volume;
cellSizeY = cellSize.y;
cellSizeZ = cellSize.z;
parfor i = 1:nTimes
    volFrac = ((density{i} * cellVolume) / 1000) / cellVolume;
    
    UpMag = Up{i}; UpMag = sqrt(UpMag(:,1).^2 + UpMag(:,2).^2 + UpMag(:,3).^2);
    UfMag = Uf{i}; UfMag = sqrt(UfMag(:,1).^2 + UfMag(:,2).^2 + UfMag(:,3).^2);
    UsMag = Up{i} - Uf{i}; UsMag = sqrt(UsMag(:,1).^2 + UsMag(:,2).^2 + UsMag(:,3).^2);
    
    Re_p = (1.269 * UsMag .* d10{i}) / 1.753758e-5;
    psi = 3 * ((0.158^(1/2) * Re_p.^(1/3)) - atan(0.158^(1/2) * Re_p.^(1/3))) ./ (0.158^(3/2) * Re_p);
    Stk = psi .* ((1000 * d10{i}.^2 * 40) / (18 * 1.753758e-5 * 1.044));
    
    Pi{i} = ((1000 * volFrac) * (cellSizeY * cellSizeZ) .* UpMag) ./ ...
            (((1.269 * (1 - volFrac)) * (cellSizeY * cellSizeZ) .* UfMag) .* (1 + Stk));
    
    % Set Empty Cells Back to Zero
    Pi{i}(isnan(Pi{i})) = 0;
    
    % Make Array Sparse
    Pi{i} = sparse(Pi{i});
    
    % Remove Unnecessary Data
    density{i} = -1;
    Up{i} = -1;
    Uf{i} = -1;
    d10{i} = -1;
    
    % Update Waitbar
    send(dQ, []);
end
clear density Up Uf d10 cellVolume cellSizeY cellSizeZ;

delete(wB);

volumeData.Pi.inst = Pi; clear Pi;

% Calculate Instantaneous Particle Kinetic Energy
disp('        Calculating Instantaneous Particle Kinetic Energy');

% Initialise Progress Bar
wB = waitbar(0, 'Calculating Instantaneous Particle Kinetic Energy', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';
dQ = parallel.pool.DataQueue;
afterEach(dQ, @parforWaitBar);
parforWaitBar(wB, nTimes);

% Perform Calculation
Ek = cell(nTimes,1); nParticles(:) = {zeros([nCells,1])}; % Mean Kinetic Energy in Cell

Up = volumeData.Up.inst;
d10 = volumeData.d10.inst;
parfor i = 1:nTimes
    mass = 1000 * ((1 / 12) * tau * (d10{i}.^3));
    
    UpMag = Up{i}; UpMag = sqrt(UpMag(:,1).^2 + UpMag(:,2).^2 + UpMag(:,3).^2);
    
    Ek{i} = 0.5 * (mass .* (UpMag.^2));
    
    % Set Empty Cells Back to Zero
    Ek{i}(isnan(Ek{i})) = 0;
    
    % Make Array Sparse
    Ek{i} = sparse(Ek{i});
    
    % Remove Unnecessary Data
    Up{i} = -1;
    d10{i} = -1;
    
    % Update Waitbar
    send(dQ, []);
end
clear Up d10;

delete(wB);

volumeData.Ek.inst = Ek; clear Ek;

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


%% Save Volume Field Data

% disp('Data Save Options');
% disp('------------------');
% 
% valid = false;
% while ~valid
%     disp(' ');
%     
%     selection = input('Save Data for Future Use? [y/n]: ', 's');
%     
%     if selection == 'n' | selection == 'N' %#ok<OR2>
%         valid = true;
%     elseif selection == 'y' | selection == 'Y' %#ok<OR2>
%         
%         % Save Data
%         switch format
%             
%             case 'A'
%                 
%                 if ~exist([saveLoc, '/Numerical/MATLAB/volumeField/', campaignID, '/', caseID, '/nearWake'], 'dir')
%                     mkdir([saveLoc, '/Numerical/MATLAB/volumeField/', campaignID, '/', caseID, '/nearWake']);
%                 end
%                 
%                 disp(['    Saving to: ', saveLoc, '/Numerical/MATLAB/volumeField/', campaignID, '/', caseID, '/nearWake/', dataID, '.mat']);
%                 
%                 save([saveLoc, '/Numerical/MATLAB/volumeField/', campaignID, '/', caseID, '/nearWake/', dataID, '.mat'], ...
%                      'campaignID', 'caseID', 'dataID', 'volumeData', 'cellSize', 'sampleInt', 'timePrecision', 'dLims', '-v7.3', '-noCompression');
%                 
%                 disp('        Success');
%                 
%             case 'B'
%                 
%                 if ~exist([saveLoc, '/Numerical/MATLAB/volumeField/', campaignID, '/', caseID, '/midWake'], 'dir')
%                     mkdir([saveLoc, '/Numerical/MATLAB/volumeField/', campaignID, '/', caseID, '/midWake']);
%                 end
%                 
%                 disp(['    Saving to: ', saveLoc, '/Numerical/MATLAB/volumeField/', campaignID, '/', caseID, '/midWake/', dataID, '.mat']);
%                 
%                 save([saveLoc, '/Numerical/MATLAB/volumeField/', campaignID, '/', caseID, '/midWake/', dataID, '.mat'], ...
%                      'campaignID', 'caseID', 'dataID', 'volumeData', 'cellSize', 'sampleInt', 'timePrecision', 'dLims', '-v7.3', '-noCompression');
%                 
%                 disp('        Success');
%                 
%             case 'C'
%                 
%                 if ~exist([saveLoc, '/Numerical/MATLAB/volumeField/', campaignID, '/', caseID, '/farWake'], 'dir')
%                     mkdir([saveLoc, '/Numerical/MATLAB/volumeField/', campaignID, '/', caseID, '/farWake']);
%                 end
%                 
%                 disp(['    Saving to: ', saveLoc, '/Numerical/MATLAB/volumeField/', campaignID, '/', caseID, '/farWake/', dataID, '.mat']);
%                 
%                 save([saveLoc, '/Numerical/MATLAB/volumeField/', campaignID, '/', caseID, '/farWake/', dataID, '.mat'], ...
%                      'campaignID', 'caseID', 'dataID', 'volumeData', 'cellSize', 'sampleInt', 'timePrecision', 'dLims', '-v7.3', '-noCompression');
%                 
%                 disp('        Success');
%                 
%         end
%         
%         valid = true;
%     else
%         disp('    WARNING: Invalid Entry');
%     end
% 
% end
% clear valid;
% 
% disp(' ');
% disp(' ');


%% Select Presentation Options

disp('Presentation Options');
disp('---------------------');

valid = false;
while ~valid
    disp(' ');
    
    selection = input('Plot Instantaneous MCP? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        plotMCP = false;
        
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        plotMCP = true;
        
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;

valid = false;
while ~valid
    disp(' ');
    
    selection = input('Plot Instantaneous Uslip? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        plotSlip = false;
        
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        plotSlip = true;
        
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;

valid = false;
while ~valid
    disp(' ');
    
    selection = input('Plot Instantaneous Kinetic Energy? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        plotEk = false;
        
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        plotEk = true;
        
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;


if plotMCP || plotSlip || plotEk
    
    valid = false;
    while ~valid
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
    end
    clear valid;
    
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

        volumeData.positionGrid = volumeData.positionGrid / normLength;
    end

end

disp(' ');
disp(' ');


%% Present Volume Fields

disp('Volume Field Presentation');
disp('--------------------------');

disp(' ');

if plotMCP || plotSlip || plotEk
    gridShape = [height(unique(volumeData.positionGrid(:,1))), ...
                 height(unique(volumeData.positionGrid(:,2))), ...
                 height(unique(volumeData.positionGrid(:,3)))];
             
    spatialRes = cellSize.target / 2;
    xInit = reshape(volumeData.positionGrid(:,1), gridShape);
    yInit = reshape(volumeData.positionGrid(:,2), gridShape);
    zInit = reshape(volumeData.positionGrid(:,3), gridShape);
    POD = false;
    
    if strcmp(campaignID, 'Windsor_fullScale')
        
        if strcmp(caseID, 'Windsor_SB_fullScale_multiPhase_uncoupled')
            cMap = graphColours(4);
        elseif strcmp(caseID, 'Windsor_SB_fullScale_multiPhase_coupled')
            cMap = graphColours(5);
        elseif strcmp(caseID, 'Windsor_SB_fullScale_multiPhase_halfTread')
            cMap = graphColours(6);
        elseif strcmp(caseID, 'Windsor_SB_fullScale_multiPhase_20deg')
            cMap = graphColours(7);
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
    
    switch format

        case 'A' % 1 L
            xLimsPlot = [-0.637116858237548; 1.562883141762452];
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

if plotMCP
    disp('Presenting Instantaneous MCP...');
    
    nSurfaces = 1;
    surfaceNo = 1;
    isoValue = 0.01; % 0.01;
    figTitle = '{ }'; % Leave Blank ('{ }') for Formatting Purposes
    viewAngle = [30, 30];
    multiView = false;
    
    for i = 1:height(isoValue)
        figHold = fig;
    
        for j = startFrame:endFrame

            if j ~= startFrame
                clf(fig);
                fig = figHold;
            end
            
            fieldData = reshape(full(volumeData.Pi.inst{j}), gridShape);
            figTime = num2str(volumeData.time(j), ['%.', num2str(timePrecision), 'f']);
            figName = ['Local_Mean_Momentum_Coupling_Parameter_', num2str(isoValue(i)), '_T'...
                       erase(figTime, '.'), '_', caseID];
            figTitle = ['{', figTime, ' \it{s}}'];
            
            [fig, surfaceNo] = plotVolumeField(xLimsData, yLimsData, zLimsData, spatialRes, ...
                                               xInit, yInit, zInit, POD, fieldData, nSurfaces, surfaceNo, ...
                                               fig, figName, geometry, isoValue(i), cMap, figTitle, viewAngle, ...
                                               multiView, xLimsPlot, yLimsPlot, zLimsPlot, figSave);

        end
        clear j;
        
    end
    clear i;
    
    disp(' ');
end

if plotSlip
    disp('Presenting Instantaneous Uslip...');
    
    nSurfaces = 1;
    surfaceNo = 1;
    isoValue = 0.125;
    figTitle = '{ }'; % Leave Blank ('{ }') for Formatting Purposes
    viewAngle = [30, 30];
    multiView = false;
    
    for i = 1:height(isoValue)
        figHold = fig;
    
        for j = startFrame:endFrame

            if j ~= startFrame
                clf(fig);
                fig = figHold;
            end
            
            fieldData = full(volumeData.Us.inst{j}) / 22.22;
            fieldData = sqrt(fieldData(:,1).^2 + fieldData(:,2).^2 + fieldData(:,3).^2);
            fieldData = reshape(fieldData, gridShape);
            
            figTime = num2str(volumeData.time(j), ['%.', num2str(timePrecision), 'f']);
            figName = ['Local_Mean_Slip_Velocity_', num2str(isoValue(i)), '_T'...
                       erase(figTime, '.'), '_', caseID];
            figTitle = ['{', figTime, ' \it{s}}'];
            
            [fig, surfaceNo] = plotVolumeField(xLimsData, yLimsData, zLimsData, spatialRes, ...
                                               xInit, yInit, zInit, POD, fieldData, nSurfaces, surfaceNo, ...
                                               fig, figName, geometry, isoValue(i), cMap, figTitle, viewAngle, ...
                                               multiView, xLimsPlot, yLimsPlot, zLimsPlot, figSave);

        end
        clear j;
        
    end
    clear i;
    
    disp(' ');
end

if plotEk
    disp('Presenting Instantaneous Kinetic Energy...');
    
    nSurfaces = 1;
    surfaceNo = 1;
    isoValue = 1e-6;
    figTitle = '{ }'; % Leave Blank ('{ }') for Formatting Purposes
    viewAngle = [30, 30];
    multiView = false;
    
    for i = 1:height(isoValue)
        figHold = fig;
    
        for j = startFrame:endFrame

            if j ~= startFrame
                clf(fig);
                fig = figHold;
            end
            
            fieldData = reshape(full(volumeData.Ek.inst{j}), gridShape);            
            figTime = num2str(volumeData.time(j), ['%.', num2str(timePrecision), 'f']);
            figName = ['Local_Mean_Kinetic_Energy_', num2str(isoValue(i)), '_T'...
                       erase(figTime, '.'), '_', caseID];
            figTitle = ['{', figTime, ' \it{s}}'];
            
            [fig, surfaceNo] = plotVolumeField(xLimsData, yLimsData, zLimsData, spatialRes, ...
                                               xInit, yInit, zInit, POD, fieldData, nSurfaces, surfaceNo, ...
                                               fig, figName, geometry, isoValue(i), cMap, figTitle, viewAngle, ...
                                               multiView, xLimsPlot, yLimsPlot, zLimsPlot, figSave);

        end
        clear j;
        
    end
    clear i;
    
    disp(' ');
end

if ~plotMCP && ~plotSlip && ~plotEk
    disp('Skipping Volume Field Presentation...');

    disp(' ');
end


%% Local Functions

function D = inputD(type)

    D = str2double(input(['    ', type, ' Diameter of Interest [', char(956), 'm]: '], 's'));
    
    if isnan(D) || length(D) > 1 || D < 1
        disp('        WARNING: Invalid Entry');
        
        D = -1;
    end
    
end


function frameNo = inputFrames(Nt, type)

    frameNo = str2double(input(['    Input Desired ', type, ' Frame [1-', num2str(Nt), ']: '], 's'));
    
    if isnan(frameNo) || frameNo < 1 || frameNo > Nt
        disp('        WARNING: Invalid Entry');
        
        frameNo = -1;
    end

end        
