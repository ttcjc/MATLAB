%% Planar Pressure POD Calculator v1.0

clear variables;
close all;
clc;
evalc('delete(gcp(''nocreate''));');

saveLocation = '/mnt/Processing/Data';
% saveLocation = '~/Data';

normalise = true; % Normalisation of Dimensions

nProc = maxNumCompThreads - 2; % Number of Processors Used for Parallel Collation

fig = 0; % Initialise Figure Tracking
figHold = 0; % Enable Overwriting of Figures

disp('===================================');
disp('Planar Pressure POD Calculator v1.0');
disp('===================================');

disp(' ');
disp(' ');


%% Changelog

% v1.0 - Initial Commit


%% Initialisation

[caseName, dataID, probeData, sampleInterval, timePrecision, geometry, ...
 xDims, yDims, zDims, spacePrecision] = initialisePressureProbeData(saveLocation, normalise, nProc);    

if normalise
    dataID = [dataID, '_Norm'];
end

disp(' ');
disp(' ');


%% Perform Planar POD (Snapshot Method)

disp('Planar Proper Orthogonal Decomposition');
disp('---------------------------------------');

disp(' ');

disp('***********');
disp('  RUNNING ');

tic;

disp(' ');

disp('    Initialising...');

% Initialise POD Variables
PODdata.positionGrid = probeData.position;
PODdata.time = probeData.time;
PODdata.p.mean = probeData.pMean;
PODdata.p.inst = probeData.p;
PODdata.p.prime = probeData.pPrime;

clear probeData;

% Adjust Data Origin
if contains(caseName, 'Run_Test') || (contains(caseName, 'Windsor') && contains(caseName, 'Upstream'))
    
    if normalise
        PODdata.positionGrid(:,1) = PODdata.positionGrid(:,1) + round((1.325 / 1.044), spacePrecision);
    else
        PODdata.positionGrid(:,1) = PODdata.positionGrid(:,1) + 1.325; %#ok<*UNRCH>
    end
    
end

% Specify Map Boundaries
parts = fieldnames(geometry);
for i = 1:height(parts)
    
    if max(geometry.(parts{i}).vertices(:,1)) == xDims(2)
        break;
    end
    
    if i == height(parts)
        error('Mismatch Between ''xDims'' and Geometry Bounding Box')
    end
    
end

geoPoints = geometry.(parts{i}).vertices;
basePoints = geoPoints((geoPoints(:,1) == xDims(2)),:);

mapPerim = boundary(basePoints(:,2), basePoints(:,3), 0.95);
mapPerim = basePoints(mapPerim,:);
basePoly = polyshape(mapPerim(:,2), mapPerim(:,3), 'keepCollinearPoints', true);
basePoly = polybuffer(basePoly, -0.0025, 'jointType', 'square');
mapPerim = ones(height(basePoly.Vertices),3) * mapPerim(1,1);
mapPerim(:,[2,3]) = basePoly.Vertices(:,[1,2]);

if ~all(mapPerim(1,:) == mapPerim(end,:))
    mapPerim = [mapPerim; mapPerim(1,:)]; % Close Boundary
end

clear basePoints basePoly;

xLimsData = xDims(2);
yLimsData = [min(mapPerim(:,2)); max(mapPerim(:,2))];
zLimsData = [min(mapPerim(:,3)); max(mapPerim(:,3))];

% Shift Data off Base
PODdata.positionGrid(:,1) = PODdata.positionGrid(:,1) + 1e-3;

disp(' ');

% Perform Planar Snapshot POD
[fig, PODdata, modesEnergetic, modes80percent, Ns, Nt] = performPOD(fig, PODdata, 'p', ...
                                                                    'scalar', 'Base');

executionTime = toc;

disp(' ');

disp(['    Run Time: ', num2str(executionTime), 's']);

disp(' ');

disp('  SUCCESS  ');
disp('***********');

disp(' ');
disp(' ');


%% Select Mode Presentation Options

disp('Mode Presentation Options');
disp('--------------------------');

valid = false;
while ~valid
    disp(' ');
    selection = input('Plot Individual Modes? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        plotModes = false;
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        plotModes = true;
        nModes = inputModes(Nt);
        
        if nModes == -1
            continue;
        else
            nModes = sort(nModes);
        end
        
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;

disp(' ');
disp(' ');


%% Present POD Modes

disp('Mode Presentation');
disp('------------------');

disp(' ');

if ~isempty(plotModes)
    % Define Plot Limits
    orientation = 'YZ';
    
    if contains(caseName, ["Run_Test", "Windsor"])
        xLimsPlot = [0.31875; 1.52725];
        yLimsPlot = [-0.2445; 0.2445];
        zLimsPlot = [0; 0.389];
    end

    if normalise
        xLimsPlot = round((xLimsPlot / 1.044), spacePrecision);
        yLimsPlot = round((yLimsPlot / 1.044), spacePrecision);
        zLimsPlot = round((zLimsPlot / 1.044), spacePrecision);
    end

    positionData = PODdata.positionGrid;
    cMap = cool2warm(24);
    figTitle = '-'; % Leave Blank ('-') for Formatting Purposes
    cLims = [-1; 1];

    for i = nModes
        disp(['    Presenting Mode #', num2str(i), '...']);

        scalarData = rescale(PODdata.phi_mode(:,i), -1, 1);
        figName = ['Pressure_Planar_POD_M', num2str(i)];
        CoM = [];
        figSubtitle = [num2str(round(PODdata.modeEnergy(i), 2), '%.2f'), '\it{%}'];

        fig = planarScalarPlots(orientation, xLimsData, yLimsData, zLimsData, positionData, scalarData, ...
                                mapPerim, fig, figName, cMap, geometry, xDims, yDims, zDims, ...
                                CoM, figTitle, figSubtitle, cLims, xLimsPlot, yLimsPlot, zLimsPlot, normalise);
    end
    
else
    disp('    Skipping Mode Presentation');
end

disp(' ');
disp(' ');


%% Save POD Data

disp('Data Save Options');
disp('------------------');

valid = false;
while ~valid
    disp(' ');
    selection = input('Save Data for Future Use? [y/n]: ', 's');
    
    if selection == 'n' | selection == 'N' %#ok<OR2>
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        
        if ~exist([saveLocation, '/Numerical/MATLAB/planarPressurePOD/', caseName], 'dir')
            mkdir([saveLocation, '/Numerical/MATLAB/planarPressurePOD/', caseName]);
        end
        
        disp(['    Saving to: ', saveLocation, '/Numerical/MATLAB/planarPressurePOD/', caseName, '/', dataID, '.mat']);
        save([saveLocation, '/Numerical/MATLAB/planarPressurePOD/', caseName, '/', dataID, '.mat'], ...
             'dataID', 'PODdata', 'sampleInterval', 'normalise', '-v7.3', '-noCompression');
        disp('        Success');
        
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;

disp(' ');
disp(' ');


%% Select Reconstruction Options

disp('Reconstruction Options');
disp('-----------------------');

valid = false;
while ~valid
    disp(' ');
    selection = input('Perform Field Reconstruction Using N Modes? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        return;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        nModes = inputModes(Nt);
        
        if nModes == -1
            continue;
        else
            nModes = sort(nModes);
        end
        
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;

disp(' ');
disp(' ');


%% Perform Field Reconstruction

disp('Field Reconstruction');
disp('---------------------');

disp(' ');

disp('***********');
disp('  RUNNING ');

tic;

disp(' ');

disp('    Initialising...');

% Initialise Reconstruction Variables
reconData.positionGrid = PODdata.positionGrid;
reconData.time = PODdata.time;
reconData.p.mean = PODdata.p.mean;
reconData.p.inst = cell(Nt,1);

for i = 1:Nt
    reconData.p.inst{i} = reconData.p.mean;
end

disp(' ');

% Perform Field Reconstruction
reconData = reconstructPOD(reconData, PODdata, 'p', nModes, Ns, Nt, 'scalar', true);

% Calculate Pressure Coefficient
if contains(caseName, 'Run_Test') || (contains(caseName, 'Windsor') && contains(caseName, 'Upstream'))
    U = 40; % m/s
    rho = 1.269; % kg/m^3
    pRef = 0 * rho; % Pa
    
    reconData.p.mean = (reconData.p.mean - pRef) / (0.5 * rho * U^2);
    
    for i = 1:height(reconData.time)
        reconData.p.inst{i} = (reconData.p.inst{i} - pRef) / (0.5 * rho * U^2);
    end
    
end

% Perform Blockage Correction
if contains(caseName, 'Run_Test') || (contains(caseName, 'Windsor') && contains(caseName, 'Upstream'))
    Am = (0.289 * 0.389) + (2 * (0.046 * 0.055));
    At = (2 * (0.9519083 + (3.283 * tan(atan(0.0262223 / 9.44)))) * 1.32);
    
    reconData.p.mean = (reconData.p.mean + (2 * (Am / At))) / (1 + (2 * (Am / At)));
    reconData.p.mean = (reconData.p.mean + (2 * (Am / At))) / (1 + (2 * (Am / At)));
    
    for i = 1:height(reconData.time)
        reconData.p.inst{i} = (reconData.p.inst{i} + (2 * (Am / At))) / (1 + (2 * (Am / At)));
    end
    
end

% Initialise Progress Bar
wB = waitbar(0, 'Calculating Reconstructed Centre of Pressure', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';

% Calculate Reconstructed CoP
reconData.p.CoP = cell(Nt,1);

for i = 1:Nt
    reconData.p.CoP{i} = zeros(1,3);

    reconData.p.CoP{i}(1) = reconData.positionGrid(1,1);
    reconData.p.CoP{i}(2) = sum(reconData.p.inst{i} .* reconData.positionGrid(:,2)) / sum(reconData.p.inst{i});
    reconData.p.CoP{i}(3) = sum(reconData.p.inst{i} .* reconData.positionGrid(:,3)) / sum(reconData.p.inst{i});
    
    waitbar((i / Nt), wB);
end

delete(wB);

executionTime = toc;

disp(' ');

disp(['    Run Time: ', num2str(executionTime), 's']);

disp(' ');

disp('  SUCCESS  ');
disp('***********');

disp(' ');
disp(' ');


%% Select Reconstruction Presentation Options

disp('Reconstruction Presentation Options');
disp('------------------------------------');

valid = false;
while ~valid
    disp(' ');
    selection = input('Plot Reconstructed Field? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        plotRecon = false;
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        plotRecon = true;
        nFrames = inputFrames(Nt);
        
        if nFrames == -1
            continue;
        end
        
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;

disp(' ');
disp(' ');


%% Present Reconstruction

disp('Reconstruction Presentation');
disp('----------------------------');

disp(' ');

if plotRecon
    disp('    Presenting Reconstructed Field...');
    
    % Define Plot Limits
    orientation = 'YZ';
    
    if contains(caseName, ["Run_Test", "Windsor"])
        xLimsPlot = [0.31875; 1.52725];
        yLimsPlot = [-0.2445; 0.2445];
        zLimsPlot = [0; 0.389];
    end

    if normalise
        xLimsPlot = round((xLimsPlot / 1.044), spacePrecision);
        yLimsPlot = round((yLimsPlot / 1.044), spacePrecision);
        zLimsPlot = round((zLimsPlot / 1.044), spacePrecision);
    end

    positionData = reconData.positionGrid;
    cMap = viridis(24);
    figTitle = '-'; % Leave Blank ('-') for Formatting Purposes
    cLims = [-0.235; -0.023]; % Base Cp

    figHold = fig;

    for i = 1:nFrames
        
        if i ~= 1
            clf(fig);
            fig = figHold;
        end
        
        scalarData = reconData.p.inst{i};
        figTime = num2str(reconData.time(i), ['%.', num2str(timePrecision), 'f']);
        figName = ['Cp_Reconstruction_T', erase(figTime, '.')];
        CoM = reconData.p.CoP{i};
        figSubtitle = [num2str(reconData.time(i), ['%.', num2str(timePrecision), 'f']), ' \it{s}'];
        
        fig = planarScalarPlots(orientation, xLimsData, yLimsData, zLimsData, positionData, scalarData, ...
                                mapPerim, fig, figName, cMap, geometry, xDims, yDims, zDims, ...
                                CoM, figTitle, figSubtitle, cLims, xLimsPlot, yLimsPlot, zLimsPlot, normalise);
    end
    
else
    disp('    Skipping Reconstruction Presentation');
end

disp(' ');
disp(' ');


%% Save Reconstruction Data

disp('Data Save Options');
disp('------------------');

valid = false;
while ~valid
    disp(' ');
    selection = input('Save Data for Future Use? [y/n]: ', 's');
    
    if selection == 'n' | selection == 'N' %#ok<OR2>
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        
        if ~exist([saveLocation, '/Numerical/MATLAB/planarPressureReconstruction/', caseName], 'dir')
            mkdir([saveLocation, '/Numerical/MATLAB/planarPressureReconstruction/', caseName]);
        end
        
        disp(['    Saving to: ', saveLocation, '/Numerical/MATLAB/planarPressureReconstruction/', caseName, '/', dataID, '_', mat2str(nModes), '.mat']);
        save([saveLocation, '/Numerical/MATLAB/planarPressureReconstruction/', caseName, '/', dataID, '_', mat2str(nModes), '.mat'], ...
             'dataID', 'reconData', 'nModes', 'sampleInterval', 'normalise', '-v7.3', '-noCompression');
         disp('        Success');
        
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;


%% Local Functions

function modes = inputModes(nModes)

    modes = str2num(input('    Input Desired Modes [Row Vector Form]: ', 's')); %#ok<ST2NM>
    
    if isempty(modes) || any(isnan(modes)) || ~isrow(modes) > 1 || any(modes <= 0) || any(modes > nModes)
        disp('        WARNING: Invalid Entry');
        modes = -1;
    end

end


function nFrames = inputFrames(Nt)

    nFrames = str2double(input(['    Input Desired Frame Count [1-', num2str(Nt), ']: '], 's'));
    
    if isnan(nFrames) || nFrames <= 0 || nFrames > Nt
        disp('        WARNING: Invalid Entry');
        nFrames = -1;
    end

end