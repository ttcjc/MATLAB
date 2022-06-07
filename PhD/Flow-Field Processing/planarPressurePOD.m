%% Planar Pressure POD Calculator v1.0

clear variables;
close all;
clc;
evalc('delete(gcp(''nocreate''));');

normalise = true; % Normalisation of Dimensions

nProc = maxNumCompThreads - 2; % Number of Processors Used for Parallel Collation

fig = 0; % Initialise Figure Tracking
figHold = 0; % Enable Overwriting of Figures

disp ('===================================');
disp ('Planar Pressure POD Calculator v1.0');
disp ('===================================');

disp (' ');
disp (' ');


%% Changelog

% v1.0 - Initial Commit


%% Initialisation

[caseName, probeData, timePrecision, geometry, ...
 xDims, yDims, zDims, spacePrecision] = initialisePressureProbeData(normalise, nProc);    

disp(' ');
disp(' ');


%% Perform Planar POD (Snapshot Method)

disp('Planar Proper Orthogonal Decomposition');
disp('---------------------------------------');

disp(' ');

disp('***********');
disp('  RUNNING ');

tic;
evalc('parpool(nProc);');

disp(' ');

disp('    Initialising...');

% Specify Map Boundaries
parts = fieldnames(geometry);
for i = 1:height(parts)
    
    if max(geometry.(parts{i}).vertices(:,1)) == xDims(2)
        break
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
    mapPerim = vertcat(mapPerim, mapPerim(1,:)); % Close Boundary
end

clear basePoints basePoly;

xLimsData = xDims(2);
yLimsData = [min(mapPerim(:,2)); max(mapPerim(:,2))];
zLimsData = [min(mapPerim(:,3)); max(mapPerim(:,3))];

% Initialise POD Variables
PODdata.positionGrid = probeData.position;
PODdata.time = probeData.time;
PODdata.p.inst = probeData.p;
PODdata.p.mean = probeData.pMean;
PODdata.p.prime = probeData.pPrime;

clear probeData

% Adjust Data Origin
if contains(caseName, 'Run_Test') || (contains(caseName, 'Windsor') && contains(caseName, 'Upstream'))
    
    if normalise
        PODdata.positionGrid(:,1) = PODdata.positionGrid(:,1) + round((1.325 / 1.044), spacePrecision);
    else
        PODdata.positionGrid(:,1) = PODdata.positionGrid(:,1) + 1.325; %#ok<*UNRCH>
    end
    
end

% Shift Data off Base
PODdata.positionGrid(:,1) = PODdata.positionGrid(:,1) + 1e-3;

disp(' ');

disp('    Performing Planar POD Using the Snapshot Method...');

Ns = height(PODdata.positionGrid); % Number of Spatial Points
Nt = height(PODdata.time); % Number of Time Instances

% Initialise Progress Bar
wB = waitbar(0, 'Assembling Snapshot Matrix', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';
dQ = parallel.pool.DataQueue;
afterEach(dQ, @parforWaitBar);

parforWaitBar(wB, Nt);

% Assemble Snapshot Matrix
snapshotMatrix = zeros(Nt,Ns);

pPrime = PODdata.p.prime;
parfor i = 1:Nt
    
    for j = 1:Ns
        snapshotMatrix(i,j) = pPrime{i}(j);
    end
    
    send(dQ, []);
end
clear varPrime;

delete(wB);

PODdata.snapshotMatrix = snapshotMatrix;
clear snapshotMatrix;

% Generate Correlation Matrix
PODdata.C = (PODdata.snapshotMatrix * PODdata.snapshotMatrix') / (Nt - 1);

% Solve Eigenvalue Problem
[PODdata.A_mode, PODdata.lambda] = eig(PODdata.C, 'vector');

% Sort Eigenvalues and Eigenvalues in Descending Order
[PODdata.lambda, index] = sort(PODdata.lambda, 'descend');
PODdata.A_mode = PODdata.A_mode(:,index); % Temporal Modes

% Calculate Spatial Coefficients
PODdata.phi_coeff = PODdata.snapshotMatrix' * PODdata.A_mode;

% Normalisation to Match Direct Method
PODdata.phi_mode = normc(PODdata.phi_coeff); % Spatial Modes
PODdata.A_coeff = PODdata.snapshotMatrix * PODdata.phi_mode; % Temporal Coefficients

% Identify Mode Energy Content
PODdata.modeEnergy = (PODdata.lambda / sum(PODdata.lambda)) * 100;
modesEnergetic = height(find(PODdata.modeEnergy > 1));
modes80percent = find(cumsum(PODdata.modeEnergy) > 80, 1);

disp(' ');

disp(['    First ', num2str(modesEnergetic), ' Modes Each Contain Greater Than 1% of Total Energy']);
disp(['    First ', num2str(modes80percent), ' Modes Contain Approximately 80% of Total Energy']);

% Figure Setup
fig = fig + 1;
figName = 'Base_Pressure_Planar_POD_Energy_Content';
set(figure(fig), 'outerPosition', [25, 25, 1275, 850], 'name', figName);
set(gca, 'lineWidth', 2, 'fontName', 'LM Mono 12', ...
         'fontSize', 20, 'layer', 'top');
hold on;

% Plot
plot(PODdata.modeEnergy(1:((ceil(modesEnergetic / 10) * 10) - 1)), 'lineWidth', 1.5, 'marker', 'o', 'color', ([74, 24, 99] / 255));

% Figure Formatting
axis on;
box on;
grid off;
xlim([0; (ceil(modesEnergetic / 10) * 10)]);
ylim([0; (ceil(max(PODdata.modeEnergy)/10) * 10)]);
tickData = (0:(((ceil(modesEnergetic / 10) * 10) - 0) / 5):(ceil(modesEnergetic / 10) * 10));
xticks(tickData(2:(end - 1)));
tickData = (0:(((ceil(max(PODdata.modeEnergy)/10) * 10) - 0) / 5):(ceil(max(PODdata.modeEnergy)/10) * 10));
yticks(tickData(2:(end - 1)));
xL = xlabel({' ', '{\bf{Mode}}'}, 'fontName', 'LM Roman 12');
yL = ylabel({'{\bf{Energy Content (\it{%})}}', ' '}, 'fontName', 'LM Roman 12');
set(gca, 'outerPosition', [0.05, 0.05, 0.9, 0.9]);
hold off;

pause(2);
exportgraphics(gcf, ['~/MATLAB/Output/Figures/', figName, '.png'], 'resolution', 300);

evalc('delete(gcp(''nocreate''));');
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
        plotModes = [];
        
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        plotModes = inputModes(Nt);
        
        if plotModes == -1
            continue
        end
        
        plotModes = sort(plotModes);
        
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
    cMap = turbo(24);
    figTitle = '-'; % Leave Blank ('-') for Formatting Purposes

    for i = plotModes
        disp(['    Presenting Mode #', num2str(i), '...']);

        scalarData = rescale(PODdata.phi_mode(:,i), -1, 1);
        figName = ['Base_Pressure_Planar_POD_M', num2str(i)];
        CoP = [];
        figSubtitle = [num2str(round(PODdata.modeEnergy(i), 2), '%.2f'), '\it{%}'];
        cLims = [-1; 1];

        fig = planarScalarPlots(orientation, xLimsData, yLimsData, zLimsData, positionData, scalarData, ...
                                mapPerim, fig, figName, cMap, geometry, xDims, yDims, zDims, ...
                                CoP, figTitle, figSubtitle, cLims, xLimsPlot, yLimsPlot, zLimsPlot, normalise);
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
        
        if ~exist(['/mnt/Processing/Data/Numerical/MATLAB/planarPressurePOD/', caseName], 'dir')
            mkdir(['/mnt/Processing/Data/Numerical/MATLAB/planarPressurePOD/', caseName]);
        end
        
        save(['/mnt/Processing/Data/Numerical/MATLAB/planarPressurePOD/', caseName, '/', fileName], ...
              'PODdata', 'sampleInterval', 'normalise', '-v7.3', '-noCompression');
        disp(['    Saving to: ~/Data/Numerical/MATLAB/planarPressurePOD/', caseName, '/', fileName]);
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
        reconModes = [];
        
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        reconModes = inputModes(Nt);
        
        if reconModes == -1
            continue
        end
        
        reconModes = sort(reconModes);
        
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
evalc('parpool(nProc);');

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

% Perform Reconstruction
disp('    Performing Field Reconstruction...');

for i = reconModes
    % Initialise Progress Bar
    wB = waitbar(0, ['Adding Mode #', num2str(i), ' to Reconstruction'], 'name', 'Progress');
    wB.Children.Title.Interpreter = 'none';
    dQ = parallel.pool.DataQueue;
    afterEach(dQ, @parforWaitBar);

    parforWaitBar(wB, Nt);
    
    % Identify Mode Contribution
    mode = ['M_', num2str(i)];
    
    modeMatrix = PODdata.A_coeff(:,i) * PODdata.phi_mode(:,i)';
    pPrime = cell(Nt,1);
    
    parfor j = 1:Nt
        pPrime{j} = zeros(Ns,1);
        
        for k = 1:Ns
            pPrime{j}(k) = modeMatrix(j,k);
        end
        
        send(dQ, []);
    end
    
    delete(wB);
    
    reconData.(mode).modeMatrix = modeMatrix;
    reconData.(mode).prime = pPrime;
    clear modeMatrix varPrime;
    
    % Add Mode to Reconstruction
    for j = 1:Nt
        reconData.p.inst{j} = reconData.p.inst{j} + reconData.(mode).prime{j};
    end
    
end

% Calculate Pressure Coefficient
if contains(caseName, 'Run_Test') || (contains(caseName, 'Windsor') && contains(caseName, 'Upstream'))
    U = 40; % m/s
    rho = 1.269; % kg/m^3
    pRef = 0 * rho; % Pa
    
    reconData.Cp.mean = (reconData.p.mean - pRef) / (0.5 * rho * U^2);
    reconData.Cp.inst = cell(height(reconData.time),1);
    
    for i = 1:height(reconData.time)
        reconData.Cp.inst{i} = (reconData.p{i} - pRef) / (0.5 * rho * U^2);
    end
    
end

% Perform Blockage Correction
if contains(caseName, 'Run_Test') || (contains(caseName, 'Windsor') && contains(caseName, 'Upstream'))
    Am = (0.289 * 0.389) + (2 * (0.046 * 0.055));
    At = (2 * (0.9519083 + (3.283 * tan(atan(0.0262223 / 9.44)))) * 1.32);
    
    reconData.p.mean = (reconData.Cp.mean + (2 * (Am / At))) / (1 + (2 * (Am / At)));
    reconData.Cp.mean = (reconData.Cp.mean + (2 * (Am / At))) / (1 + (2 * (Am / At)));
    
    for i = 1:height(reconData.time)
        reconData.p.inst{i} = (reconData.Cp.inst{i} + (2 * (Am / At))) / (1 + (2 * (Am / At)));
        reconData.Cp.inst{i} = (reconData.Cp.inst{i} + (2 * (Am / At))) / (1 + (2 * (Am / At)));
    end
    
end

% Initialise Progress Bar
wB = waitbar(0, 'Calculating Reconstructed Centre of Pressure', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';
dQ = parallel.pool.DataQueue;
afterEach(dQ, @parforWaitBar);

parforWaitBar(wB, Nt);

% Calculate Reconstructed CoP
CoP = cell(Nt,1);

p = reconData.Cp.inst;
positionGrid = reconData.positionGrid;
parfor i = 1:Nt
    CoP{i} = zeros(1,3);
    CoP{i}(1) = positionGrid(1,1); %#ok<PFBNS>
    
    for j = 1:height(positionGrid)
        CoP{i}(2) = CoP{i}(2) + (p{i}(j) * positionGrid(j,2));
        CoP{i}(3) = CoP{i}(3) + (p{i}(j) * positionGrid(j,3));
    end
    
    CoP{i}(2) = CoP{i}(2) / sum(p{i});
    CoP{i}(3) = CoP{i}(3) / sum(p{i});
    send(dQ, []);
end
clear p positionGrid;

delete(wB);

reconData.p.CoP = CoP;
reconData.Cp.CoP = reconData.p.CoP;
clear CoP;

evalc('delete(gcp(''nocreate''));');
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

    figHold = fig;

    for i = 1:Nt
        
        if i ~= 1
            clf(fig);
            fig = figHold;
        end
        
        scalarData = reconData.Cp.inst{i};
        figTime = num2str(reconData.time(i), ['%.', num2str(timePrecision), 'f']);
        figName = ['Cp_Reconstruction_T', erase(figTime, '.')];
        CoM = reconData.p.CoP{i};
        figSubtitle = [num2str(reconData.time(i), ['%.', num2str(timePrecision), 'f']), ' \it{s}'];
        cLims = [-0.235; -0.023]; % Base Cp
        
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
        
        if ~exist(['/mnt/Processing/Data/Numerical/MATLAB/planarPressureReconstruction/', caseName], 'dir')
            mkdir(['/mnt/Processing/Data/Numerical/MATLAB/planarPressureReconstruction/', caseName]);
        end
        
        save(['/mnt/Processing/Data/Numerical/MATLAB/planarPressureReconstruction/', caseName, '/', fileName], ...
              'reconData', 'sampleInterval', 'normalise', '-v7.3', '-noCompression');
        disp(['    Saving to: ~/Data/Numerical/MATLAB/planarPressureReconstruction/', caseName, '/', fileName]);
        disp('        Success');
        
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;


%% Local Functions

function modes = inputModes(nModes)

    modes = str2num(input('    Input Desired Modes (Row Vector Form) [s]: ', 's')); %#ok<ST2NM>
    
    if isempty(modes) || any(isnan(modes)) || ~isrow(modes) > 1 || any(modes <= 0) || any(modes >= nModes)
        disp('        WARNING: Invalid Entry');
        modes = -1;
    end

end