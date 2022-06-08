%% Planar Velocity POD Calculator v2.0

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

% v1.0 - Initial Commit (Base Contamination and Far-Field Extraction Plane)
% v2.0 - Rewrite, Accommodating New OpenFOAM Data Formats


%% Initialisation

[caseName, dataID, probeData, sampleInterval, timePrecision, geometry, ...
 xDims, yDims, zDims, spacePrecision] = initialiseVelocityProbeData('planarPOD', normalise, nProc);

planeName = cell2mat(fieldnames(probeData));
probeData = probeData.(cell2mat(fieldnames(probeData)));

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

% Initialise POD Variables
PODdata.planeOrientation = probeData.planeOrientation;
PODdata.xLims = probeData.xLims;
PODdata.yLims = probeData.yLims;
PODdata.zLims = probeData.zLims;

PODdata.positionGrid = probeData.position;
PODdata.time = probeData.time;
PODdata.u.mean = probeData.uMean;
PODdata.v.mean = probeData.vMean;
PODdata.w.mean = probeData.wMean;
PODdata.u.inst = probeData.u;
PODdata.v.inst = probeData.v;
PODdata.w.inst = probeData.w;
PODdata.u.prime = probeData.uPrime;
PODdata.v.prime = probeData.vPrime;
PODdata.w.prime = probeData.wPrime;

clear probeData

% Adjust Data Origin
if contains(caseName, 'Run_Test') || (contains(caseName, 'Windsor') && contains(caseName, 'Upstream'))
    
    if normalise
        PODdata.positionGrid(:,1) = PODdata.positionGrid(:,1) + round((1.325 / 1.044), spacePrecision);
        PODdata.xLims = PODdata.xLims + round((1.325 / 1.044), spacePrecision);
    else
        PODdata.positionGrid(:,1) = PODdata.positionGrid(:,1) + 1.325; %#ok<*UNRCH>
        PODdata.xLims = PODdata.xLims + 1.325;
    end
    
end

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
uSnapshotMatrix = zeros(Nt,Ns);
vSnapshotMatrix = uSnapshotMatrix;
wSnapshotMatrix = uSnapshotMatrix;

uPrime = PODdata.u.prime;
vPrime = PODdata.v.prime;
wPrime = PODdata.w.prime;
parfor i = 1:Nt
    
    for j = 1:Ns
        uSnapshotMatrix(i,j) = uPrime{i}(j);
        vSnapshotMatrix(i,j) = vPrime{i}(j);
        wSnapshotMatrix(i,j) = wPrime{i}(j);
    end
    
    send(dQ, []);    
end
clear uPrime vPrime wPrime

delete(wB);

PODdata.snapshotMatrix = [uSnapshotMatrix, vSnapshotMatrix, wSnapshotMatrix];
clear uSnapshotMatrix vSnapshotMatrix wSnapshotMatrix;

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
    % Specify Default Axes Limits
    orientation = PODdata.planeOrientation;
    
    switch orientation
        
        case 'YZ'
            if contains(caseName, ["Run_Test", "Windsor"])
                xLimsPlot = [0.31875; 4.65925]; % [m]
                yLimsPlot = [-0.5945; 0.5945];
                zLimsPlot = [0; 0.739];
            end
            
        case {'XZ', 'XY'}
            
            if contains(caseName, ["Run_Test", "Windsor"])
                xLimsPlot = [0.31875; 1.075]; % [m]
                yLimsPlot = [-0.3445; 0.3445];
                zLimsPlot = [0; 0.489];
            end
            
    end
    
    if normalise
        xLimsPlot = round((xLimsPlot / 1.044), spacePrecision);
        yLimsPlot = round((yLimsPlot / 1.044), spacePrecision);
        zLimsPlot = round((zLimsPlot / 1.044), spacePrecision);
    end
    
    % Modify Axes and Data Limits Based on Format
    xLimsData = PODdata.xLims;
    yLimsData = PODdata.yLims;
    zLimsData = PODdata.zLims;
    
    xLimsPlot = [min(min(xLimsPlot), min(xLimsData)); max(max(xLimsPlot), max(xLimsData))];
    yLimsPlot = [min(min(yLimsPlot), min(yLimsData)); max(max(yLimsPlot), max(yLimsData))];
    zLimsPlot = [min(min(zLimsPlot), min(zLimsData)); max(max(zLimsPlot), max(zLimsData))];
    
    positionData = PODdata.positionGrid;
    nComponents = 1;
    
    switch orientation
        
        case 'YZ'
            component = 'u';
            
        case 'XZ'
            component = 'v';
            
        case 'XY'
            component = 'w';
    
    end
    
    cMap = turbo(24);
    streamlines = true;
    figTitle = '-'; % Leave Blank ('-') for Formatting Purposes
    cLims = [-1; 1];
    
    for i = plotModes
        disp(['    Presenting Mode #', num2str(i), '...']);
        
        vectorData = rescale([PODdata.phi_mode((1:Ns),i), ...
                              PODdata.phi_mode(((Ns + 1):(2 * Ns)),i), ...
                              PODdata.phi_mode((((2 * Ns) + 1):end),i)], -1, 1);
        figName = ['Velocity_Planar_POD_M', num2str(i)];
        figSubtitle = [num2str(round(PODdata.modeEnergy(i), 2), '%.2f'), '\it{%}'];
        
        fig = planarVectorPlots(orientation, xLimsData, yLimsData, zLimsData, positionData, ...
                                vectorData, nComponents, component, fig, figName, cMap, geometry, ...
                                streamlines, xDims, yDims, zDims, figTitle, figSubtitle, cLims, ...
                                xLimsPlot, yLimsPlot, zLimsPlot, normalise);
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
        
        if ~exist(['/mnt/Processing/Data/Numerical/MATLAB/planarVelocityPOD/', caseName, '/', dataID], 'dir')
            mkdir(['/mnt/Processing/Data/Numerical/MATLAB/planarVelocityPOD/', caseName, '/', dataID]);
        end
        
        save(['/mnt/Processing/Data/Numerical/MATLAB/planarVelocityPOD/', caseName, '/', dataID, '/', planeName], ...
              'PODdata', 'sampleInterval', 'normalise', '-v7.3', '-noCompression');
        disp(['    Saving to: ~/Data/Numerical/MATLAB/planarVelocityPOD/', caseName, '/', dataID, '/', planeName]);
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
        return
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
reconData.u.mean = PODdata.u.mean;
reconData.v.mean = PODdata.v.mean;
reconData.w.mean = PODdata.w.mean;
reconData.u.inst = cell(Nt,1);
reconData.v.inst = reconData.u.inst;
reconData.w.inst = reconData.u.inst;

for i = 1:Nt
    reconData.u.inst{i} = reconData.u.mean;
    reconData.v.inst{i} = reconData.v.mean;
    reconData.w.inst{i} = reconData.w.mean;
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








%% Local Functions

function modes = inputModes(nModes)

    modes = str2num(input('    Input Desired Modes (Row Vector Form) [s]: ', 's')); %#ok<ST2NM>
    
    if isempty(modes) || any(isnan(modes)) || ~isrow(modes) > 1 || any(modes <= 0) || any(modes >= nModes)
        disp('        WARNING: Invalid Entry');
        modes = -1;
    end

end