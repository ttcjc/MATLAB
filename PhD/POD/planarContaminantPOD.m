%% Planar Lagrangian Contaminant POD Calculator v2.0

clear variables;
close all;
clc;
evalc('delete(gcp(''nocreate''));');

nProc = maxNumCompThreads - 2; % Number of Processors Used for Parallel Collation

fig = 0; % Initialise Figure Tracking
figHold = 0; % Enable Overwriting of Figures

disp('======================================');
disp('Planar Contaminant POD Calculator v2.0');
disp('======================================');

disp(' ');
disp(' ');


%% Changelog

% v1.0 - Initial Commit (Base Contamination and Far-Field Extraction Plane)
% v1.1 - Rename and Restructure to Account for Changes to 'mapData' Output
% v2.0 - Rewrite, Accommodating New OpenFOAM Data Formats


%% Select Mapping Location

disp('Mapping Location');
disp('-----------------');

disp(' ');

disp('Possible Mapping Locations:');
disp('    A: Base Contamination');
disp('    B: Far-Field Spray Transport');

valid = false;
while ~valid
    disp(' ');
    selection = input('Select Mapping Location [A/B]: ', 's');

    if selection == 'a' | selection == 'A' %#ok<OR2>
        format = 'A';
        
        valid = true;
    elseif selection == 'b' | selection == 'B' %#ok<OR2>
        format = 'B';
        
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;

disp(' ');
disp(' ');


%% Acquire Contaminant Map

disp('Contaminant Map Acquisition');
disp('----------------------------');

valid = false;
while ~valid
    disp(' ');
    [fileName, filePath] = uigetfile('/mnt/Processing/Data/Numerical/MATLAB/contaminantMap/*.mat', ...
                                      'Select Map Data');
    
    switch format
        
        case 'A'
            
            if contains(filePath, '/base')
                disp(['Loading ''', fileName, '''...']);
                dataID = load([filePath, fileName], 'dataID').dataID;
                mapData = load([filePath, fileName], 'mapData').mapData;
                sampleInterval = load([filePath, fileName], 'sampleInterval').sampleInterval;
                dLims = load([filePath, fileName], 'dLims').dLims;
                normalise = load([filePath, fileName], 'normalise').normalise;
                disp('    Success');
                
                valid = true;
            else
                disp('WARNING: Invalid File Selection');
                clear fileName filePath;
            end
            
        case 'B'
            
            if contains(filePath, '/X_')
                disp(['Loading ''', fileName, '''...']);
                dataID = load([filePath, fileName], 'dataID').dataID;
                mapData = load([filePath, fileName], 'mapData').mapData;
                sampleInterval = load([filePath, fileName], 'sampleInterval').sampleInterval;
                dLims = load([filePath, fileName], 'dLims').dLims;
                normalise = load([filePath, fileName], 'normalise').normalise;
                disp('    Success');
                
                valid = true;
            else
                disp('WARNING: Invalid File Selection');
                clear fileName filePath;
            end
            
    end
end
clear valid;

namePos = strfind(filePath, '/');
caseName = filePath((namePos(end - 2) + 1):(namePos(end - 1) - 1));

timePrecision = strfind(fileName, '_T') - 3;

switch format
    
    case 'B'
        planePos = filePath((namePos(end - 1) + 1):(namePos(end) - 1));

end

disp(' ');
disp(' ');


%% Select Relevant Geometry and Define Bounding Box

[geometry, xDims, yDims, zDims, spacePrecision, normalise] = selectGeometry(normalise);

disp(' ');
disp(' ');


%% Select POD Options

disp('POD Options');
disp('------------');

% Select Variable of Interest
PODvar = fieldnames(mapData.mean);
PODvar = PODvar(1:(end - 1));

valid = false;
while ~valid
    disp(' ');
    [index, valid] = listdlg('listSize', [300, 300], ...
                             'selectionMode', 'single', ...
                             'name', 'Select Variable for Decomposition', ...
                             'listString', PODvar);

    if ~valid
        disp('WARNING: No Mapping Variable Selected');
    end
end
clear valid;

PODvar = PODvar{index};

disp(['Variable of Interest: ', PODvar]);

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
switch format
    
    case 'A'
        % Identify Model Base
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
            mapPerim = [mapPerim; mapPerim(1,:)]; % Close Boundary
        end
        
        clear basePoints basePoly;
        
        xLimsData = xDims(2);
        yLimsData = [min(mapPerim(:,2)); max(mapPerim(:,2))];
        zLimsData = [min(mapPerim(:,3)); max(mapPerim(:,3))];
        
    case 'B'
    
        if contains(caseName, ["Run_Test", "Windsor"])
            mapPerim = [];
            
            xLimsData = mapData.positionGrid(1,1);
            yLimsData = [-0.5945; 0.5945];
            zLimsData = [0; 0.739];
            
            if normalise
                yLimsData = round((yLimsData / 1.044), spacePrecision);
                zLimsData = round((zLimsData / 1.044), spacePrecision);
            end
            
        end
        
end

% Initialise POD Variables
PODdata.positionGrid = mapData.positionGrid;
PODdata.time = mapData.inst.time;
PODdata.(PODvar).mean = mapData.mean.(PODvar);
PODdata.(PODvar).inst = mapData.inst.(PODvar);

clear mapData;

% Initialise Progress Bar
wB = waitbar(0, ['Calculating Instantaneous ''', PODvar, ''' Fluctuations'], 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';

% Calculate Instantaneous Variable Fluctuations
PODdata.(PODvar).prime = cell(height(PODdata.time),1);

for i = 1:height(PODdata.time)
    PODdata.(PODvar).prime{i} = PODdata.(PODvar).inst{i} - PODdata.(PODvar).mean;
    
    waitbar((i / height(PODdata.time)), wB);
end

delete(wB);

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

varPrime = PODdata.(PODvar).prime;
parfor i = 1:Nt
    
    for j = 1:Ns
        snapshotMatrix(i,j) = varPrime{i}(j);
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

switch format
    
    case 'A'
        figName = ['Base_', PODvar, '_Planar_POD_Energy_Content'];
        
    case 'B'
        figName = [planePos, '_', PODvar, '_Planar_POD_Energy_Content'];
        
end
        
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
xlabel({' ', '{\bf{Mode}}'}, 'fontName', 'LM Roman 12');
ylabel({'{\bf{Energy Content (\it{%})}}', ' '}, 'fontName', 'LM Roman 12');
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
    switch format

        case 'A'
            orientation = 'YZ';

            if contains(caseName, ["Run_Test", "Windsor"])
                xLimsPlot = [0.31875; 1.52725];
                yLimsPlot = [-0.2445; 0.2445];
                zLimsPlot = [0; 0.389];
            end

        case 'B'
            orientation = 'YZ';

            if contains(caseName, ["Run_Test", "Windsor"])
                xLimsPlot = [0.31875; 4.65925];
                yLimsPlot = [-0.5945; 0.5945];
                zLimsPlot = [0; 0.739];
            end

    end

    if normalise
        xLimsPlot = round((xLimsPlot / 1.044), spacePrecision);
        yLimsPlot = round((yLimsPlot / 1.044), spacePrecision);
        zLimsPlot = round((zLimsPlot / 1.044), spacePrecision);
    end

    positionData = PODdata.positionGrid;
    cMap = turbo(24);
    figTitle = '-'; % Leave Blank ('-') for Formatting Purposes
    cLims = [-1; 1];

    for i = plotModes
        disp(['    Presenting Mode #', num2str(i), '...']);

        scalarData = rescale(PODdata.phi_mode(:,i), -1, 1);

        switch format

            case 'A'
                figName = ['Base_', PODvar, '_Planar_POD_M', num2str(i)];

            case 'B'
                figName = [planePos, '_', PODvar, '_Planar_POD_M', num2str(i)];

        end

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
        
        switch format
            
            case 'A'
                
                if ~exist(['/mnt/Processing/Data/Numerical/MATLAB/planarContaminantPOD/', caseName, '/base/', PODvar], 'dir')
                    mkdir(['/mnt/Processing/Data/Numerical/MATLAB/planarContaminantPOD/', caseName, '/base/', PODvar]);
                end
                
            case 'B'
                
                if ~exist(['/mnt/Processing/Data/Numerical/MATLAB/planarContaminantPOD/', caseName, '/', planePos, '/', PODvar], 'dir')
                    mkdir(['/mnt/Processing/Data/Numerical/MATLAB/planarContaminantPOD/', caseName, '/', planePos, '/', PODvar]);
                end
                
        end
        
        switch format
            
            case 'A'
                disp(['    Saving to: /mnt/Processing/Data/Numerical/MATLAB/planarContaminantPOD/', caseName, '/base/', PODvar, '/', dataID, '.mat']);
                save(['/mnt/Processing/Data/Numerical/MATLAB/planarContaminantPOD/', caseName, '/base/', PODvar, '/', dataID, '.mat'], ...
                     'dataID', 'PODdata', 'sampleInterval', 'dLims', 'normalise', '-v7.3', '-noCompression');
                disp('        Success');
                 
            case 'B'
                disp(['    Saving to: /mnt/Processing/Data/Numerical/MATLAB/planarContaminantPOD/', caseName, '/', planePos, '/', PODvar, '/', dataID, '.mat']);
                save(['/mnt/Processing/Data/Numerical/MATLAB/planarContaminantPOD/', caseName, '/', planePos, '/', PODvar, '/', dataID, '.mat'], ...
                     'dataID', 'PODdata', 'sampleInterval', 'dLims', 'normalise', '-v7.3', '-noCompression');
                disp('        Success');
        
        end
        
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
reconData.(PODvar).mean = PODdata.(PODvar).mean;
reconData.(PODvar).inst = cell(Nt,1);

for i = 1:Nt
    reconData.(PODvar).inst{i} = reconData.(PODvar).mean;
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
    mode = ['M', num2str(i)];
    
    modeMatrix = PODdata.A_coeff(:,i) * PODdata.phi_mode(:,i)';
    varPrime = cell(Nt,1);
    
    parfor j = 1:Nt
        varPrime{j} = zeros(Ns,1);
        
        for k = 1:Ns
            varPrime{j}(k) = modeMatrix(j,k);
        end
        
        send(dQ, []);
    end
    
    delete(wB);
    
    reconData.(mode).modeMatrix = modeMatrix;
    reconData.(mode).prime = varPrime;
    clear modeMatrix varPrime;
    
    % Add Mode to Reconstruction
    for j = 1:Nt
        reconData.(PODvar).inst{j} = reconData.(PODvar).inst{j} + reconData.(mode).prime{j};
    end
    
end

if any(strcmp(PODvar, {'mass', 'massNorm'}))
    % Initialise Progress Bar
    wB = waitbar(0, 'Calculating Reconstructed Centre of Mass', 'name', 'Progress');
    wB.Children.Title.Interpreter = 'none';
    dQ = parallel.pool.DataQueue;
    afterEach(dQ, @parforWaitBar);
    
    parforWaitBar(wB, Nt);
    
    % Calculate Reconstructed CoM
    CoM = cell(Nt,1);
    
    switch format
        
        case {'A', 'B'}
            mass = reconData.(PODvar).inst;
            positionGrid = reconData.positionGrid;
            parfor i = 1:Nt
                CoM{i} = zeros(1,3);
                CoM{i}(1) = positionGrid(1,1); %#ok<PFBNS>
                
                for j = 1:height(positionGrid)
                    CoM{i}(2) = CoM{i}(2) + (mass{i}(j) * positionGrid(j,2));
                    CoM{i}(3) = CoM{i}(3) + (mass{i}(j) * positionGrid(j,3));
                end
                
                CoM{i}(2) = CoM{i}(2) / sum(mass{i});
                CoM{i}(3) = CoM{i}(3) / sum(mass{i});
                
                send(dQ, []);
            end
            clear mass positionGrid;
            
    end
    
    delete(wB);
    
    reconData.(PODvar).CoM = CoM;
    clear CoM;

end

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
    switch format

        case 'A'
            orientation = 'YZ';

            if contains(caseName, ["Run_Test", "Windsor"])
                xLimsPlot = [0.31875; 1.52725];
                yLimsPlot = [-0.2445; 0.2445];
                zLimsPlot = [0; 0.389];
            end

        case 'B'
            orientation = 'YZ';

            if contains(caseName, ["Run_Test", "Windsor"])
                xLimsPlot = [0.31875; 4.65925];
                yLimsPlot = [-0.5945; 0.5945];
                zLimsPlot = [0; 0.739];
            end

    end

    if normalise
        xLimsPlot = round((xLimsPlot / 1.044), spacePrecision);
        yLimsPlot = round((yLimsPlot / 1.044), spacePrecision);
        zLimsPlot = round((zLimsPlot / 1.044), spacePrecision);
    end

    positionData = reconData.positionGrid;
    cMap = viridis(24);
    figTitle = '-'; % Leave Blank ('-') for Formatting Purposes
    
    if any(strcmp(PODvar, {'d10', 'd20', 'd30', 'd32'}))
        cLims = dLims;
    elseif strcmp(PODvar, 'massNorm')
        cLims = [0; 1];
    else
        cLims = [0; max(cellfun(@max, reconData.(PODvar).inst))];
    end

    figHold = fig;

    for i = 1:Nt
        
        if i ~= 1
            clf(fig);
            fig = figHold;
        end
        
        scalarData = reconData.(PODvar).inst{i};
        figTime = num2str(reconData.time(i), ['%.', num2str(timePrecision), 'f']);
        
        switch format
            
            case 'A'
                figName = ['Base_', PODvar, '_Reconstruction_T', erase(figTime, '.')];
            
            case 'B'
                figName = [planePos, '_', PODvar, '_Reconstruction_T', erase(figTime, '.')];
        
        end
        
        if any(strcmp(PODvar, {'mass', 'massNorm'}))
            CoM = reconData.(PODvar).CoM{i};
        else
            CoM = [];
        end
        
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
        
        switch format
            
            case 'A'
                
                if ~exist(['/mnt/Processing/Data/Numerical/MATLAB/planarContaminantReconstruction/', caseName, '/base/', PODvar], 'dir')
                    mkdir(['/mnt/Processing/Data/Numerical/MATLAB/planarContaminantReconstruction/', caseName, '/base/', PODvar]);
                end
                
            case 'B'
                
                if ~exist(['/mnt/Processing/Data/Numerical/MATLAB/planarContaminantReconstruction/', caseName, '/', planePos, '/', PODvar], 'dir')
                    mkdir(['/mnt/Processing/Data/Numerical/MATLAB/planarContaminantReconstruction/', caseName, '/', planePos, '/', PODvar]);
                end
                
        end
        
        switch format
            
            case 'A'
                disp(['    Saving to: /mnt/Processing/Data/Numerical/MATLAB/planarContaminantReconstruction/', caseName, '/base/', PODvar, '/', dataID, '.mat']);
                save(['/mnt/Processing/Data/Numerical/MATLAB/planarContaminantReconstruction/', caseName, '/base/', PODvar, '/', dataID, '.mat'], ...
                     'dataID', 'reconData', 'sampleInterval', 'dLims', 'normalise', '-v7.3', '-noCompression');
                disp('        Success');
                 
            case 'B'
                disp(['    Saving to: /mnt/Processing/Data/Numerical/MATLAB/planarContaminantReconstruction/', caseName, '/', planePos, '/', PODvar, '/', dataID, '.mat']);
                save(['/mnt/Processing/Data/Numerical/MATLAB/planarContaminantReconstruction/', caseName, '/', planePos, '/', PODvar, '/', dataID, '.mat'], ...
                     'dataID', 'reconData', 'sampleInterval', 'dLims', 'normalise', '-v7.3', '-noCompression');
                disp('        Success');
        
        end
        
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid


%% Local Functions

function modes = inputModes(nModes)

    modes = str2num(input('    Input Desired Modes (Row Vector Form) [s]: ', 's')); %#ok<ST2NM>
    
    if isempty(modes) || any(isnan(modes)) || ~isrow(modes) > 1 || any(modes <= 0) || any(modes >= nModes)
        disp('        WARNING: Invalid Entry');
        modes = -1;
    end

end
