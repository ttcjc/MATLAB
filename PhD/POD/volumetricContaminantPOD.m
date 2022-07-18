%% Volumetric Lagrangian Contaminant POD Calculator v2.0

clear variables;
close all;
clc;
evalc('delete(gcp(''nocreate''));');

saveLocation = '/mnt/Processing/Data';
% saveLocation = '~/Data';

fig = 0; % Initialise Figure Tracking
figHold = 0; % Enable Overwriting of Figures

disp('=====================================================');
disp('Volumetric Lagrangian Contaminant POD Calculator v2.0');
disp('=====================================================');

disp(' ');
disp(' ');


%% Changelog

% v1.0 - Initial Commit
% v2.0 - Rewrite, Accommodating New OpenFOAM Data Formats


%% Select Region of Interest

disp('Region of Interest');
disp('-------------------');

disp(' ');

disp('Possible Regions of Interest:');
disp('    A: Near-Field');
disp('    B: Far-Field');

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
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;

disp(' ');
disp(' ');


%% Acquire Volume Field

disp('Volume Field Acquisition');
disp('-------------------------');

valid = false;
while ~valid
    disp(' ');
    [fileName, filePath] = uigetfile([saveLocation, '/Numerical/MATLAB/volumeField/*.mat'], ...
                                      'Select Map Data');
    
    switch format
        
        case 'A'
            
            if contains(filePath, '/nearField')
                disp(['Loading ''', fileName, '''...']);
                dataID = load([filePath, fileName], 'dataID').dataID;
                volumeData = load([filePath, fileName], 'volumeData').volumeData;
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
            
            if contains(filePath, '/farField')
                disp(['Loading ''', fileName, '''...']);
                dataID = load([filePath, fileName], 'dataID').dataID;
                volumeData = load([filePath, fileName], 'volumeData').volumeData;
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
PODvar = fieldnames(volumeData.mean);

valid = false;
while ~valid
    disp(' ');
    [index, valid] = listdlg('listSize', [300, 300], ...
                             'selectionMode', 'single', ...
                             'name', 'Select Variable for Decomposition', ...
                             'listString', PODvar);

    if ~valid
        disp('WARNING: No Mapping Variable Selected');
        continue
    end
    
    if ~strcmp(PODvar{index}, 'volFraction')
        disp('WARNING: Unsupported POD Variable');
        valid = false;
    end
        
end
clear valid;

PODvar = PODvar{index};

disp(['Variable of Interest: ', PODvar]);

disp(' ');
disp(' ');


%% Volumetric Planar POD (Snapshot Method)

disp('Volumetric Proper Orthogonal Decomposition');
disp('-------------------------------------------');

disp(' ');

disp('***********');
disp('  RUNNING ');

tic;

disp(' ');

disp('    Initialising...');

% Initialise POD Variables
PODdata.positionGrid = volumeData.positionGrid;
PODdata.time = volumeData.inst.time;
PODdata.(PODvar).mean = volumeData.mean.(PODvar);
PODdata.(PODvar).inst = volumeData.inst.(PODvar);

clear volumeData;

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

disp('    Performing Volumetric POD Using the Snapshot Method...');

Ns = height(PODdata.positionGrid); % Number of Spatial Points
Nt = height(PODdata.time); % Number of Time Instances

% Initialise Progress Bar
wB = waitbar(0, 'Assembling Snapshot Matrix', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';

% Assemble Snapshot Matrix
PODdata.snapshotMatrix = zeros(Nt,Ns);

varPrime = PODdata.(PODvar).prime;
for i = 1:Nt
    
    for j = 1:Ns
        PODdata.snapshotMatrix(i,j) = PODdata.(PODvar).prime{i}(j);
    end
    
    waitbar((i / Nt), wB);
end

delete(wB);

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
        figName = ['Near_Field_', PODvar, '_Planar_POD_Energy_Content'];
        
    case 'B'
        figName = ['Far_Field_', PODvar, '_Planar_POD_Energy_Content'];
        
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
            continue
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

if plotModes
    % Specify Region Boundaries
    switch format

        case 'A'

            if contains(caseName, ["Run_Test", "Windsor"])
                xLimsData = [0.31875; 1.38325];
                yLimsData = [-0.3445; 0.3445];
                zLimsData = [0; 0.489];
            end

        case 'B'

            if contains(caseName, ["Run_Test", "Windsor"])
                xLimsData = [0.31875; 3.61525];
                yLimsData = [-0.5945; 0.5945];
                zLimsData = [0; 0.739];
            end

    end

    if normalise
        xLimsData = round((xLimsData / 1.044), spacePrecision);
        yLimsData = round((yLimsData / 1.044), spacePrecision);
        zLimsData = round((zLimsData / 1.044), spacePrecision);
    end
    
    % Format Position Grid
    gridShape = [height(unique(PODdata.positionGrid(:,1))), ...
                 height(unique(PODdata.positionGrid(:,2))), ...
                 height(unique(PODdata.positionGrid(:,3)))];
             
    xInit = reshape(PODdata.positionGrid(:,1), gridShape);
    yInit = reshape(PODdata.positionGrid(:,2), gridShape);
    zInit = reshape(PODdata.positionGrid(:,3), gridShape);
    POD = true;
    isoValue = 5e-7;
    cMap = cool2warm(24);
    figTitle = '-'; % Leave Blank ('-') for Formatting Purposes
    xLimsPlot = xLimsData;
    yLimsPlot = yLimsData;
    zLimsPlot = zLimsData;
    
    for i = nModes
        disp(['    Presenting Mode #', num2str(i), '...']);
        
        modeMatrix = PODdata.A_coeff(:,i) * PODdata.phi_mode(:,i)';
        
        [~, indexA] = min(PODdata.A_coeff(:,i));
        [~, indexB] = max(PODdata.A_coeff(:,i));
        
        primeFieldA = zeros(Ns,1);
        primeFieldB = zeros(Ns,1);
        
        for j = 1:Ns
            primeFieldA(j) = modeMatrix(indexA,j);
            primeFieldB(j) = modeMatrix(indexB,j);
        end
        
        fieldDataA = reshape(primeFieldA, gridShape);
        fieldDataB = reshape(primeFieldB, gridShape);
        
        fieldData = {fieldDataA; fieldDataB};
        
        clear indexA indexB primeFieldA primeFieldB fieldDataA fieldDataB
        
        switch format

            case 'A'
                figName = ['Near_Field_', PODvar, '_Volumetric_POD_M', num2str(i)];

            case 'B'
                figName = ['Far_Field_', PODvar, '_Volumetric_POD_M', num2str(i)];

        end

        figSubtitle = [num2str(round(PODdata.modeEnergy(i), 2), '%.2f'), '\it{%}'];

        fig = volumeFieldPlots(xLimsData, yLimsData, zLimsData, xInit, yInit, zInit, POD, fieldData, ...
                               fig, figName, geometry, isoValue, cMap, figTitle, figSubtitle, ...
                               xLimsPlot, yLimsPlot, zLimsPlot);
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
                
                if ~exist([saveLocation, '/Numerical/MATLAB/volumetricContaminantPOD/', caseName, '/nearField/', PODvar], 'dir')
                    mkdir([saveLocation, '/Numerical/MATLAB/volumetricContaminantPOD/', caseName, '/nearField/', PODvar]);
                end
                
            case 'B'
                
                if ~exist([saveLocation, '/Numerical/MATLAB/volumetricContaminantPOD/', caseName, '/farField/', PODvar], 'dir')
                    mkdir([saveLocation, '/Numerical/MATLAB/volumetricContaminantPOD/', caseName, '/farField/', PODvar]);
                end
                
        end
        
        switch format
            
            case 'A'
                disp(['    Saving to: ', saveLocation, '/Numerical/MATLAB/volumetricContaminantPOD/', caseName, '/nearField/', PODvar, '/', dataID, '.mat']);
                save([saveLocation, '/Numerical/MATLAB/volumetricContaminantPOD/', caseName, '/nearField/', PODvar, '/', dataID, '.mat'], ...
                     'dataID', 'PODdata', 'sampleInterval', 'dLims', 'normalise', '-v7.3', '-noCompression');
                disp('        Success');
                 
            case 'B'
                disp(['    Saving to: ', saveLocation, '/Numerical/MATLAB/volumetricContaminantPOD/', caseName, '/farField/', PODvar, '/', dataID, '.mat']);
                save([saveLocation, '/Numerical/MATLAB/volumetricContaminantPOD/', caseName, '/farField/', PODvar, '/', dataID, '.mat'], ...
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
        reconstruct = false;
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        reconstruct = true;
        nModes = inputModes(Nt);
        
        if nModes == -1
            continue
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
reconData.(PODvar).mean = PODdata.(PODvar).mean;
reconData.(PODvar).inst = cell(Nt,1);

for i = 1:Nt
    reconData.(PODvar).inst{i} = reconData.(PODvar).mean;
end

disp(' ');

% Perform Reconstruction
disp('    Performing Field Reconstruction...');

for i = nModes
    % Initialise Progress Bar
    wB = waitbar(0, ['Adding Mode #', num2str(i), ' to Reconstruction'], 'name', 'Progress');
    wB.Children.Title.Interpreter = 'none';
    
    % Identify Mode Contribution
    mode = ['M', num2str(i)];
    
    reconData.(mode).modeMatrix = PODdata.A_coeff(:,i) * PODdata.phi_mode(:,i)';
    reconData.(mode).prime = cell(Nt,1);
    
    for j = 1:Nt
        reconData.(mode).prime{j} = zeros(Ns,1);
        
        for k = 1:Ns
            reconData.(mode).prime{j}(k) = reconData.(mode).modeMatrix(j,k);
        end
        
        waitbar((j / Nt), wB);
    end
    
    delete(wB);
    
    % Add Mode to Reconstruction
    for j = 1:Nt
        reconData.(PODvar).inst{j} = reconData.(PODvar).inst{j} + reconData.(mode).prime{j};
    end
    
end

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
            continue
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
    % Specify Region Boundaries
    switch format

        case 'A'

            if contains(caseName, ["Run_Test", "Windsor"])
                xLimsData = [0.31875; 1.38325];
                yLimsData = [-0.3445; 0.3445];
                zLimsData = [0; 0.489];
            end

        case 'B'

            if contains(caseName, ["Run_Test", "Windsor"])
                xLimsData = [0.31875; 3.61525];
                yLimsData = [-0.5945; 0.5945];
                zLimsData = [0; 0.739];
            end

    end

    if normalise
        xLimsData = round((xLimsData / 1.044), spacePrecision);
        yLimsData = round((yLimsData / 1.044), spacePrecision);
        zLimsData = round((zLimsData / 1.044), spacePrecision);
    end
    
    % Format Position Grid
    gridShape = [height(unique(reconData.positionGrid(:,1))), ...
                 height(unique(reconData.positionGrid(:,2))), ...
                 height(unique(reconData.positionGrid(:,3)))];
             
    xInit = reshape(reconData.positionGrid(:,1), gridShape);
    yInit = reshape(reconData.positionGrid(:,2), gridShape);
    zInit = reshape(reconData.positionGrid(:,3), gridShape);
    POD = false;
    isoValue = 1e-6;

    if strcmp(caseName, 'Windsor_SB_wW_Upstream_SC')
        cMap = ([74, 24, 99] / 255);
    elseif strcmp(caseName, 'Windsor_ST_20D_wW_Upstream_SC')
        cMap = ([230, 0, 126] / 255);
    elseif strcmp(caseName, 'Windsor_RSST_16D_U50_wW_Upstream_SC')
        cMap = ([34, 196, 172] / 255);
    else
        cMap = ([252, 194, 29] / 255);
    end

    figTitle = '-'; % Leave Blank ('-') for Formatting Purposes
    xLimsPlot = xLimsData;
    yLimsPlot = yLimsData;
    zLimsPlot = zLimsData;
    
    figHold = fig;
    
    for i = 1:nFrames
        
        if i ~= 1
            clf(fig);
            fig = figHold;
        end

        fieldData = reshape(reconData.(PODvar).inst{i}, gridShape);
        figTime = num2str(reconData.time(i), ['%.', num2str(timePrecision), 'f']);

        switch format

            case 'A'
                figName = ['Near_Field_', PODvar, '_Reconstruction_T', erase(figTime, '.')];

            case 'B'
                figName = ['Far_Field_', PODvar, '_Reconstruction_T', erase(figTime, '.')];

        end
        
        figSubtitle = [figTime, ' \it{s}'];
        
        fig = volumeFieldPlots(xLimsData, yLimsData, zLimsData, xInit, yInit, zInit, POD, fieldData, ...
                               fig, figName, geometry, isoValue, cMap, figTitle, figSubtitle, ...
                               xLimsPlot, yLimsPlot, zLimsPlot);

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
                
                if ~exist([saveLocation, '/Numerical/MATLAB/volumetricContaminantReconstruction/', caseName, '/nearField/', PODvar], 'dir')
                    mkdir([saveLocation, '/Numerical/MATLAB/volumetricContaminantReconstruction/', caseName, '/nearField/', PODvar]);
                end
                
            case 'B'
                
                if ~exist([saveLocation, '/Numerical/MATLAB/volumetricContaminantReconstruction/', caseName, '/farField/', PODvar], 'dir')
                    mkdir([saveLocation, '/Numerical/MATLAB/volumetricContaminantReconstruction/', caseName, '/farField/', PODvar]);
                end
                
        end
        
        switch format
            
            case 'A'
                disp(['    Saving to: ', saveLocation, '/Numerical/MATLAB/volumetricContaminantReconstruction/', caseName, '/nearField/', PODvar, '/', dataID, '.mat']);
                save([saveLocation, '/Numerical/MATLAB/volumetricContaminantReconstruction/', caseName, '/nearField/', PODvar, '/', dataID, '.mat'], ...
                     'dataID', 'reconData', 'sampleInterval', 'dLims', 'normalise', '-v7.3', '-noCompression');
                disp('        Success');
                 
            case 'B'
                disp(['    Saving to: ', saveLocation, '/Numerical/MATLAB/volumetricContaminantReconstruction/', caseName, '/farField/', PODvar, '/', dataID, '.mat']);
                save([saveLocation, '/Numerical/MATLAB/volumetricContaminantReconstruction/', caseName, '/farField/', PODvar, '/', dataID, '.mat'], ...
                     'dataID', 'reconData', 'sampleInterval', 'dLims', 'normalise', '-v7.3', '-noCompression');
                disp('        Success');
        
        end
        
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

    nFrames = str2double(input(['    Input Desired Frame Count [1 - ', num2str(Nt), ']: '], 's'));
    
    if isnan(nFrames) || nFrames <= 0 || nFrames > Nt
        disp('        WARNING: Invalid Entry');
        nFrames = -1;
    end

end