%% Preamble

clear variables;
close all;
clc;
evalc('delete(gcp(''nocreate''));');

fig = 0; % Initialise Figure Tracking
figHold = 0; % Enable Overwriting of Figures

saveLocation = '/mnt/Processing/Data';


%% PIV Data Processor v1.0

normalise = false; % Normalisation of Dimensions
nProc = maxNumCompThreads - 2; % Number of Processors Used for Parallel Collation

disp('=======================');
disp('PIV Data Processor v1.0');
disp('=======================');

disp(' ');
disp(' ');


%% Changelog

% v1.0 - Initial Commit


%% Initialise Case

% Select Relevant Geometry and Define Bounding Box
[geometry, xDims, yDims, zDims, spacePrecision, normalise] = selectGeometry(normalise);

disp(' ');
disp(' ');

% Load PIV Data
[caseFolder, caseName, PIVdata, format, samplingFrequency] = initialisePIVdata(saveLocation, nProc);


%% Perform PIV Processing

disp('PIV Processing');
disp('---------------');

disp(' ');

disp('***********');
disp('  RUNNING ');

tic;

disp(' ');

disp('    Initialising...');

switch format
    
    case 'planar'
        
    case 'stereo'
        
    case 'tomo'
        % Adjust Origin
        if contains(caseName, 'Cyclist_H_50')
            PIVdata.positionGrid(:,1) = PIVdata.positionGrid(:,1) + (0.04 + 0.235 + 0.1545);
            PIVdata.positionGrid(:,3) = PIVdata.positionGrid(:,3) + (0.622 - 0.1545);
        end
    
end

disp(' ');

% Remove Invalid Vectors
disp('    Removing Invalid Vectors...');

switch format
    
    case 'planar'
        
    case 'stereo'
        
    case 'tomo'

        for i = 1:height(PIVdata.time)
            index = find(PIVdata.isValid{i} == 0);

            if i == 1
                invalidPos = index;
            else
                invalidPos = union(invalidPos, index);
            end

        end
        clear i;

        PIVdata.positionGrid(index,:) = [];

        for i = 1:height(PIVdata.time)
            PIVdata.u{i}(index) = [];
            PIVdata.v{i}(index) = [];
            PIVdata.w{i}(index) = [];
            PIVdata.corrCoeff{i}(index) = [];
        end
        clear i;

        PIVdata = rmfield(PIVdata, 'isValid');
        
end

disp(' ');

% Remove Frames With Excessive Interpolation
disp('    Removing Poorly Seeded Frames...');

switch format
    
    case 'planar'
        
    case 'stereo'
        
    case 'tomo'
        i = 1;
        while i <= height(PIVdata.time)

            if (height(PIVdata.corrCoeff{i}(PIVdata.corrCoeff{i} < 0.1)) / height(PIVdata.positionGrid)) > 0.1
                PIVdata.time(i) = [];
                PIVdata.u(i) = [];
                PIVdata.v(i) = [];
                PIVdata.w(i) = [];
                PIVdata.corrCoeff(i) = [];
            else
                i = i + 1;
            end

        end
        clear i;  
end

disp(' ');

% Normalise Velocity
disp('    Normalising Velocities...');

switch format
    
    case 'planar'
        
    case 'stereo'
        
    case 'tomo'
        U = str2double(cell2mat(extractBetween(caseName, 'U_', '_Yaw')));

        for i = 1:height(PIVdata.time)
            PIVdata.u{i} = PIVdata.u{i} / U;
            PIVdata.v{i} = PIVdata.v{i} / U;
            PIVdata.w{i} = PIVdata.w{i} / U;
        end
        clear i;
        
end

disp(' ');

% Calculate Ensemble Averages
disp('    Calculating Ensamble Averages...');

switch format
    
    case 'planar'
        
    case 'stereo'
        
    case 'tomo'
        PIVdata.uMean = zeros(height(PIVdata.positionGrid),1);
        PIVdata.vMean = zeros(height(PIVdata.positionGrid),1);
        PIVdata.wMean = zeros(height(PIVdata.positionGrid),1);
        PIVdata.corrCoeffMean = zeros(height(PIVdata.positionGrid),1);

        for i = 1:height(PIVdata.time)
            PIVdata.uMean = PIVdata.uMean + PIVdata.u{i};
            PIVdata.vMean = PIVdata.vMean + PIVdata.v{i};
            PIVdata.wMean = PIVdata.wMean + PIVdata.w{i};
            PIVdata.corrCoeffMean = PIVdata.corrCoeffMean + PIVdata.corrCoeff{i};
        end
        clear i;

        PIVdata.uMean = PIVdata.uMean / height(PIVdata.time);
        PIVdata.vMean = PIVdata.vMean / height(PIVdata.time);
        PIVdata.wMean = PIVdata.wMean / height(PIVdata.time);
        PIVdata.corrCoeffMean = PIVdata.corrCoeffMean / height(PIVdata.time);
        
end

executionTime = toc;

disp(' ');

disp(['    Run Time: ', num2str(executionTime), 's']);

disp(' ');

disp('  SUCCESS  ');
disp('***********');

disp(' ');
disp(' ');


%% Plotting

% Extract Centreline Plane
if any(PIVdata.positionGrid(:,2) == min(abs(unique(PIVdata.positionGrid(:,2)))))
    index = (PIVdata.positionGrid(:,2) == min(abs(unique(PIVdata.positionGrid(:,2)))));
else
    index = (PIVdata.positionGrid(:,2) == -min(abs(unique(PIVdata.positionGrid(:,2)))));
end

positionData = PIVdata.positionGrid(index,:);
vectorData = [PIVdata.uMean(index), PIVdata.vMean(index), PIVdata.wMean(index)];
% vectorData = [PIVdata.u{1}(index), PIVdata.v{1}(index), PIVdata.w{1}(index)];

orientation = 'XZ';
xLimsData = [min(positionData(:,1)), max(positionData(:,1))];
yLimsData = positionData(1,2);
zLimsData = [min(positionData(:,3)), max(positionData(:,3))];
nComponents = 3;
component = [];
figName = 'Tomo_Test';
cMap = parula(32);
streamlines = true;
figTitle = '-'; % Leave Blank ('-') for Formatting Purposes
figSubtitle = ' ';
cLims = [0.2; 1];
xLimsPlot = [-0.6; 0.6];
yLimsPlot = [-0.2; 0.2];
zLimsPlot = [0; 0.8];

fig = planarVectorPlots(orientation, xLimsData, yLimsData, zLimsData, positionData, ...
                        vectorData, nComponents, component, fig, figName, cMap, geometry, ...
                        streamlines, xDims, yDims, zDims, figTitle, figSubtitle, cLims, ...
                        xLimsPlot, yLimsPlot, zLimsPlot, normalise);


