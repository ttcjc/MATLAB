%% Lagrangian Line of Sight Calculator v1.0

clear variables;
close all;
clc;
evalc('delete(gcp(''nocreate''));');

% saveLocation = '/mnt/Processing/Data';
saveLocation = '~/Data';

nProc = maxNumCompThreads - 2; % Number of Processors Used for Parallel Collation

fig = 0; % Initialise Figure Tracking
figHold = 0; % Enable Overwriting of Figures

disp('==============================');
disp('Depth of Field Calculator v1.0');
disp('==============================');

disp(' ');
disp(' ');


%% Changelog

% v1.0 - Initial Commit


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
                                      'Select Volumetric Data');
    
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


%% Select Target Point

disp('Target Point Definition');
disp('------------------------');

disp(' ');

% Restore Original Dimensions
if normalise

    if contains(caseName, ["Run_Test", "Windsor"])
        volumeData.positionGrid = round((volumeData.positionGrid * 1.044), spacePrecision);
    end

end

% Identify Data Limits
xLimsData = [min(volumeData.positionGrid(:,1)); max(volumeData.positionGrid(:,1))];
yLimsData = [min(volumeData.positionGrid(:,2)); max(volumeData.positionGrid(:,2))];
zLimsData = [min(volumeData.positionGrid(:,3)); max(volumeData.positionGrid(:,3))];

disp('Volume Boundaries:');
disp(['    X: ', num2str(xLimsData(1), ['%+.', num2str(spacePrecision), 'f']), ' [m] -> ', num2str(xLimsData(2), ['%+.', num2str(spacePrecision), 'f']), ' [m]']);
disp(['    Y: ', num2str(yLimsData(1), ['%+.', num2str(spacePrecision), 'f']), ' [m] -> ', num2str(yLimsData(2), ['%+.', num2str(spacePrecision), 'f']), ' [m]']);
disp(['    Z: ', num2str(zLimsData(1), ['%+.', num2str(spacePrecision), 'f']), ' [m] -> ', num2str(zLimsData(2), ['%+.', num2str(spacePrecision), 'f']), ' [m]']);

% Select Target Point
LOSdata.targetPoint = zeros(1,3);

valid = false;
while ~valid
    disp(' ')

    disp('Specify Target Point:')
    
    LOSdata.targetPoint(1) = inputPos('X');
    LOSdata.targetPoint(2) = inputPos('Y');
    LOSdata.targetPoint(3) = inputPos('Z');

    if (LOSdata.targetPoint(1) < xLimsData(1) || LOSdata.targetPoint(1) > xLimsData(2)) || ...
       (LOSdata.targetPoint(2) < yLimsData(1) || LOSdata.targetPoint(2) > yLimsData(2)) || ...
       (LOSdata.targetPoint(3) < zLimsData(1) || LOSdata.targetPoint(3) > zLimsData(2))
        disp('        WARNING: Point Lies Outside Volume');
        disp(' ');
    else
        valid = true;
    end

end

disp(' ');

% Shift Target Point To Nearest Volume Node
[offset, index] = min(abs(LOSdata.targetPoint - volumeData.positionGrid));

if any(offset ~= 0)
    disp('Shifting Point To Nearest Volume Node:');
        
    LOSdata.targetPoint(1) = volumeData.positionGrid(index(1),1);
    LOSdata.targetPoint(2) = volumeData.positionGrid(index(2),2);
    LOSdata.targetPoint(3) = volumeData.positionGrid(index(3),3);
    
    disp(['    X: ', num2str(LOSdata.targetPoint(1), ['%+.', num2str(spacePrecision), 'f']), ' [m]']);
    disp(['    Y: ', num2str(LOSdata.targetPoint(2), ['%+.', num2str(spacePrecision), 'f']), ' [m]']);
    disp(['    Z: ', num2str(LOSdata.targetPoint(3), ['%+.', num2str(spacePrecision), 'f']), ' [m]']);
end

disp(' ');
disp(' ');


%% Select Reference Plane

disp('Reference Plane Definition');
disp('---------------------------');

% Select Reference Plane
LOSdata.orientation = 'XY';

valid = false;
while ~valid
    disp(' ')

    disp('Specify Downstream Reference Plane:');
    
    LOSdata.position = inputPos('X');

    if LOSdata.position <= LOSdata.targetPoint(1)
        disp('        WARNING: Plane Lies Upstream of Target Point');
    elseif LOSdata.position > xLimsData(2)
        disp('        WARNING: Plane Lies Outside Volume');
        disp(' ');
    else
        valid = true;
    end

end

disp(' ');

% Shift Reference Plane To Nearest Volume Slice
[offset, index] = min(abs(LOSdata.position - volumeData.positionGrid(:,1)));

if offset ~= 0
    disp('Shifting Plane To Nearest Volume Slice');

    LOSdata.position = volumeData.positionGrid(index,1);
end

% Re-Normalise Dimensions
if normalise

    if contains(caseName, ["Run_Test", "Windsor"])
        volumeData.positionGrid = round((volumeData.positionGrid / 1.044), spacePrecision);
        LOSdata.targetPoint = round((LOSdata.targetPoint / 1.044), spacePrecision);
        LOSdata.position = round((LOSdata.position / 1.044), spacePrecision);
    end

end

% Store Reference Plane
index = find(volumeData.positionGrid(:,1) == LOSdata.position);
LOSdata.positionGrid = volumeData.positionGrid(index,:);

disp(' ');
disp(' ');


%% Calculate Line of Sight

disp('Line of Slight Calculation');
disp('---------------------------');

disp(' ');

disp('***********');
disp('  RUNNING ');

% tic;

disp(' ');

disp('    Initialising...');

% Reduce Volume to Required Limits
index = find((volumeData.positionGrid(:,1) >= LOSdata.targetPoint(1)) & (volumeData.positionGrid(:,1) <= LOSdata.position));
volumeData.positionGrid = volumeData.positionGrid(index,:);

fields = fieldnames(volumeData.mean);

tic;
for i = 1:height(fields)
    volumeData.mean.(fields{i}) = volumeData.mean.(fields{i})(index,:);

    for j = 1:height(volumeData.inst.time)
        volumeData.inst.(fields{i}){j} = volumeData.inst.(fields{i}){j}(index,:);
    end

end
toc;

disp(' ');

% Calculate True Line of Sight
disp('    Calculating True Line of Sight...');

slices = unique(volumeData.positionGrid(:,1));
LOSdata.LOSmatrix = zeros(height(LOSdata.positionGrid), width(LOSdata.positionGrid), height(slices));

tic;
for i = 1:height(LOSdata.LOSmatrix) % https://math.stackexchange.com/questions/576137/finding-a-point-on-a-3d-line
    t = (slices - LOSdata.targetPoint(1)) / (LOSdata.positionGrid(i,1) - LOSdata.targetPoint(1));
    LOSdata.LOSmatrix(i,:,:) = (LOSdata.targetPoint + t*(LOSdata.positionGrid(i,:) - LOSdata.targetPoint))';
end
toc;

disp(' ');

% Shift Line of Sight to Pass Through Local Volume Nodes
disp('    Shifting Line of Sight to Pass Through Local Volume Nodes...');

LOSdata.LOSmatrixShifted = LOSdata.LOSmatrix;
LOSdata.LOSindex = zeros(height(LOSdata.LOSmatrix),depth(LOSdata.LOSmatrix));

% Initialise Progress Bar
wB = waitbar(0, 'Shifting Line of Sight to Pass Through Local Volume Nodes', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';

tic;
for i = 1:depth(LOSdata.LOSmatrix)
    slicePositionGrid = volumeData.positionGrid((volumeData.positionGrid(:,1) == slices(i)),:);

    for j = 1:height(LOSdata.LOSmatrix)
        [~, index] = min(abs(slicePositionGrid - LOSdata.LOSmatrix(j,:,i)));
        LOSdata.LOSindex(j,i) = find((volumeData.positionGrid(:,1) == slicePositionGrid(index(1),1)) & ...
                                      (volumeData.positionGrid(:,2) == slicePositionGrid(index(2),2)) & ...
                                      (volumeData.positionGrid(:,3) == slicePositionGrid(index(3),3)));

        LOSdata.LOSmatrixShifted(j,:,i) = volumeData.positionGrid(LOSdata.LOSindex(j,i),:);
    end

    waitbar((i / depth(LOSdata.LOSmatrix)), wB);
end
clear slicePositionGrid index;
toc;

delete(wB);

% Verify Line of Sight Shift
testLine = randi(height(LOSdata.positionGrid),[10,1]);
figure;
hold on;

for i = 1:10
    plot3(squeeze(LOSdata.LOSmatrix(testLine(i),1,:)), squeeze(LOSdata.LOSmatrix(testLine(i),2,:)), squeeze(LOSdata.LOSmatrix(testLine(i),3,:)), 'r-o');
    plot3(squeeze(LOSdata.LOSmatrixShifted(testLine(i),1,:)), squeeze(LOSdata.LOSmatrixShifted(testLine(i),2,:)), squeeze(LOSdata.LOSmatrixShifted(testLine(i),3,:)), 'b-x');
end

xLimsData = [min(volumeData.positionGrid(:,1)); max(volumeData.positionGrid(:,1))];
yLimsData = [min(volumeData.positionGrid(:,2)); max(volumeData.positionGrid(:,2))];
zLimsData = [min(volumeData.positionGrid(:,3)); max(volumeData.positionGrid(:,3))];
xlim(xLimsData);
ylim(yLimsData);
zlim(zLimsData);
view(30,30)
hold off;

% % Calculate Obstructing Mass Along Line of Site
% disp('    Calculating Obstructing Mass Along Line of Site...');
% 
% LOSdata.inst.time = volumeData.inst.time;
% LOSdata.inst.mass = cell(height(LOSdata.inst.time),1);
% LOSdata.mean.mass = zeros(height(LOSdata.LOSmatrixShifted),1);
% 
% % Initialise Progress Bar
% wB = waitbar(0, 'Shifting Line of Sight to Pass Through Local Volume Nodes', 'name', 'Progress');
% wB.Children.Title.Interpreter = 'none';
% 
% tic;
% for i = 1:height(LOSdata.LOSmatrixShifted)
%     obstructingMass(i) = sum(volumeData.mean.mass(index(i,:)));
% 
%     waitbar((i / height(LOSdata.LOSmatrixShifted)), wB);
% end
% toc;
% 
% delete(wB);

% executionTime = toc;

disp(' ');

% disp(['    Run Time: ', num2str(executionTime), 's']);

disp(' ');

disp('  SUCCESS  ');
disp('***********');

disp(' ');
disp(' ');


%% Interpolation Test Code

clc;

clear interp cellSize sX i t samplePoints;

disp(' ');

%

tic;

interp = scatteredInterpolant(volumeData.positionGrid(:,1), ...
                   volumeData.positionGrid(:,2), ...
                   volumeData.positionGrid(:,3), ...
                   volumeData.mean.mass, 'linear', 'none');

toc;

disp(' ');

%

tic;

cellSize = 8e-3;
dX = (LOSdata.targetPoint(1):(cellSize / 4):LOSdata.position)';

mass = zeros(height(LOSdata.positionGrid),1);

for i = 1:height(LOSdata.positionGrid)
        t = (dX - LOSdata.targetPoint(1)) / (LOSdata.positionGrid(i,1) - LOSdata.targetPoint(1));
        samplePoints = (LOSdata.targetPoint + t*(LOSdata.positionGrid(i,:) - LOSdata.targetPoint));
        mass(i) = sum(interp(samplePoints(:,1), samplePoints(:,2), samplePoints(:,3)));
end

toc;


%% Local Functions

function pos = inputPos(orientation)

    pos = str2double(input(['    ', orientation, '-Position [m]: '], 's'));
    
    if isnan(pos) || length(pos) > 1
        disp('        WARNING: Invalid Entry');
        pos = -1;
    end
    
end