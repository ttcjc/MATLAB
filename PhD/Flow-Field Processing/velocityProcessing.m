%% Velocity Processing v4.0

clear variables;
close all;
clc;

fig = 0;
figHold = 0;

nProc = 4; % Number of Processors Used for Parallel Collation

disp ('========================');
disp ('Velocity Processing v4.0');
disp ('========================');

disp (' ');
disp (' ');


%% Changelog

% v1.0 - Initial Commit
% v2.0 - Rewrite to Improve Versatility
% v3.0 - Moved Plotting to velocityPlots.m
% v4.0 - Rewrite to Support ParaView, Probe and Experimental Planar Data


%% Select Data Format

disp('Data Format');
disp('------------');

disp(' ');

disp('Compatible Data Formats:');
disp('    A: ParaView Planar Data');
disp('    B: Probe Planar Data');
disp('    C: Experimental Planar Data (Varney)');

valid = false;
while ~valid
    disp(' ');
    selection = input('Select Data Format [A/B/C]: ', 's');

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

disp (' ');
disp (' ');


%% Initialisation

switch format

    case 'A'
        [caseFolder, data, geometry, xDims, yDims, zDims, precision] = initialisePVdata('U', true);

    case 'B'
%         [caseFolder, data, geometry, xDims, yDims, zDims, precision] = initialiseProbeData('U', true, nProc);

    case 'C'
        [campaign, data, geometry, xDims, yDims, zDims, precision] = initialiseExpData('U', true);

end


%% Data Formatting

planes = fieldnames(data);

switch format

    case 'A'        
        % Adjust Data Origin
        if contains(caseFolder, ["Run_Test", "Windsor"]) && contains(caseFolder, 'Upstream')
            
            for i = 1:height(planes)
                data.(planes{i}).position(:,1) = data.(planes{i}).position(:,1) + round(1.325 / 1.044, precision);

                if strcmp(data.(planes{i}).planeOrientation, 'X')
                    data.(planes{i}).planePosition(:,1) = data.(planes{i}).planePosition(:,1) + round(1.325 / 1.044, precision);
                end

            end
            
        end
        
        % Normalise Velocity
        if contains(caseFolder, ["Run_Test", "Windsor"])
            
            for i = 1:height(planes)
                data.(planes{i}).u = data.(planes{i}).u / 40;
                data.(planes{i}).v = data.(planes{i}).v / 40;
                data.(planes{i}).w = data.(planes{i}).w / 40;
            end
            
        end
        
    case 'B'
%         % Adjust Data Origin
%         if contains(caseFolder, ["Run_Test", "Windsor"]) && contains(caseFolder, 'Upstream')
%             data.position(:,1) = data.position(:,1) + 1.325;
%             data..planePosition(:,1) = data..planePosition(:,1) + 1.325;
%         end
%         
%         % Normalise Velocity
%         if contains(caseFolder, ["Run_Test", "Windsor"])
%             data.uMean = data.uMean / 40;
%             data.vMean = data.vMean / 40;
%             data.wMean = data.wMean / 40;
% 
%             for i = 1:height(data.time)
%                 data.u{i} = data.u{i} / 40;
%                 data.v{i} = data.v{i} / 40;
%                 data.w{i} = data.w{i} / 40;
%                 data.uPrime{i} = data.uPrime{i} / 40;
%                 data.vPrime{i} = data.vPrime{i} / 40;
%                 data.wPrime{i} = data.wPrime{i} / 40;
%             end
%             
%         end
    
    case 'C'
        % Normalise Velocity
        if contains(campaign, 'Varney')
            
            for i = 1:height(planes)
                    data.(planes{i}).uMean = data.(planes{i}).uMean / 40;
                    data.(planes{i}).vMean = data.(planes{i}).vMean / 40;
                    data.(planes{i}).wMean = data.(planes{i}).wMean / 40;
                    data.(planes{i}).uRMS = data.(planes{i}).uRMS / 40;
                    data.(planes{i}).vRMS = data.(planes{i}).vRMS / 40;
                    data.(planes{i}).wRMS = data.(planes{i}).wRMS / 40;
            end
            
        end
        
end

disp (' ');
disp (' ');


%% Data Presentation

disp('Data Presentation');
disp('------------------');

disp(' ');

disp('Available Planes:');

for i = 1:height(planes)
    disp(['    ', num2str(i), '. ', planes{i}]);
end

valid = false;
while ~valid
    disp(' ');
    selection = input('Plot All Available Planes? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        plotPlanes = inputPlanes;
        
        if (height(plotPlanes) > 1 && min(plotPlanes) == 0) || ...
           min(plotPlanes) < 0 || ...
           max(plotPlanes) > height(planes)
            disp('        WARNING: Invalid Plane Selection');
        else
            valid = true;
        end
        
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        plotPlanes = 1:height(planes);
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end  
    
end

disp(' ');

if ~plotPlanes
    disp('    Skipping Plane Presentation');
    return
end

for i = plotPlanes
    disp(['    Presenting ', planes{i}, '...']);
    
    planeOrientation = data.(planes{i}).planeOrientation;
    caseType = 'Windsor';
    normalise = true;
%     precision = precision;
    xLimsData = [min(data.(planes{i}).position(:,1)); max(data.(planes{i}).position(:,1))];
    yLimsData = [min(data.(planes{i}).position(:,2)); max(data.(planes{i}).position(:,2))];
    zLimsData = [min(data.(planes{i}).position(:,3)); max(data.(planes{i}).position(:,3))];
    planePosition = data.(planes{i}).planePosition;
    positionData = data.(planes{i}).position;
    vectorData = [data.(planes{i}).u, data.(planes{i}).v, data.(planes{i}).w];
    nComponents = 3;
    component = [];
%     geometry = geometry;
%     fig = fig;
    figName = planes{i};
    cMap = viridis(24);
    streamlines = true;
%     xDims = xDims;
%     yDims = yDims;
%     zDims = zDims;
    figTitle = ' ';
    cLims = [0, 1];
    
    fig = vectorPlots(planeOrientation, caseType, normalise, precision, ...
                      xLimsData, yLimsData, zLimsData, planePosition, positionData, vectorData, ...
                      nComponents, component, geometry, fig, figName, ...
                      cMap, streamlines, xDims, yDims, zDims, figTitle, cLims);
end    
    

%% Local Functions

function [P] = inputPlanes

    valid = false;
    while ~valid
        P = str2num(input('    List Desired Planes (Row Vector Form): ', 's')); %#ok<ST2NM>

        if any(isnan(P)) || ~isrow(P)
            disp('        WARNING: Invalid Entry');
        else
            valid = true;
        end

    end

end
