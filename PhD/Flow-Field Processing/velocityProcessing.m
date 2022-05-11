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
        [caseFolder, data, geometry, xDims, yDims, zDims, precision] = initialiseProbeData('U', true, true, nProc);    

    case 'C'
        [campaign, data, geometry, xDims, yDims, zDims, precision] = initialiseExpData('U', true);

end

switch format
    
    case {'A', 'B'}
        
        if contains(caseFolder, 'Run_Test')
            caseType = 'Run_Test';
        elseif contains(caseFolder, 'Windsor')
            caseType = 'Windsor';
        end
        
    case 'C'
        caseType = campaign;
        
end        


%% Data Formatting

planes = fieldnames(data);

switch format

    case 'A'       
        % Adjust Data Origin
        if (strcmp(caseType, 'Run_Test') || strcmp(caseType, 'Windsor')) && contains(caseFolder, 'Upstream')
            
            for i = 1:height(planes)
                
                if normalise
                    data.(planes{i}).position(:,1) = data.(planes{i}).position(:,1) + round(1.325 / 1.044, precision);

                    if strcmp(data.(planes{i}).planeOrientation, 'X')
                        data.(planes{i}).planePosition(:,1) = data.(planes{i}).planePosition(:,1) + round(1.325 / 1.044, precision);
                    end
                    
                else
                    data.(planes{i}).position(:,1) = data.(planes{i}).position(:,1) + 1.325; %#ok<*UNRCH>

                    if strcmp(data.(planes{i}).planeOrientation, 'X')
                        data.(planes{i}).planePosition(:,1) = data.(planes{i}).planePosition(:,1) + 1.325;
                    end
                    
                end

            end
            
        end
        
        % Normalise Velocity
        if strcmp(caseType, 'Run_Test') || strcmp(caseType, 'Windsor')
            U = 40;
            
            for i = 1:height(planes)
                data.(planes{i}).uMean = data.(planes{i}).uMean / U;
                data.(planes{i}).vMean = data.(planes{i}).vMean / U;
                data.(planes{i}).wMean = data.(planes{i}).wMean / U;
            end
            
        end
        
    case 'B'
        % Adjust Data Origin
        if (strcmp(caseType, 'Run_Test') || strcmp(caseType, 'Windsor')) && contains(caseFolder, 'Upstream')
            
            for i = 1:height(planes)
                
                if normalise
                    data.(planes{i}).position(:,1) = data.(planes{i}).position(:,1) + round(1.325 / 1.044, precision);

                    if strcmp(data.(planes{i}).planeOrientation, 'X')
                        data.(planes{i}).planePosition(:,1) = data.(planes{i}).planePosition(:,1) + round(1.325 / 1.044, precision);
                    end
                    
                else
                    data.(planes{i}).position(:,1) = data.(planes{i}).position(:,1) + 1.325;

                    if strcmp(data.(planes{i}).planeOrientation, 'X')
                        data.(planes{i}).planePosition(:,1) = data.(planes{i}).planePosition(:,1) + 1.325;
                    end
                    
                end

                % Normalise Velocity
                if strcmp(caseType, 'Run_Test') || strcmp(caseType, 'Windsor')
                    U = 40;
                    
                    data.(planes{i}).uMean = data.(planes{i}).uMean / U;
                    data.(planes{i}).vMean = data.(planes{i}).vMean / U;
                    data.(planes{i}).wMean = data.(planes{i}).wMean / U;

                    for j = 1:height(data.(planes{i}).time)
                        data.(planes{i}).u{j} = data.(planes{i}).u{j} / U;
                        data.(planes{i}).v{j} = data.(planes{i}).v{j} / U;
                        data.(planes{i}).w{j} = data.(planes{i}).w{j} / U;
                        data.(planes{i}).uPrime{j} = data.(planes{i}).uPrime{j} / U;
                        data.(planes{i}).vPrime{j} = data.(planes{i}).vPrime{j} / U;
                        data.(planes{i}).wPrime{j} = data.(planes{i}).wPrime{j} / U;
                    end

                end

            end
            
        end
        
    case 'C'
        % Normalise Velocity
        if strcmp(campaign, 'Varney')
            U = 40;
            
            for i = 1:height(planes)
                    data.(planes{i}).uMean = data.(planes{i}).uMean / U;
                    data.(planes{i}).vMean = data.(planes{i}).vMean / U;
                    data.(planes{i}).wMean = data.(planes{i}).wMean / U;
                    data.(planes{i}).uRMS = data.(planes{i}).uRMS / U;
                    data.(planes{i}).vRMS = data.(planes{i}).vRMS / U;
                    data.(planes{i}).wRMS = data.(planes{i}).wRMS / U;
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

if strcmp(format, 'B')
    
    valid = false;
    while ~valid
        disp(' ');
        selection = input('Plot Instantaneous Data? [y/n]: ', 's');

        if selection == 'n' | selection == 'N' %#ok<OR2>
            valid = true;
        elseif selection == 'y' | selection == 'Y' %#ok<OR2>
            plotInst = true;
            valid = true;
        else
            disp('    WARNING: Invalid Entry');
        end  

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
    
    switch planeOrientation

        case 'X'

            if contains(caseType, ["Run_Test", "Windsor", "Varney"])
                xLimsPlot = [0.31875; 4.65925]; % [m]
                yLimsPlot = [-0.4945; 0.4945];
                zLimsPlot = [0; 0.639];
            end

        case {'Y', 'Z'}

            if contains(caseType, ["Run_Test", "Windsor", "Varney"])
                xLimsPlot = [0.31875; 1.075]; % [m]
                yLimsPlot = [-0.3445; 0.3445];
                zLimsPlot = [0; 0.489];
            end

    end

    if normalise
        xLimsPlot = round(xLimsPlot / 1.044, precision);
        yLimsPlot = round(yLimsPlot / 1.044, precision);
        zLimsPlot = round(zLimsPlot / 1.044, precision);
    end
    
    switch format
        
        case 'A'
            xLimsData = xLimsPlot;
            yLimsData = yLimsPlot;
            zLimsData = zLimsPlot;
            
        case 'B'
            xLimsData = data.(planes{i}).xLims;
            yLimsData = data.(planes{i}).yLims;
            zLimsData = data.(planes{i}).zLims;
            
            xLimsPlot = [min(xLimsPlot(1), xLimsData(1)); max(xLimsPlot(2), xLimsData(2))];
            yLimsPlot = [min(yLimsPlot(1), yLimsData(1)); max(yLimsPlot(2), yLimsData(2))];
            zLimsPlot = [min(zLimsPlot(1), zLimsData(1)); max(zLimsPlot(2), zLimsData(2))];
        
        case 'C'
            xLimsData = [min(data.(planes{i}).position(:,1)); max(data.(planes{i}).position(:,1))];
            yLimsData = [min(data.(planes{i}).position(:,2)); max(data.(planes{i}).position(:,2))];
            zLimsData = [min(data.(planes{i}).position(:,3)); max(data.(planes{i}).position(:,3))];
            
    end
    
    planePosition = data.(planes{i}).planePosition;
    positionData = data.(planes{i}).position;
    vectorData = [data.(planes{i}).uMean, data.(planes{i}).vMean, data.(planes{i}).wMean];
    
    if strcmp(caseType, 'Varney') && ~strcmp(planeOrientation, 'X')
        nComponents = 2;
    else
        nComponents = 3;
    end
    
    component = [];
    figName = [caseType, '_', planes{i}];
    cMap = viridis(24);
    streamlines = true;
    figTitle = '-'; % Leave Blank ('-') for Formatting Purposes
    figSubtitle = ' ';
    cLims = [0, 1];
    
    fig = vectorPlots(xLimsPlot, yLimsPlot, zLimsPlot, xLimsData, yLimsData, zLimsData, ...
                      planeOrientation, planePosition, positionData, vectorData, ...
                      nComponents, component, fig, figName, cMap, geometry, streamlines, ...
                      xDims, yDims, zDims, figTitle, figSubtitle, cLims, normalise);
    
    if exist('plotInst', 'var')
        figHold = fig;
        
        for j = 1:3 % height(data.(planes{i}).time)
            
            if j ~= 1
                clf(fig)
            end
            
            figTime = num2str(data.(planes{i}).time(j), '%.4f');
            
            vectorData = [data.(planes{i}).u{j}, data.(planes{i}).v{j}, data.(planes{i}).w{j}];
            fig = figHold;
            figName = [caseType, '_', planes{i}, '_T', erase(figTime, '.')];
            figSubtitle = [figTime, ' \it{s}'];
            

            fig = vectorPlots(xLimsPlot, yLimsPlot, zLimsPlot, xLimsData, yLimsData, zLimsData, ...
                  planeOrientation, planePosition, positionData, vectorData, ...
                  nComponents, component, fig, figName, cMap, geometry, streamlines, ...
                  xDims, yDims, zDims, figTitle, figSubtitle, cLims, normalise);
        end
        
    end
    
end
    

%% Local Functions

function P = inputPlanes

    valid = false;
    while ~valid
        P = str2num(input('    List Desired Planes [Row Vector]: ', 's')); %#ok<ST2NM>

        if any(isnan(P)) || ~isrow(P)
            disp('        WARNING: Invalid Entry');
        else
            valid = true;
        end

    end

end