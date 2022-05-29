%% Velocity Processing v4.0

clear variables;
close all;
clc;
evalc('delete(gcp(''nocreate''));');

normalise = true; % Normalisation of Dimensions

nProc = maxNumCompThreads - 2; % Number of Processors Used for Parallel Collation

fig = 0; % Initialise Figure Tracking
figHold = 0; % Enable Overwriting of Figures

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
        [caseName, data, geometry, xDims, yDims, zDims, spacePrecision] = initialisePVdata('U', normalise);

    case 'B'
        [caseName, data, geometry, xDims, yDims, zDims, spacePrecision] = initialiseProbeData('U', true, normalise, nProc);    

    case 'C'
        [caseName, data, geometry, xDims, yDims, zDims, spacePrecision] = initialiseExpData('U', normalise);

end


%% Data Formatting

planes = fieldnames(data);

switch format

    case 'A'       
        % Adjust Data Origin
        if contains(caseName, 'Run_Test') || (contains(caseName, 'Windsor') && contains(caseName, 'Upstream'))
            
            for i = 1:height(planes)
                
                if normalise
                    data.(planes{i}).position(:,1) = data.(planes{i}).position(:,1) + round((1.325 / 1.044), spacePrecision);

                    if strcmp(data.(planes{i}).planeOrientation, 'YZ')
                        data.(planes{i}).planePosition(:,1) = data.(planes{i}).planePosition(:,1) + round((1.325 / 1.044), spacePrecision);
                    end
                    
                else
                    data.(planes{i}).position(:,1) = data.(planes{i}).position(:,1) + 1.325; %#ok<*UNRCH>

                    if strcmp(data.(planes{i}).planeOrientation, 'YZ')
                        data.(planes{i}).planePosition(:,1) = data.(planes{i}).planePosition(:,1) + 1.325;
                    end
                    
                end

            end
            
        end
        
        % Normalise Velocity
        if contains(caseName, ["Run_Test", "Windsor"])
            U = 40;
            
            for i = 1:height(planes)
                data.(planes{i}).uMean = data.(planes{i}).uMean / U;
                data.(planes{i}).vMean = data.(planes{i}).vMean / U;
                data.(planes{i}).wMean = data.(planes{i}).wMean / U;
            end
            
        end
        
    case 'B'
        % Adjust Data Origin
        if contains(caseName, 'Run_Test') || (contains(caseName, 'Windsor') && contains(caseName, 'Upstream'))
            
            for i = 1:height(planes)
                
                if normalise
                    data.(planes{i}).position(:,1) = data.(planes{i}).position(:,1) + round((1.325 / 1.044), spacePrecision);

                    if strcmp(data.(planes{i}).planeOrientation, 'YZ')
                        data.(planes{i}).planePosition(:,1) = data.(planes{i}).planePosition(:,1) + round((1.325 / 1.044), spacePrecision);
                    end
                    
                else
                    data.(planes{i}).position(:,1) = data.(planes{i}).position(:,1) + 1.325;

                    if strcmp(data.(planes{i}).planeOrientation, 'YZ')
                        data.(planes{i}).planePosition(:,1) = data.(planes{i}).planePosition(:,1) + 1.325;
                    end
                    
                end

                % Normalise Velocity
                if contains(caseName, ["Run_Test", "Windsor"])
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
        if strcmp(caseName, 'Varney')
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

% Select Plane(s) of Interest
valid = false;
while ~valid
    [index, valid] = listdlg('listSize', [300, 300], ...
                             'selectionMode', 'multiple', ...
                             'name', 'Select Variable(s) to Plot', ...
                             'listString', planes);

    if ~valid
        disp(' ');
        disp('WARNING: No Planes Selected');
    end
end
clear valid;

plotPlanes = planes(index);

switch format
    
    case 'B'
        valid = false;
        while ~valid
            disp(' ');
            selection = input('Plot Instantaneous Data? [y/n]: ', 's');

            if selection == 'n' | selection == 'N' %#ok<OR2>
                plotInst = false;
                
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

for i = 1:height(plotPlanes)
    disp(['    Presenting ', plotPlanes{i}, '...']);
    
    planeOrientation = data.(plotPlanes{i}).planeOrientation;
    
    switch planeOrientation

        case 'YZ'

            if contains(caseName, ["Run_Test", "Windsor", "Varney"])
                xLimsPlot = [0.31875; 4.65925]; % [m]
                yLimsPlot = [-0.5945; 0.5945];
                zLimsPlot = [0; 0.739];
            end

        case {'XZ', 'XY'}

            if contains(caseName, ["Run_Test", "Windsor", "Varney"])
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
    
    switch format
        
        case 'A'
            xLimsData = xLimsPlot;
            yLimsData = yLimsPlot;
            zLimsData = zLimsPlot;
            
        case 'B'
            xLimsData = data.(plotPlanes{i}).xLims;
            yLimsData = data.(plotPlanes{i}).yLims;
            zLimsData = data.(plotPlanes{i}).zLims;
            
            xLimsPlot = [min(xLimsPlot(1), xLimsData(1)); max(xLimsPlot(2), xLimsData(2))];
            yLimsPlot = [min(yLimsPlot(1), yLimsData(1)); max(yLimsPlot(2), yLimsData(2))];
            zLimsPlot = [min(zLimsPlot(1), zLimsData(1)); max(zLimsPlot(2), zLimsData(2))];
        
        case 'C'
            xLimsData = [min(data.(planes{i}).position(:,1)); max(data.(planes{i}).position(:,1))];
            yLimsData = [min(data.(planes{i}).position(:,2)); max(data.(planes{i}).position(:,2))];
            zLimsData = [min(data.(planes{i}).position(:,3)); max(data.(planes{i}).position(:,3))];
            
    end
    
    planePosition = data.(plotPlanes{i}).planePosition;
    positionData = data.(plotPlanes{i}).position;
    vectorData = [data.(plotPlanes{i}).uMean, data.(planes{i}).vMean, data.(plotPlanes{i}).wMean];
    
    if strcmp(caseName, 'Varney') && ~strcmp(planeOrientation, 'YZ')
        nComponents = 2;
    else
        nComponents = 3;
    end
    
    component = [];
    figName = [caseName, '_', plotPlanes{i}];
    cMap = viridis(24);
    streamlines = true;
    figTitle = '-'; % Leave Blank ('-') for Formatting Purposes
    figSubtitle = ' ';
    cLims = [0, 1];
    
    fig = vectorPlots(xLimsPlot, yLimsPlot, zLimsPlot, xLimsData, yLimsData, zLimsData, ...
                      planeOrientation, planePosition, positionData, vectorData, ...
                      nComponents, component, fig, figName, cMap, geometry, streamlines, ...
                      xDims, yDims, zDims, figTitle, figSubtitle, cLims, normalise);
    
%     switch format
%         
%         case 'B'
%             
%             if plotInst
%                 figHold = fig;
%                 
%                 for j = 1:height(data.(plotPlanes{i}).time)
%                     
%                     if j ~= 1
%                         clf(fig)
%                     end
%                     
%                     figTime = num2str(mapData.inst.time(j), ['%.', num2str(timePrecision), 'f']);
%                     
%                     vectorData = [data.(plotPlanes{i}).u{j}, data.(plotPlanes{i}).v{j}, data.(plotPlanes{i}).w{j}];
%                     fig = figHold;
%                     figName = [caseType, '_', plotPlanes{i}, '_T', erase(figTime, '.')];
%                     figSubtitle = [figTime, ' \it{s}'];
%                     
%                     fig = vectorPlots(xLimsPlot, yLimsPlot, zLimsPlot, xLimsData, yLimsData, zLimsData, ...
%                   planeOrientation, planePosition, positionData, vectorData, ...
%                   nComponents, component, fig, figName, cMap, geometry, streamlines, ...
%                   xDims, yDims, zDims, figTitle, figSubtitle, cLims, normalise);
%                 end
%                 
%             end
%             
%     end
    
end