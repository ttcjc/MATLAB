%% Planar Velocity Processing v4.0

clear variables;
close all;
clc;
evalc('delete(gcp(''nocreate''));');

normalise = true; % Normalisation of Dimensions

nProc = maxNumCompThreads - 2; % Number of Processors Used for Parallel Collation

fig = 0; % Initialise Figure Tracking
figHold = 0; % Enable Overwriting of Figures

disp('===============================');
disp('Planar Velocity Processing v4.0');
disp('===============================');

disp(' ');
disp(' ');


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
clear valid;

disp(' ');
disp(' ');


%% Initialisation

switch format

    case 'A'
        [caseName, velData, geometry, ...
         xDims, yDims, zDims, spacePrecision] = initialisePVdata('U', normalise);

    case 'B'
        [caseName, dataID, velData, timePrecision, geometry, ...
         xDims, yDims, zDims, spacePrecision] = initialiseVelocityProbeData('planar', normalise, nProc);    

    case 'C'
        [caseName, velData, geometry, ...
         xDims, yDims, zDims, spacePrecision] = initialiseExpData('U', normalise);

end


%% Data Formatting

planes = fieldnames(velData);

% Adjust Data Origin
switch format

    case 'A'
        
        if contains(caseName, 'Run_Test') || (contains(caseName, 'Windsor') && contains(caseName, 'Upstream'))
            
            for i = 1:height(planes)
                
                if normalise
                    velData.(planes{i}).position(:,1) = velData.(planes{i}).position(:,1) + round((1.325 / 1.044), spacePrecision);

                    if strcmp(velData.(planes{i}).planeOrientation, 'YZ')
                        velData.(planes{i}).planePosition = velData.(planes{i}).planePosition + round((1.325 / 1.044), spacePrecision);
                    end
                    
                else
                    velData.(planes{i}).position(:,1) = velData.(planes{i}).position(:,1) + 1.325; %#ok<*UNRCH>

                    if strcmp(velData.(planes{i}).planeOrientation, 'YZ')
                        velData.(planes{i}).planePosition = velData.(planes{i}).planePosition + 1.325;
                    end
                    
                end

            end
            
        end
        
    case 'B'
        
        if contains(caseName, 'Run_Test') || (contains(caseName, 'Windsor') && contains(caseName, 'Upstream'))
            
            for i = 1:height(planes)
                
                if normalise
                    velData.(planes{i}).position(:,1) = velData.(planes{i}).position(:,1) + round((1.325 / 1.044), spacePrecision);
                    
                    switch velData.(planes{i}).planeOrientation
                        
                        case 'YZ'
                            velData.(planes{i}).planePosition = velData.(planes{i}).planePosition + round((1.325 / 1.044), spacePrecision);
                            velData.(planes{i}).xLims = velData.(planes{i}).xLims + round((1.325 / 1.044), spacePrecision);
                            
                        case {'XZ', 'XY'}
                            velData.(planes{i}).xLims = velData.(planes{i}).xLims + round((1.325 / 1.044), spacePrecision);
                            
                    end
                    
                else
                    velData.(planes{i}).position(:,1) = velData.(planes{i}).position(:,1) + 1.325; %#ok<*UNRCH>
                    
                    switch velData.(planes{i}).planeOrientation
                        
                        case 'YZ'
                            velData.(planes{i}).planePosition = velData.(planes{i}).planePosition + 1.325;
                            velData.(planes{i}).xLims = velData.(planes{i}).xLims + 1.325;
                        
                        case {'XZ', 'XY'}
                            velData.(planes{i}).xLims = velData.(planes{i}).xLims + 1.325;
                    
                    end
                    
                end

            end
            
        end
        
end

% Normalise Velocity
switch format
    
    case 'A'
        
        if contains(caseName, ["Run_Test", "Windsor"])
            U = 40; % m/s
            
            for i = 1:height(planes)
                velData.(planes{i}).uMean = velData.(planes{i}).uMean / U;
                velData.(planes{i}).vMean = velData.(planes{i}).vMean / U;
                velData.(planes{i}).wMean = velData.(planes{i}).wMean / U;
            end
            
        end
        
    case 'B'
        
        if contains(caseName, ["Run_Test", "Windsor"])
            U = 40; % m/s
            
            for i = 1:height(planes)
                velData.(planes{i}).uMean = velData.(planes{i}).uMean / U;
                velData.(planes{i}).vMean = velData.(planes{i}).vMean / U;
                velData.(planes{i}).wMean = velData.(planes{i}).wMean / U;
                
                for j = 1:height(velData.(planes{i}).time)
                    velData.(planes{i}).u{j} = velData.(planes{i}).u{j} / U;
                    velData.(planes{i}).v{j} = velData.(planes{i}).v{j} / U;
                    velData.(planes{i}).w{j} = velData.(planes{i}).w{j} / U;
                    velData.(planes{i}).uPrime{j} = velData.(planes{i}).uPrime{j} / U;
                    velData.(planes{i}).vPrime{j} = velData.(planes{i}).vPrime{j} / U;
                    velData.(planes{i}).wPrime{j} = velData.(planes{i}).wPrime{j} / U;
                end
                
            end
            
        end
        
    case 'C'
        
        if strcmp(caseName, 'Varney')
            U = 40; % m/s
            
            for i = 1:height(planes)
                    velData.(planes{i}).uMean = velData.(planes{i}).uMean / U;
                    velData.(planes{i}).vMean = velData.(planes{i}).vMean / U;
                    velData.(planes{i}).wMean = velData.(planes{i}).wMean / U;
                    velData.(planes{i}).uRMS = velData.(planes{i}).uRMS / U;
                    velData.(planes{i}).vRMS = velData.(planes{i}).vRMS / U;
                    velData.(planes{i}).wRMS = velData.(planes{i}).wRMS / U;
            end
            
        end
        
end

% Perform Blockage Correction
switch format
    
    case 'A'
        
        if contains(caseName, 'Run_Test') || (contains(caseName, 'Windsor') && contains(caseName, 'Upstream'))
            Am = (0.289 * 0.389) + (2 * (0.046 * 0.055));
            At = (2 * (0.9519083 + (3.283 * tan(atan(0.0262223 / 9.44)))) * 1.32);
            
            for i = 1:height(planes)
                velData.(planes{i}).uMean = velData.(planes{i}).uMean * (At / (At - Am));
                velData.(planes{i}).vMean = velData.(planes{i}).vMean * (At / (At - Am));
                velData.(planes{i}).wMean = velData.(planes{i}).wMean * (At / (At - Am));
            end
            
        end
        
    case 'B'
        
        if contains(caseName, 'Run_Test') || (contains(caseName, 'Windsor') && contains(caseName, 'Upstream'))
            Am = (0.289 * 0.389) + (2 * (0.046 * 0.055));
            At = (2 * (0.9519083 + (3.283 * tan(atan(0.0262223 / 9.44)))) * 1.32);
            
            for i = 1:height(planes)
                velData.(planes{i}).uMean = velData.(planes{i}).uMean * (At / (At - Am));
                velData.(planes{i}).vMean = velData.(planes{i}).vMean * (At / (At - Am));
                velData.(planes{i}).wMean = velData.(planes{i}).wMean * (At / (At - Am));
                
                for j = 1:height(velData.(planes{i}).time)
                    velData.(planes{i}).u{j} = velData.(planes{i}).u{j} * (At / (At - Am));
                    velData.(planes{i}).v{j} = velData.(planes{i}).v{j} * (At / (At - Am));
                    velData.(planes{i}).w{j} = velData.(planes{i}).w{j} * (At / (At - Am));
                    velData.(planes{i}).uPrime{j} = velData.(planes{i}).uPrime{j} * (At / (At - Am));
                    velData.(planes{i}).vPrime{j} = velData.(planes{i}).vPrime{j} * (At / (At - Am));
                    velData.(planes{i}).wPrime{j} = velData.(planes{i}).wPrime{j} * (At / (At - Am));
                end
                
            end
            
        end
        
end

disp(' ');
disp(' ');


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
                nFrames = inputFrames(height(velData.(planes{1}).time));

                if nFrames == -1
                    continue;
                end

                valid = true;
            else
                disp('    WARNING: Invalid Entry');
            end

        end
        clear valid;

end

disp(' ');

for i = 1:height(plotPlanes)
    disp(['    Presenting ', plotPlanes{i}, '...']);
    
    orientation = velData.(plotPlanes{i}).planeOrientation;
    
    % Specify Default Axes Limits
    switch orientation

        case 'YZ'

            if contains(caseName, ["Run_Test", "Windsor", "Varney"])
                xLimsPlot = [0.31875; 4.65925];
                yLimsPlot = [-0.5945; 0.5945];
                zLimsPlot = [0; 0.739];
            end

        case {'XZ', 'XY'}

            if contains(caseName, ["Run_Test", "Windsor", "Varney"])
                xLimsPlot = [0.31875; 1.08325];
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
    switch format
        
        case 'A'
            xLimsData = xLimsPlot;
            yLimsData = yLimsPlot;
            zLimsData = zLimsPlot;
            
        case 'B'
            xLimsData = velData.(plotPlanes{i}).xLims;
            yLimsData = velData.(plotPlanes{i}).yLims;
            zLimsData = velData.(plotPlanes{i}).zLims;
            
            xLimsPlot = [min(min(xLimsPlot), min(xLimsData)); max(max(xLimsPlot), max(xLimsData))];
            yLimsPlot = [min(min(yLimsPlot), min(yLimsData)); max(max(yLimsPlot), max(yLimsData))];
            zLimsPlot = [min(min(zLimsPlot), min(zLimsData)); max(max(zLimsPlot), max(zLimsData))];
        
        case 'C'
            xLimsData = [min(velData.(planes{i}).position(:,1)); max(velData.(planes{i}).position(:,1))];
            yLimsData = [min(velData.(planes{i}).position(:,2)); max(velData.(planes{i}).position(:,2))];
            zLimsData = [min(velData.(planes{i}).position(:,3)); max(velData.(planes{i}).position(:,3))];
            
    end
    
    % Restrict Data Limits to Plane Position
    switch orientation
        
        case 'YZ'
            xLimsData = velData.(plotPlanes{i}).planePosition;
            
        case 'XZ'
            yLimsData = velData.(plotPlanes{i}).planePosition;
            
        case 'XY'
            zLimsData = velData.(plotPlanes{i}).planePosition;
            
    end
    
    positionData = velData.(plotPlanes{i}).position;
    vectorData = [velData.(plotPlanes{i}).uMean, velData.(plotPlanes{i}).vMean, velData.(plotPlanes{i}).wMean];
    
    if strcmp(caseName, 'Varney') && ~strcmp(orientation, 'YZ')
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
    
    fig = planarVectorPlots(orientation, xLimsData, yLimsData, zLimsData, ...
                            positionData, vectorData, nComponents, component, ...
                            fig, figName, cMap, geometry, streamlines, ...
                            xDims, yDims, zDims, figTitle, figSubtitle, cLims, ...
                            xLimsPlot, yLimsPlot, zLimsPlot, normalise);
end

disp(' ');

switch format

    case 'B'

        if plotInst
            for i = 1:height(plotPlanes)
                figHold = fig;
                
                for j = 1:nFrames
                    
                    if j ~= 1
                        clf(fig)
                        fig = figHold;
                    end

                    vectorData = [velData.(plotPlanes{i}).u{j}, velData.(plotPlanes{i}).v{j}, velData.(plotPlanes{i}).w{j}];
                    figTime = num2str(velData.(plotPlanes{i}).time(j), ['%.', num2str(timePrecision), 'f']);
                    figName = [caseName, '_', plotPlanes{i}, '_T', erase(figTime, '.')];
                    figSubtitle = [figTime, ' \it{s}'];
                    
                    fig = planarVectorPlots(orientation, xLimsData, yLimsData, zLimsData, ...
                                            positionData, vectorData, nComponents, component, ...
                                            fig, figName, cMap, geometry, streamlines, ...
                                            xDims, yDims, zDims, figTitle, figSubtitle, cLims, ...
                                            xLimsPlot, yLimsPlot, zLimsPlot, normalise);
                end
                
            end
            
        end
        
end


%% Local Functions

function nFrames = inputFrames(Nt)

    nFrames = str2double(input(['    Input Desired Frame Count [1-', num2str(Nt), ']: '], 's'));
    
    if isnan(nFrames) || nFrames <= 0 || nFrames > Nt
        disp('        WARNING: Invalid Entry');
        nFrames = -1;
    end

end