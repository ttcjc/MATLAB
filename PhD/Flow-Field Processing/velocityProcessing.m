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
        [caseFolder, data, geometry, xDims, yDims, zDims, precision] = initialiseProbeData('U', true, nProc);

    case 'C'
        [data, geometry, xDims, yDims, zDims, precision] = initialiseExpData('U', true);

end


%% Formatting

switch format

    case 'A'
        planes = fieldnames(data);
        
        % Adjust Data Origin
        if contains(caseFolder, ["Run_Test", "Windsor"]) && contains(caseFolder, 'Upstream')
            
            for i = 1:height(planes)
                data.(planes{i}).position(:,1) = data.(planes{i}).position(:,1) + 1.325;
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
        % Adjust Data Origin
        if contains(caseFolder, ["Run_Test", "Windsor"]) && contains(caseFolder, 'Upstream')
            data.position(:,1) = data.position(:,1) + 1.325;
        end
        
        % Normalise Velocity
        if contains(caseFolder, ["Run_Test", "Windsor"])
            data.uMean = data.uMean / 40;
            data.vMean = data.vMean / 40;
            data.wMean = data.wMean / 40;

            for i = 1:height(data.time)
                data.u{i} = data.u{i} / 40;
                data.v{i} = data.v{i} / 40;
                data.w{i} = data.w{i} / 40;
                data.uPrime{i} = data.uPrime{i} / 40;
                data.vPrime{i} = data.vPrime{i} / 40;
                data.wPrime{i} = data.wPrime{i} / 40;
            end
            
        end
    
    case 'C'
        planes = fieldnames(data);
        
        % Normalise Velocity
        for i = 1:height(planes)
                data.(planes{i}).uMean = data.(planes{i}).uMean / 40;
                data.(planes{i}).vMean = data.(planes{i}).vMean / 40;
                data.(planes{i}).wMean = data.(planes{i}).wMean / 40;
                data.(planes{i}).uRMS = data.(planes{i}).uRMS / 40;
                data.(planes{i}).vRMS = data.(planes{i}).vRMS / 40;
                data.(planes{i}).wRMS = data.(planes{i}).wRMS / 40;
        end
        
end

disp (' ');
disp (' ');


%% Presentation

disp('Data Presentation');
disp('------------------');