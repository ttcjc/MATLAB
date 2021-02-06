%% Particle Deposition Mapper v1.0

clearvars;
close all;
clc;

fig = 0; %#ok<*NASGU>
figHold = 0; %#ok<*NASGU>

disp ('================================');
disp ('Particle Deposition Mapper v1.0');
disp ('================================');
disp (' ');


%% Changelog

% v1.0 - Initial Commit (Base Contamination Only)


%% Case Selection

disp('CASE SELECTION');
disp('--------------');
disp(' ');

% caseFolder = uigetdir('~/OpenFOAM', 'Select Case');
% caseFolder = ('/home/cjc96/Mount/Uni/OpenFOAM/ttcjc-7/run/Test_Block');
caseFolder = ('/home/cjc96/Mount/Athena/OpenFOAM/ttcjc-7/results/Windsor_Square_MP_Test');

if contains(caseFolder, ["Windsor_Square", "Test_Block"])
    xDims = [-0.56075; 0.48325];
    yDims = [-(0.389 / 2); (0.389 / 2)];
    zDims = [0.05; 0.339];
else
    error('Unsupported Case')
end

[timeDirs, deltaT] = timeDirectories(caseFolder);
disp(' ');

disp(' ');
disp('GEOMETRY SELECTION');
disp('------------------');
disp(' ');

% [file, path] = uigetfile('~/CAD/CFD Geometries/*.*', 'Select Subject Geometry', 'multiSelect', 'on');
% file = 'Test_Block.stl';
% path = '/home/cjc96/CAD/CFD Geometries/';
file = {'Windsor_Axles.stl', 'Windsor_Body_Square.stl', 'Windsor_Wheels_Flat.stl'};
path = '/home/cjc96/CAD/CFD Geometries/Windsor Model/Output/Presentation Quality/';

if isa(file, 'cell')
    disp('Geometries: ');
    
    for i = 1:size(file,2)
        part = file{1,i}(1:end-4);
        geometry.(part) = stlread([path,file{1,i}]);
        disp(['    ', part]);
    end
    
else
    part = file(1:end-4);
    geometry.(part) = stlread([path, file]);
    disp(['Geometry: ', part]);
end

disp(' ');


%% Data Acquisition

disp(' ');
disp('LAGRANGIAN DATA COLLECTION');
disp('--------------------------');

valid = false;
while ~valid
    disp(' ');
    selection = input('Load Saved Particle Data? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        disp(' ');
        disp(' ');
        [particleData, particleProps] = lagrangianData_v2(caseFolder, timeDirs);
        valid = true;
        disp(' ');
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        [particleData, path] = uigetfile('~/Documents/Engineering/PhD/Data/Numerical/MATLAB/particleData/*.*', 'Select Lagrangian Data');
        
        if contains(path, '/MATLAB/particleData/')
            disp(['    Loading: ', particleData]);
            load([path, particleData])
            valid = true;
            disp(' ');
        else
            clear lagrangianData path;
            disp('    WARNING: Invalid File Selection');
            disp(' ');
        end
        
    else
        disp('    WARNING: Invalid Entry');
        disp(' ');
    end
    
end

%% Mapping Options

disp(' ');
disp('MAPPING OPTIONS');
disp('---------------');

valid = false;
while ~valid
    disp(' ');
    disp('Possible Mapping Locations:');
    disp('    A: Base');
    disp('    B: Body (NYI)');
    disp('    C: Far-Field Extraction Plane (NYI)');
    disp(' ');
    selection = input('    Select Method [A/B/C]: ', 's');

    if selection == 'a' | selection == 'A' %#ok<OR2>
        method = 'A';
        valid = true;
        disp(' ');
    elseif selection == 'b' | selection == 'B' %#ok<OR2>
        method = 'B';
        valid = true;
        disp(' ');
    elseif selection == 'c' | selection == 'C' %#ok<OR2>
        method = 'C';
        valid = true;
        disp(' ');
    else
        disp('        WARNING: Invalid Entry');
        disp(' ');
    end

end

% Temporary Implementation Control
if method == 'B' | method == 'C' %#ok<OR2>
    eror('Mapping Method Not Yet Implemented');
end

valid = false;
while ~valid
    disp(' ');
    selection = input('Filter Particle Diameters? [y/n]: ', 's');
    
    if selection == 'n' | selection == 'N' %#ok<OR2>
        minD = round(min(particleData.d{end,1}) * 1e6);
        maxD = round(max(particleData.d{end,1}) * 1e6);
        valid = true;
        disp(' ');
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        minD = inputD('Minimum');
        maxD = inputD('Maximum');
        
        if maxD < round(min(particleData.d{end,1}) * 1e6) || minD > round(max(particleData.d{end,1}) * 1e6)
            disp('        WARNING: No Lagrangian Data in Selected Data Range');
        elseif maxD < minD
            disp('        WARNING: Invalid Entry');
            disp(' ');
        else
            valid = true;
            disp(' ');
        end
        
    else
        disp('    WARNING: Invalid Entry');
        disp(' ');
    end

end

switch method
    
    case 'A'
        index = find(~particleData.active{end,1} & (particleData.d{end,1} * 1e6) >= minD & (particleData.d{end,1} * 1e6) <= maxD & round(particleData.globalPos{end,1}(:,1),5) == max(xDims));
        
	case 'B'
		index = find(~particleData.active{end,1} & (particleData.d{end,1} * 1e6) >= minD & (particleData.d{end,1} * 1e6) <= maxD);

	case 'C'
        % Will use contents of caseFolder/extractionData

end

for i = 1:size(particleProps,1)
    prop = particleProps{i,1};
    contaminantData.(prop) = particleData.(prop){end,1}(index,:);
end

disp(' ');

switch method
    
    case {'A', 'B'}
        disp([num2str(size(contaminantData.active,1)), ' Particles Recorded on Surface']);
		disp(' ');
        
    case 'C'
        disp([num2str(size(contaminantData.active,1)), ' Particles Recorded on Far-Field Extraction Plane']);
		disp(' ');
        
end


%% Mapping

disp(' ');
disp('***********');
disp('  Running  ');
disp(' ');

tic;

switch method
    
    case 'A'
        mappingData.pos = contaminantData.globalPos;
		cellSize = 0.01;
        
        for i = 1:size(contaminantData.active,1)
            mappingData.pos(i,2) = (round(contaminantData.globalPos(i,2) / cellSize)) * cellSize;
            mappingData.pos(i,3) = (round(contaminantData.globalPos(i,3) / cellSize)) * cellSize;
        end
        
        mappingData.occupiedCells = unique(mappingData.pos, 'stable', 'rows');
        
        mappingData.mass = zeros(size(mappingData.occupiedCells,1),1);

        
		for i = 1:size(mappingData.occupiedCells,1)
            index = find(ismember(mappingData.pos, mappingData.occupiedCells(i,:), 'rows'));
           
			for j = 1:size(index,1)
                mappingData.mass(i,1) = mappingData.mass(i,1) + (1000 * (4 / 3) * pi * (contaminantData.d(index(j,1),1) / 2) ^ 3);
			end
			
		end
        
		mappingData.massNorm = mappingData.mass / max(mappingData.mass);
		
        massInterp = scatteredInterpolant(mappingData.occupiedCells(:,2), mappingData.occupiedCells(:,3), mappingData.massNorm(:,1));
        [y, z] = meshgrid(min(yDims):0.001:max(yDims), min(zDims):0.001:max(zDims));
        mass = massInterp(y, z);
        
    % Figure Setup
    fig = fig + 1;
    figure('name', 'Base Contamination Map (CFD)');
    hold on;
    set(figure(fig), 'outerPosition', [25, 25, 800, 800]);

    % Plot
    contourf(y, z, mass, 16, 'edgeColor', 'none');

    % Figure Formatting
    tickData = min(yDims):((max(yDims) - min(yDims)) / 6):max(yDims);
    xticks(tickData(2:end-1));
    tickData = min(zDims):((max(zDims) - min(zDims)) / 6):max(zDims);
    yticks(tickData(2:end-1));
    xtickformat('%.3f');
    ytickformat('%.3f');
    caxis([0, 1]);
    xlabel({' ', 'y (\it{l})'});
    ylabel({'z (\it{l})', ' '});
    box on;
    colormap viridis;
    set(gca, 'units', 'normalized', 'position', [0.1275, 0.1275, 0.745, 0.745], ...
             'fontName', 'LM Roman 12', 'fontSize', 12, 'layer', 'top');
    hold off;

	namePos = max(strfind(caseFolder, '/'));
	savefig(fig, ['~/MATLAB/Figures/', caseFolder(namePos(end):end), '_Contamination_Map_Base_', num2str(minD), '_', num2str(maxD)]);
    print(fig, ['~/MATLAB/Figures/', caseFolder(namePos(end):end), '_Contamination_Map_Base_', num2str(minD), '_', num2str(maxD)], '-dpng', '-r300');
            
end

executionTime = toc;

disp(['    Execution Time: ', num2str(executionTime), 's']);
disp(' ');
disp('  Success  ');
disp('***********');


%% Colour Bar

% Figure Setup
fig = fig + 1;
figure('name', 'Colour Bar');
hold on;
set(figure(fig), 'outerPosition', [25, 75, 800, 800]);

% Figure Formatting
caxis([0, 1]);
axis off;
colormap viridis;
c = colorbar('ticks', (0.1:0.2:0.9), 'location', 'south', 'axisLocation', 'out');
c.Label.String = {' ', 'Normalised Contaminant Mass'};
set(gca, 'units', 'normalized', 'position', [0.1275, 0.1275, 0.745, 0.745], ...
		 'fontName', 'LM Roman 12', 'fontSize', 12, 'layer', 'top');
hold off;


namePos = max(strfind(caseFolder, '/'));
savefig(fig, ['~/MATLAB/Figures/', caseFolder(namePos(end):end), '_Contamination_Map_Mass_Bar_', num2str(minD), '_', num2str(maxD)]);
print(fig, ['~/MATLAB/Figures/', caseFolder(namePos(end):end), '_Contamination_Map_Mass_Bar_', num2str(minD), '_', num2str(maxD)], '-dpng', '-r300');

% clearvars -except particleData contaminantData mappingData;


%% Local Functions
    
function D = inputD(type)

    valid = false;
    while ~valid
        disp(' ');
        D = str2double(input(['    ', type, ' Diameter of Interest [', char(956), 'm]: '], 's'));

        if isnan(D) || length(D) > 1
            disp('        WARNING: Invalid Entry');
        else
            valid = true;
        end

    end

end