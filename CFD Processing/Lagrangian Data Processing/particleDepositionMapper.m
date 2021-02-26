%% Particle Deposition Mapper v1.0

clearvars;
close all;
clc;

fig = 0; %#ok<*NASGU>
figHold = 0; %#ok<*NASGU>

disp ('===============================');
disp ('Particle Deposition Mapper v1.0');
disp ('===============================');
disp (' ');


%% Changelog

% v1.0 - Initial Commit (Base Contamination Only)


%% Case Initialisation

[caseFolder, xDims, yDims, zDims, timeDirs, deltaT, geometry] = initialiseCase;


%% Lagrangian Data Acquisition

disp('LAGRANGIAN DATA ACQUISITION');
disp('---------------------------');

valid = false;
while ~valid
	disp(' ');
	selection = input('Load Saved Particle Data? [y/n]: ', 's');

	if selection == 'n' | selection == 'N' %#ok<OR2>
		disp(' ');
		disp(' ');
		[particleData, particleProps] = lagrangianData(caseFolder, timeDirs);
		valid = true;
		disp(' ');
	elseif selection == 'y' | selection == 'Y' %#ok<OR2>
		[particleData, path] = uigetfile('~/Documents/Engineering/PhD/Data/Numerical/MATLAB/particleData/*.*', 'Select Lagrangian Data');

		if contains(path, '/MATLAB/particleData/')
			disp(' ');
			disp(' ');
			disp(['Loading: ', particleData]);
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

disp(' ');


%% Mapping Options

disp('MAPPING OPTIONS');
disp('---------------');

disp(' ');
disp('Possible Mapping Locations:');
disp('    A: Base');
disp('    B: Body (NYI)');
disp('    C: Far-Field Extraction Plane (NYI)');

valid = false;
while ~valid
	disp(' ');
	selection = input('    Select Mapping Location [A/B/C]: ', 's');

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
if method ~= 'A'
	error('Mapping Location Not Yet Implemented');
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

disp(' ');
disp('Data Available for the Following Time Instances:');

for i = 1:size(particleData.time,1)
	disp(['    ', num2str(particleData.time{i,1}), 's']);
end

valid = false;
while ~valid
	disp(' ');
	selection = input('Plot Only Specific Time Instances? [y/n]: ', 's');

	if selection == 'n' | selection == 'N' %#ok<OR2>
		times = 1:size(particleData.time,1);
		valid = true;
		disp(' ');
	elseif selection == 'y' | selection == 'Y' %#ok<OR2>
		times = inputTimes(particleData)';

		if isempty(times)
			disp('        WARNING: No Lagrangian Data in Selected Data Range');
		else
			valid = true;
			disp(' ');
		end

	else
		disp('    WARNING: Invalid Entry');
		disp(' ');
	end

end

disp(' ');


%% Mapping

tic;
disp('***********');
disp('  Running  ');
disp(' ');

for i = 1:size(times,2)

	% Identify Particles of Interest
	switch method

		case 'A'
			index = find(~particleData.active{times(1,i),1} & (particleData.d{times(1,i),1} * 1e6) >= minD & (particleData.d{times(1,i),1} * 1e6) <= maxD & ...
					round(particleData.positionCartesian{times(1,i),1}(:,1),5) == max(xDims));

		case 'B'
			index = find(~particleData.active{times(1,i),1} & (particleData.d{times(1,i),1} * 1e6) >= minD & (particleData.d{times(1,i),1} * 1e6) <= maxD & ...
					round(particleData.positionCartesian{times(1,i),1}(:,1),5) >= min(xDims) & round(particleData.positionCartesian{times(1,i),1}(:,1),5) <= max(xDims) & ...
					round(particleData.positionCartesian{times(1,i),1}(:,2),3) >= min(yDims) & round(particleData.positionCartesian{times(1,i),1}(:,2),3) <= max(yDims) & ...
					round(particleData.positionCartesian{times(1,i),1}(:,3),3) >= min(zDims) & round(particleData.positionCartesian{times(1,i),1}(:,3),3) <= max(zDims));

		case 'C'
			% Will use caseFolder/extractionData

	end

	% Collate Particles of Interest
	for j = 1:size(particleProps,1)
		prop = particleProps{j,1};
		contaminantData.(prop){i,1} = particleData.(prop){times(1,i),1}(index,:);
	end
	
	if isempty(contaminantData.active{i,1})
		disp(['    Time = ', num2str(particleData.time{times(i),1}), 's: No Particles Recorded on Surface']);
		continue
	end
	
	switch method

		case 'A'
			disp(['    Time = ', num2str(particleData.time{times(i),1}), 's: ', num2str(size(contaminantData.active{i,1},1)), ' Particles Recorded on Surface']);
			
			% Initialise Contaminant Map
			cellSize = 0.01;
			mappingData.pos{i,1} = contaminantData.positionCartesian{i,1};
		
			% Assign Contaminant Data to Map Nodes
			for j = 1:size(contaminantData.active{i,1},1)
				mappingData.pos{i,1}(j,2) = (round(contaminantData.positionCartesian{i,1}(j,2) / cellSize)) * cellSize;
				mappingData.pos{i,1}(j,3) = (round(contaminantData.positionCartesian{i,1}(j,3) / cellSize)) * cellSize;
			end

			mappingData.occupiedCells{i,1} = unique(mappingData.pos{i,1}, 'stable', 'rows');

			% Calculate Contaminant Mass in Each Mesh Node
			mappingData.mass{i,1} = zeros(size(mappingData.occupiedCells{i,1},1),1);

			for j = 1:size(mappingData.occupiedCells{i,1},1)
				index = find(ismember(mappingData.pos{i,1}, mappingData.occupiedCells{i,1}(j,:), 'rows'));

				for k = 1:size(index,1)
					mappingData.mass{i,1}(j,1) = mappingData.mass{i,1}(j,1) + (1000 * (4 / 3) * pi * (contaminantData.d{i,1}(index(k,1),1) / 2)^3);
				end

			end

			% Normalise Mass Values
			mappingData.massNorm{i,1} = mappingData.mass{i,1} / max(mappingData.mass{i,1});

			massInterp = scatteredInterpolant(mappingData.occupiedCells{i,1}(:,2), mappingData.occupiedCells{i,1}(:,3), mappingData.massNorm{i,1}(:,1));
			[y, z] = meshgrid(min(yDims):0.001:max(yDims), min(zDims):0.001:max(zDims));
			mass = massInterp(y, z);

			% Figure Setup
			fig = fig + 1;
			figure('name', ['Base Contamination Map (', num2str(particleData.time{times(i),1}), 's)']);
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

			% Save Plot
			namePos = max(strfind(caseFolder, '/'));
			timeInstance = erase(num2str(particleData.time{times(i),1}), '.');
			savefig(fig, ['~/MATLAB/Output/Figures/', caseFolder(namePos(end):end), '_Contamination_Map_Base_D', num2str(minD), '_D', num2str(maxD), '_', timeInstance, 's']);
			print(fig, ['~/MATLAB/Output/Figures/', caseFolder(namePos(end):end), '_Contamination_Map_Base_D', num2str(minD), '_D', num2str(maxD), '_', timeInstance, 's'], '-dpng', '-r300');

	end
	
end

executionTime = toc;

disp(' ');
disp(['    Execution Time: ', num2str(executionTime), 's']);
disp(' ');
disp('  Success  ');
disp('***********');

disp(' ');
disp(' ');


%% Data Key

disp('LAGRANGIAN DATA KEY');
disp('-------------------');
disp(' ');

disp('ParticleData:');
disp('    All Recorded Particles');
disp(' ');

disp('contaminantData:');
disp('    Particles of Interest');
disp(' ');

disp('mappingData:');
disp('    Surface Map of contaminantData');


%% Cleaning

clearvars -except particleData contaminantData mappingData;
disp(' ');


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

function T = inputTimes(particleData)

	valid = false;
	while ~valid
		disp(' ');
		T = str2num(input('    List Desired Time Instances (Row Vector Form): ', 's')); %#ok<ST2NM>

		if any(isnan(T)) || ~isrow(T)
			disp('        WARNING: Invalid Entry');
		else
			valid = true;
		end

	end

	T = find(ismember(particleData.time, strsplit(num2str(T))));

end
