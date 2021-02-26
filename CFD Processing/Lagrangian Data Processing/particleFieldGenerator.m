%% Particle Volume Field Generator v1.0

clearvars;
close all;
clc;

fig = 0; %#ok<*NASGU>
figHold = 0; %#ok<*NASGU>

disp ('====================================');
disp ('Particle Volume Field Generator v1.0');
disp ('====================================');
disp (' ');


%% Changelog

% v1.0 - Initial Commit


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


%% Field Options

disp('FIELD OPTIONS');
disp('-------------');

disp(' ');
disp('Possible Field Locations:');
disp('    A: Downstream (Self-Contamination)');
disp('    B: Upstream (Third-Party Contamination)');
disp(' ');

valid = false;
while ~valid
	disp(' ');
	selection = input('Select Field Location [A/B]: ', 's');

	if selection == 'a' | selection == 'A' %#ok<OR2>
		method = 'A';
		valid = true;
		disp(' ');
	elseif selection == 'b' | selection == 'B' %#ok<OR2>
		method = 'B';
		valid = true;
		disp(' ');
	else
		disp('    WARNING: Invalid Entry');
		disp(' ');
	end

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


%% Field Generation

tic;
disp('***********');
disp('  Running  ');
disp(' ');


% Volume of Interest
switch method

	case 'A'
		xRange = [0.3; 2.3];
		yRange = [-0.48; 0.48];
		zRange = [0; 0.66];

	case 'B'
		xRange = [-3.68325; 3.68325];
		yRange = [-0.72; 0.72];
		zRange = [0; 0.99];

end


for i = 1:size(times,2)

	% Identify Particles of Interest
	switch method

		case 'A'
			index = find(particleData.active{times(1,i),1} & ...
			(particleData.d{times(1,i),1} * 1e6) >= minD & (particleData.d{times(1,i),1} * 1e6) <= maxD & ...
			particleData.positionCartesian{times(1,i),1}(:,1) >= min(xRange) & particleData.positionCartesian{times(1,i),1}(:,1) <= max(xRange) & ...
			particleData.positionCartesian{times(1,i),1}(:,2) >= min(yRange) & particleData.positionCartesian{times(1,i),1}(:,2) <= max(yRange) & ...
			particleData.positionCartesian{times(1,i),1}(:,3) >= min(zRange) & particleData.positionCartesian{times(1,i),1}(:,3) <= max(zRange));

		case 'B'
			index = find(particleData.active{times(1,i),1} & ...
			(particleData.d{times(1,i),1} * 1e6) >= minD & (particleData.d{times(1,i),1} * 1e6) <= maxD & ...
			particleData.positionCartesian{times(1,i),1}(:,1) >= min(xRange) & particleData.positionCartesian{times(1,i),1}(:,1) <= max(xRange) & ...
			particleData.positionCartesian{times(1,i),1}(:,2) >= min(yRange) & particleData.positionCartesian{times(1,i),1}(:,2) <= max(yRange) & ...
			particleData.positionCartesian{times(1,i),1}(:,3) >= min(zRange) & particleData.positionCartesian{times(1,i),1}(:,3) <= max(zRange));

	end
	
	% Collate Particles of Interest
	for j = 1:size(particleProps,1)
		prop = particleProps{j,1};
		contaminantData.(prop){i,1} = particleData.(prop){times(1,i),1}(index,:);
	end

	% Initialise Volume Field
	cellSize = 0.005;
	[fieldData.y, fieldData.x, fieldData.z] = meshgrid(min(yRange):cellSize:max(yRange), min(xRange):cellSize:max(xRange), min(zRange):cellSize:max(zRange));
	
	fieldData.nParticles{i,1} = zeros(size(fieldData.x));
	fieldData.volumeFraction{i,1} = zeros(size(fieldData.x));
	fieldData.dMean{i,1} = zeros(size(fieldData.x));

	% Assign Contaminant Data to Field Nodes
	for j = 1:size(contaminantData.active{i,1},1)
		xDist = abs(contaminantData.positionCartesian{i,1}(j,1) - fieldData.x(:,1,1));
		yDist = abs(contaminantData.positionCartesian{i,1}(j,2) - fieldData.y(1,:,1));
		zDist = abs(contaminantData.positionCartesian{i,1}(j,3) - fieldData.z(1,1,:));

		index = zeros(1,3); % Array Position of Closest Mesh Node
		[~, index(1)] = min(xDist);
		[~, index(2)] = min(yDist);
		[~, index(3)] = min(zDist);

		fieldData.nParticles{i,1}(index(1),index(2),index(3)) = fieldData.nParticles{i,1}(index(1),index(2),index(3)) + contaminantData.nParticle{i,1}(j,1);
		fieldData.volumeFraction{i,1}(index(1),index(2),index(3)) = fieldData.volumeFraction{i,1}(index(1),index(2),index(3)) + (contaminantData.nParticle{i,1}(j,1) * (4 / 3) * pi * (contaminantData.d{i,1}(j,1) / 2)^3);
		fieldData.dMean{i,1}(index(1),index(2),index(3)) = fieldData.dMean{i,1}(index(1),index(2),index(3)) + contaminantData.d{i,1}(j,1);
	end

	fieldData.volumeFraction{i,1} = smooth3((fieldData.volumeFraction{i,1} / (cellSize ^ 3)), 'box');
	fieldData.dMean{i,1} = fieldData.dMean{i,1} ./ max(fieldData.nParticles{i,1},1);
% 	fieldData.dMean{i,1} = ((6 / pi) * (fieldData.volumeFraction{i,1} ./ max(fieldData.nParticles{i,1},1))).^(1 / 3);

	% Figure Setup
	fig = fig + 1;
	figure('name', ['Lagrangian Volume Field (', num2str(particleData.time{times(i),1}), 's)']);
	hold on;
	set(figure(fig), 'outerPosition', [25, 25, 800, 800]);
	caxis([minD, maxD]);
	colormap viridis;

	% Plot
	parts = fieldnames(geometry);

	for j = 1:size(parts,1)
		part = parts{j,1};
		patch(geometry.(part), 'faceColor', [0.5, 0.5, 0.5], 'edgeColor', [0.5, 0.5, 0.5]);
	end

	isoValue = 1e-6; % Desired Volume Fraction of Isosurface
	
	[faces, vertices, colours] = isosurface(fieldData.x, fieldData.y, fieldData.z, fieldData.volumeFraction{i,1}, isoValue, (fieldData.dMean{i,1} * 1e6));
	p = patch('vertices', vertices, 'faces', faces, 'faceVertexCData', colours, 'faceColor', 'interp', 'edgeColor', 'interp');

	% Figure Formatting
	switch method
		
		case 'A'
			xlim([0.31875, 2]);
			
		case 'B'
			xlim([-2, -0.31875]);
			
	end
	
	ylim([-0.4, 0.4]);
	zlim([0, 0.5445]);
	view([0, 0]);
	axis off
	set(gca, 'units', 'normalized', 'position', [0.1275, 0.1275, 0.745, 0.745], ...
			 'fontName', 'LM Roman 12', 'fontSize', 12, 'layer', 'top', ...
			 'dataAspectRatio', [1, 1, 1]);
	hold off;

	% Save Plot
	namePos = max(strfind(caseFolder, '/'));
	timeInstance = erase(num2str(particleData.time{times(i),1}), '.');
	savefig(fig, ['~/MATLAB/Output/Figures/', caseFolder(namePos(end):end), '_Particle_Field_D', num2str(minD), '_D', num2str(maxD), '_', timeInstance, 's']);
	print(fig, ['~/MATLAB/Output/Figures/', caseFolder(namePos(end):end), '_Particle_Field_D', num2str(minD), '_D', num2str(maxD), '_', timeInstance, 's'], '-dpng', '-r300');

end

executionTime = toc;

disp(['    Execution Time: ', num2str(executionTime), 's']);
disp(' ');
disp('  Success  ');
disp('***********');

disp(' ');
disp(' ');


%% Colour Bar

% Figure Setup
fig = fig + 1;
figure('name', 'Colour Bar');
hold on;
set(figure(fig), 'outerPosition', [25, 75, 800, 800]);

% Figure Formatting
caxis([minD, maxD]);
tickData = min(caxis):((max(caxis) - min(caxis)) / 6):max(caxis);
axis off;
colormap viridis;
c = colorbar('ticks', tickData(2:end-1), 'location', 'south', 'axisLocation', 'out');
c.Label.String = {' ', 'Mean Particle Diameter (\it{\mum})'};
set(gca, 'units', 'normalized', 'position', [0.1275, 0.1275, 0.745, 0.745], ...
		 'fontName', 'LM Roman 12', 'fontSize', 12, 'layer', 'top');
hold off;

% Save Plot
namePos = max(strfind(caseFolder, '/'));
savefig(fig, ['~/MATLAB/Output/Figures/', caseFolder(namePos(end):end), '_Field_Diameter_Bar_D', num2str(minD), '_D', num2str(maxD)]);
print(fig, ['~/MATLAB/Output/Figures/', caseFolder(namePos(end):end), '_Field_Diameter_Bar_D', num2str(minD), '_D', num2str(maxD)], '-dpng', '-r300');


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

disp('fieldData:');
disp('    Volumetric Interpretation of contaminantData');


%% Cleaning

clearvars -except particleData contaminantData fieldData;
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
