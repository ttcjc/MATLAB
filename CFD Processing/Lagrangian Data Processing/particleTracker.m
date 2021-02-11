%% Particle Tracker v1.0

clearvars;
close all;
clc;

fig = 0; %#ok<*NASGU>
figHold = 0; %#ok<*NASGU>

disp ('=====================');
disp ('Particle Tracker v1.0');
disp ('=====================');
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


%% Tracking Options

disp('TRACKING OPTIONS');
disp('----------------');

disp(' ');
disp('Possible Tracking Methods:');
disp('    A: Self-Contamination');
disp('    B: Third-Party Contamination');
disp('    C: Far-Field Particle Transport');
disp(' ');

valid = false;
while ~valid
	disp(' ');
	selection = input('Select Tracking Method [A/B/C]: ', 's');

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
		disp('    WARNING: Invalid Entry');
		disp(' ');
	end

end

valid = false;
while ~valid
	disp(' ');
	selection = input('Enable Position Interpolation? [y/n]: ', 's');

	if selection == 'n' | selection == 'N' %#ok<OR2>
		interp = false;
		valid = true;
		disp(' ');
	elseif selection == 'y' | selection == 'Y' %#ok<OR2>
		interp = true;          
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

% Identify Particles of Interest
switch method

	case {'A', 'B'}
		index = find(~particleData.active{end,1} & (particleData.d{end,1} * 1e6) >= minD & (particleData.d{end,1} * 1e6) <= maxD);

	case 'C'
		index = find(particleData.active{end,1} & (particleData.d{end,1} * 1e6) >= minD & (particleData.d{end,1} * 1e6) <= maxD & particleData.globalPos{end,1}(:,1) > 1);

end

% Collate Particles of Interest
for i = 1:size(particleProps,1)
	prop = particleProps{i,1};
	contaminantData.(prop) = particleData.(prop){end,1}(index,:);
end

disp(' ');

switch method

	case {'A', 'B'}
		disp([num2str(size(contaminantData.active,1)), ' Particles Recorded on Surface']);

	case 'C'
		disp([num2str(size(contaminantData.active,1)), ' Particles Recorded in the Far-Field']);

end

if size(contaminantData.active,1) > 1000
	disp('    WARNING: Tracking Large Datasets Will Be Computationally Expensive and Difficult to Visualise');
end

disp(' ');

valid = false;
while ~valid
	disp(' ');
	selection = input('Limit Tracking Count? [y/n]: ', 's');

	if selection == 'n' | selection == 'N' %#ok<OR2>
		count = size(contaminantData.active,1);
		valid = true;
		disp(' ');
	elseif selection == 'y' | selection == 'Y' %#ok<OR2>
		count = inputCount();

		if count > size(contaminantData.active,1)
	disp('    WARNING: Value Exceeds Available Particle Count');
	else
		valid = true;
		disp(' ');
	end

	else
		disp('        WARNING: Invalid Entry');
		disp(' ');
end

end

% Select a Random Set of Particles to Track
if count ~= size(contaminantData.active,1)
	index = sort(randperm(size(contaminantData.active,1),count))';

	for i = 1:size(particleProps,1)
		prop = particleProps{i,1};
		contaminantData.(prop) = contaminantData.(prop)(index,:);
	end

end

disp(' ');


%% Tracking

tic;
disp('***********');
disp('  Running  ');
disp(' ');
disp(['    Tracking ', num2str(size(contaminantData.active,1)), ' Particles']);

trackingData = cell(size(contaminantData.origId,1), size(particleData.origId,1));

for i = 1:size(particleData.time,1)
	trackingID = intersect(horzcat(particleData.origId{i,1}, particleData.origProcId{i,1}), horzcat(contaminantData.origId, contaminantData.origProcId), 'rows', 'stable');
	index = find(ismember(horzcat(particleData.origId{i,1}, particleData.origProcId{i,1}), trackingID, 'rows'));

	for j = 1:size(trackingID,1)
		trackingData{j,i} = particleData.globalPos{i,1}(index(j),:);
	end

end

switch method

	case 'A'
		% Figure Setup
		fig = fig + 1;
		figure('name', 'Tracking: Self-Contamination');
		hold on;
		set(figure(fig), 'outerPosition', [25, 25, 800, 800]);
		caxis([minD, maxD]);
		colormap viridis;
		map = colormap;

		% Plot
		parts = fieldnames(geometry);

		for i = 1:size(parts,1)
			part = parts{i,1};
			patch(geometry.(part), 'faceColor', [0.5, 0.5, 0.5], 'edgeColor', [0.5, 0.5, 0.5]);
		end

		for i = 1:size(trackingData,1)
			particlePath = unique(cell2mat(trackingData(i,:)'), 'stable', 'rows');

			if size(particlePath,1) == 1
				continue
			elseif interp && size(particlePath,1) > 2
				particlePath = interparc((round(deltaT / 0.001) * size(particlePath,1)), particlePath(:,1), particlePath(:,2), particlePath(:,3), 'spline');
			end

			x = particlePath(:,1);
			y = particlePath(:,2);
			z = particlePath(:,3);

			index = round(1 + (size(map,1) - 1) * ((contaminantData.d(i) * 1e6) - minD) / (maxD - minD));

			plot3(x, y, z, 'color', map(index,:));
			scatter3(x(end), y(end), z(end), 10, (contaminantData.d(i) * 1e6), 'filled');
		end

		% Figure Formatting
		view([45, 15]);
		axis off
		set(gca, 'units', 'normalized', 'position', [0.1275, 0.1275, 0.745, 0.745], ...
				 'fontName', 'LM Roman 12', 'fontSize', 12, 'layer', 'top', ...
				 'dataAspectRatio', [1, 1, 1]);
		hold off;
		
		% Save Plot
		namePos = max(strfind(caseFolder, '/'));
		savefig(fig, ['~/MATLAB/Output/Figures/', caseFolder(namePos(end):end), '_Tracking_SC_D', num2str(minD), '_D', num2str(maxD), '_c', num2str(count)]);
		print(fig, ['~/MATLAB/Output/Figures/', caseFolder(namePos(end):end), '_Tracking_SC_D', num2str(minD), '_D', num2str(maxD), '_c', num2str(count)], '-dpng', '-r300');

	case 'B'
		% Figure Setup
		fig = fig + 1;
		figure('name', 'Tracking: Third-Party Contamination');
		hold on;
		set(figure(fig), 'outerPosition', [25, 25, 800, 800]);
		caxis([minD, maxD]);
		colormap viridis;
		map = colormap;

		% Plot
		parts = fieldnames(geometry);

		for i = 1:size(parts,1)
			part = parts{i,1};
			patch(geometry.(part), 'faceColor', [0.5, 0.5, 0.5], 'edgeColor', [0.5, 0.5, 0.5]);
		end

		for i = 1:size(trackingData,1)
			particlePath = unique(cell2mat(trackingData(i,:)'), 'stable', 'rows');

			if size(particlePath,1) == 1
				continue
			elseif interp && size(particlePath,1) > 2
				particlePath = interparc((round(deltaT / 0.001) * size(particlePath,1)), particlePath(:,1), particlePath(:,2), particlePath(:,3), 'spline');
			end

			x = particlePath(:,1);
			y = particlePath(:,2);
			z = particlePath(:,3);

			index = round(1 + (size(map,1) - 1) * ((contaminantData.d(i) * 1e6) - minD) / (maxD - minD));

			plot3(x, y, z, 'color', map(index,:));
			scatter3(x(end), y(end), z(end), 10, (contaminantData.d(i) * 1e6), 'filled');
		end

		% Figure Formatting
		view([-45, 15]);
		axis off
		set(gca, 'units', 'normalized', 'position', [0.1275, 0.1275, 0.745, 0.745], ...
				 'fontName', 'LM Roman 12', 'fontSize', 12, 'layer', 'top', ...
				 'dataAspectRatio', [1, 1, 1]);
		hold off;

		% Save Plot
		namePos = max(strfind(caseFolder, '/'));
		savefig(fig, ['~/MATLAB/Output/Figures/', caseFolder(namePos(end):end), '_Tracking_TP_D', num2str(minD), '_D', num2str(maxD), '_c', num2str(count)]);
		print(fig, ['~/MATLAB/Output/Figures/', caseFolder(namePos(end):end), '_Tracking_TP_D', num2str(minD), '_D', num2str(maxD), '_c', num2str(count)], '-dpng', '-r300');

	case 'C'
		% Figure Setup
		fig = fig + 1;
		figure('name', 'Tracking: Far-Field Contaminants');
		hold on;
		set(figure(fig), 'outerPosition', [25, 25, 800, 800]);
		caxis([minD, maxD]);
		colormap viridis;
		map = colormap;

		% Plot
		parts = fieldnames(geometry);

		for i = 1:size(parts,1)
			part = parts{i,1};
			patch(geometry.(part), 'faceColor', [0.5, 0.5, 0.5], 'edgeColor', [0.5, 0.5, 0.5]);
		end

		for i = 1:size(trackingData,1)
			particlePath = unique(cell2mat(trackingData(i,:)'), 'stable', 'rows');

			if size(particlePath,1) == 1
				continue
			elseif interp && size(particlePath,1) > 2
				particlePath = interparc((round(deltaT / 0.001) * size(particlePath,1)), particlePath(:,1), particlePath(:,2), particlePath(:,3), 'spline');
			end

			x = particlePath(:,1);
			y = particlePath(:,2);
			z = particlePath(:,3);

			index = round(1 + (size(map,1) - 1) * ((contaminantData.d(i) * 1e6) - minD) / (maxD - minD));

			plot3(x, y, z, 'color', map(index,:));
			scatter3(x(end), y(end), z(end), 10, (contaminantData.d(i) * 1e6), 'filled');
		end

		% Figure Formatting
		view([45, 15]);
		axis off
		set(gca, 'units', 'normalized', 'position', [0.1275, 0.1275, 0.745, 0.745], ...
				 'fontName', 'LM Roman 12', 'fontSize', 12, 'layer', 'top', ...
				 'dataAspectRatio', [1, 1, 1]);
		hold off;

		% Save Plot
		namePos = max(strfind(caseFolder, '/'));
		savefig(fig, ['~/MATLAB/Output/Figures/', caseFolder(namePos(end):end), '_Tracking_FF_D', num2str(minD), '_D', num2str(maxD), '_c', num2str(count)]);
		print(fig, ['~/MATLAB/Output/Figures/', caseFolder(namePos(end):end), '_Tracking_FF_D', num2str(minD), '_D', num2str(maxD), '_c', num2str(count)], '-dpng', '-r300');

end

executionTime = toc;

disp(' ');
disp(['    Tracking Time: ', num2str(executionTime), 's']);
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
c.Label.String = {' ', 'Particle Diameter (\it{\mum})'};
set(gca, 'units', 'normalized', 'position', [0.1275, 0.1275, 0.745, 0.745], ...
		 'fontName', 'LM Roman 12', 'fontSize', 12, 'layer', 'top');
hold off;

% Save Plot
namePos = max(strfind(caseFolder, '/'));
savefig(fig, ['~/MATLAB/Output/Figures/', caseFolder(namePos(end):end), '_Tracking_Diameter_Bar_', num2str(minD), '_', num2str(maxD)]);
print(fig, ['~/MATLAB/Output/Figures/', caseFolder(namePos(end):end), '_Tracking_Diameter_Bar_', num2str(minD), '_', num2str(maxD)], '-dpng', '-r300');


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

disp('trackingData:');
disp('    Cartesian Positions of Tracked Particles');


%% Cleaning

clearvars -except particleData contaminantData trackingData;
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

function C = inputCount()

	valid = false;
	while ~valid
		disp(' ');
		C = str2double(input('        Number of Particles to Track: ', 's'));

		if isnan(C) || length(C) > 1
			disp('            WARNING: Invalid Entry');
		else
			valid = true;
		end

	end

end
