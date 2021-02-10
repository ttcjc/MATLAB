%% Case Initialisation v1.0


%% Changelog

% v1.0 - Initial Commit


%% Case Support

% Test_Block, Windsor_Square


%% Main Function

function[caseFolder, xDims, yDims, zDims, timeDirs, deltaT, geometry] = initialiseCase

	disp('CASE SELECTION');
	disp('--------------');
	disp(' ');

	caseFolder = uigetdir('~/OpenFOAM', 'Select Case');

	disp(['Case: ', caseFolder]);
	disp(' ');
	disp(' ');

	% Confirm Support
	if ~contains(caseFolder, ["Lag_Test", "Test_Block", "Windsor_Square"])
		error('Unsupported Case');
	end
	
	% Set Dimensions of Interest
	if contains(caseFolder, "Lag_Test")
		xDims = [-14; 14];
		yDims = [-0.96; 0.96];
		zDims = [0; 1.32];
	elseif contains(caseFolder, "Test_Block")
		xDims = [-0.56075; 0.48325];
		yDims = [-(0.389 / 2); (0.389 / 2)];
		zDims = [0.05; 0.339];
	elseif contains(caseFolder, "Windsor_Square")
		xDims = [-0.56075; 0.48325];
		yDims = [-(0.389 / 2); (0.389 / 2)];
		zDims = [0.05; 0.339];
	end
	
	% Confirm Valid Case
	disp('Analysing Case Structure:');
	
	if exist([caseFolder, '/0'], 'dir')
		timeDirs = dir(caseFolder);
	else
		error('Invalid Case Directory');
	end

	% Identify Time Directories
	i = 1;
	while i <= size(timeDirs,1)

		if isnan(str2double(timeDirs(i,1).name))
			timeDirs(i,:) = [];
		else
			i = i + 1;
		end

	end

	deltaT = str2double(timeDirs(end,1).name) - str2double(timeDirs(end-1,1).name);

	disp(['    Identified ', num2str(size(timeDirs,1)), ' Time Directories']);
	disp(['    ', char(916), 'T = ' num2str(deltaT), 's']);
	disp(' ');
	disp(' ');
	
	disp('GEOMETRY SELECTION');
	disp('------------------');
	disp(' ');

	if contains(caseFolder, "Lag_Test")
		disp('No Test Geometry Required');
		geometry = 0;
	else
		[file, path] = uigetfile('~/CAD/CFD Geometries/*.*', 'Select Subject Geometry', 'multiSelect', 'on');

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
	end

	disp(' ');
	disp(' ');
	
end
