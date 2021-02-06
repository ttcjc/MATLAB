%% Third Party Case Initialisation v2.0

clearvars;
close all;
clc;

disp ('====================================');
disp ('THIRD-PARTY CASE INITIALISATION v2.0');
disp ('====================================');
disp (' ');


%% Changelog

% v1.0 - Initial Commit
% v2.0 - Added Inlet Flow Initialisation
% v2.1 


%% Case Selection

disp('CASE SELECTION');
disp('--------------');
disp(' ');

% caseFolderExtract = uigetdir('~/OpenFOAM', 'Extraction Case');
caseFolderExtract = '~/Mount/Uni/OpenFOAM/ttcjc-7/run/Test_Block';

% caseFolderInject = uigetdir('~/OpenFOAM', 'Injection Case');
caseFolderInject = '~/Mount/Uni/OpenFOAM/ttcjc-7/run/Test_Block_TP';

% Set Dimensions of Interest
wallOffset = 10e-6; % Offset to Prevent Wall Interactions
xDims = [(-9.408 + wallOffset); (9.408 - wallOffset)];
yDims = [(-0.96 + wallOffset); (0.96 - wallOffset)];
zDims = [(0 + wallOffset); (1.32 - wallOffset)];

% Confirm Valid Case Directories
disp(' ');
disp('Analysing Case Structures:');

if ~exist([caseFolderExtract, '/0'], 'dir') && ~exist([caseFolderExtract, '/extractionData/extractionData'], 'file') && ~exist([caseFolderExtract, '/postProcessing/planeSample'], 'dir')
	error([caseFolderExtract, ' Invalid']);
else
	disp(['    ', caseFolderExtract, ' Valid']);
end

if ~exist([caseFolderInject, '/0'], 'dir')
	error([caseFolderInject, 'Invalid']);
else
	disp(['    ', caseFolderInject, ' Valid']);
end

disp(' ');
disp(' ');


%% Extraction Data Acquisition

disp('EXTRACTION DATA ACQUISITION');
disp('---------------------------');
disp(' ');

% Particle Data
fileID = fopen([caseFolderExtract, '/extractionData/extractionData']);
extractionData = cell2mat(textscan(fileID, '%f %f %f %f %f %f %f %f %f', 'headerLines', 9, 'delimiter', '\n', 'collectOutput', 1));
fclose(fileID);

if isempty(extractionData)
	error('No Valid Particles Available for Third-Party Initialisation');
else
	disp(['Identified ', num2str(size(extractionData, 1)), ' Particles Available for Third-Party Initialisation']);
end

disp(' ');
disp(' ');

% Flow Sample Data
timeDirs = dir([caseFolderExtract, '/postProcessing/planeSample']);

i = 1;
while i <= size(timeDirs,1)

	if isnan(str2double(timeDirs(i,1).name))
		timeDirs(i,:) = [];
	else
		i = i + 1;
	end

end

if isempty(timeDirs)
	error('No Valid Flow Samples');
else
	disp(['Identified ', num2str(size(timeDirs,1)), ' Flow Variations']);
end

disp(' ');
disp(' ');


%% Initialisation Options

disp('THIRD-PARTY INITIALISATION');
disp('--------------------------');

% Particles Injected Relative to Start of Injection Time
timeOffsetInjection = inputTime('Injection');

% Third-Party Inlet Flow Variations Start at T = 0 [s]
timeOffsetFlow = inputTime('Flow Sampling');

% Particles Injected Upstream of Extraction Point
planeExtraction = planeInput('Extraction');
planeInjection = -planeExtraction;

disp(['    Default (Symmetric) Injection Plane [m]: ', num2str(planeInjection)]);
disp(' ');

valid = false;
while ~valid
    disp(' ');
    selection = input('Modify Default Position? [y/n]: ', 's');
    
    if selection == 'n' | selection == 'N' %#ok<OR2>
        valid = true;
        disp(' ');
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        disp(' ');
        planeInjection = planeInput('Injection');
        valid = true;
        disp(' ');
    else
        disp('    WARNING: Invalid Entry');
    end
    
end

planeOffset = planeExtraction - planeInjection;

% Remove Invalid Particles
extractionData(extractionData(:,1) < min(xDims),:) = [];
extractionData(extractionData(:,1) > max(xDims),:) = [];
extractionData(extractionData(:,2) < min(yDims),:) = [];
extractionData(extractionData(:,2) > max(yDims),:) = [];
extractionData(extractionData(:,3) < min(zDims),:) = [];
extractionData(extractionData(:,3) > max(zDims),:) = [];

% Apply Offsets
extractionData(:,9) = extractionData(:,9) - timeOffsetInjection;
extractionData(:,1) = (extractionData(:,1) - planeOffset) + wallOffset;

disp(' ');


%% Intialisation

tic;
disp('***********');
disp('  Running  ');
disp(' ');
disp(['    Writing Particle Data to ', caseFolderInject, '/constant/parcelInjectionProperties']);

% Remove Existing Particle Data
if exist([caseFolderInject, '/constant/parcelInjectionProperties'], 'file')
	delete([caseFolderInject, '/constant/parcelInjectionProperties']);
end
	
fopen([caseFolderInject, '/constant/parcelInjectionProperties'], 'w');

% OpenFOAM Header
fprintf(fileID, '%s\n', '/*--------------------------------*- C++ -*----------------------------------*\');
fprintf(fileID, '%s\n', '  =========                 |');
fprintf(fileID, '%s\n', '  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox');
fprintf(fileID, '%s\n', '   \\    /   O peration     | Website:  https://openfoam.org');
fprintf(fileID, '%s\n', '    \\  /    A nd           | Version:  7');
fprintf(fileID, '%s\n', '     \\/     M anipulation  |');
fprintf(fileID, '%s\n', '\*---------------------------------------------------------------------------*/');
fprintf(fileID, '%s\n', 'FoamFile');
fprintf(fileID, '%s\n', '{');
fprintf(fileID, '%s\n', '    version     2.0;');
fprintf(fileID, '%s\n', '    format      ascii;');
fprintf(fileID, '%s\n', '    class       dictionary;');
fprintf(fileID, '%s\n', '    location    "constant";');
fprintf(fileID, '%s\n', '    object      scalarListList;');
fprintf(fileID, '%s\n', '}');
fprintf(fileID, '%s\n', '// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //');
fprintf(fileID, '%s\n', ' ');
fprintf(fileID, '%s\n', '	// (x y z) (Ux Uy Uz) d rho mDot particleCount injectionTime');
fprintf(fileID, '%s\n', '(');

% Write Particle Data
rho = 1000;
mDot = 1;
formatSpec = '	(%g %g %g) (%g %g %g) %g %g %g %g %g\n';

for i = 1:size(extractionData,1)
	fprintf(fileID, formatSpec, extractionData(i,1), extractionData(i,2), extractionData(i,3), extractionData(i,4), extractionData(i,5), extractionData(i,6), extractionData(i,7), rho, mDot, extractionData(i,8), extractionData(i,9));
end

% OpenFOAM Footer
fprintf(fileID, '%s\n', ');');
fprintf(fileID, '%s\n', ' ');
fprintf(fileID, '%s\n', '// ************************************************************************* //');

disp(' ');
disp(['    Copying Flow Samples to', caseFolderInject, '/constant/boundaryData/inlet']);

% Remove Existing Sample Data
if ~exist([caseFolderInject, '/constant/boundaryData/inlet'], 'dir')
    mkdir([caseFolderInject, '/constant/boundaryData/inlet']);
else
    rmdir([caseFolderInject, '/constant/boundaryData/inlet'], 's');
    mkdir([caseFolderInject, '/constant/boundaryData/inlet']);
end

% Copy Sample Data
copyfile([caseFolderExtract, '/postProcessing/planeSample/', timeDirs(1,1).name, '/thirdPartyExtraction/faceCentres'], [caseFolderInject, '/constant/boundaryData/inlet/points']);

for i = 1:size(timeDirs,1)
	inletTime = num2str(str2double(timeDirs(i,1).name) - timeOffsetFlow);
	
	mkdir([caseFolderInject, '/constant/boundaryData/inlet/', inletTime]);
	copyfile([caseFolderExtract, '/postProcessing/planeSample/', timeDirs(i,1).name, '/thirdPartyExtraction/scalarField/nut'], [caseFolderInject, '/constant/boundaryData/inlet/', num2str(inletTime), '/nut']);
	copyfile([caseFolderExtract, '/postProcessing/planeSample/', timeDirs(i,1).name, '/thirdPartyExtraction/scalarField/nuTilda'], [caseFolderInject, '/constant/boundaryData/inlet/', num2str(inletTime), '/nuTilda']);
	copyfile([caseFolderExtract, '/postProcessing/planeSample/', timeDirs(i,1).name, '/thirdPartyExtraction/vectorField/U'], [caseFolderInject, '/constant/boundaryData/inlet/', num2str(inletTime), '/U']);
end

executionTime = toc;

disp(' ');
disp(['    Initialisation Time: ', num2str(executionTime), 's']);
disp(' ');
disp('  Success  ');
disp('***********');


%% Cleaning

clearvars;
disp(' ');


%% Local Functions

function T = inputTime(type)

	valid = false;
	while ~valid
		disp(' ');
		T = str2double(input(['Enter Extraction Case ', type, ' Start Time [s]: '], 's'));

		if isnan(T)
			disp('    WARNING: Invalid Time Entry');
		elseif length(T) > 1
			disp('    WARNING: Invalid Time Entry');
		else
			valid = true;
			disp(' ');
		end

	end
	
end

function P = planeInput(type)

    valid = false;
    while ~valid
        disp(' ');
        P = str2double(input(['Enter x-Coordinate of ', type, ' Plane [m]: '], 's'));

        if isnan(P) || length(P) > 1
            disp('        WARNING: Invalid Plane Position');
        else
            valid = true;
        end

    end
    
end
