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

caseFolderExtract = uigetdir('~/OpenFOAM', 'Extraction Case');
disp(['Extraction Case: ', caseFolderExtract]);
disp(' ');

disp('Select Injection Case:');
caseFolderInject = uigetdir('~/OpenFOAM', 'Injection Case');
disp(['Injection Case: ', caseFolderInject]);
disp(' ');
disp(' ');

% Confirm Valid Case Directories
disp('Analysing Case Structure:');

if exist([caseFolderExtract, '/0'], dir) && exist([caseFolderExtract, '/extractionData/extractionData'], file)
	fileID = fopen([caseFolderExtract, '/extractionData/extractionData']);
	content = textscan(fileID, '%s', 'headerLines', 15, 'delimiter', '\n', 'collectOutput', 1);
	
	if
		moo
	else
		moo
	end
	
else
	error(['Invalid Case Directory (', caseFolderExtract, ')']);
end


%% Data Acquisition

disp('EXTRACTION DATA COLLECTION');
disp('--------------------------');

% Particle Data
if exist([caseFolderExtract, '/extractionData/kinematicCloud_extractionData'], 'file')
    fileID = fopen([caseFolderExtract, '/extractionData/kinematicCloud_extractionData']);
    dataFile = textscan(fileID, '%s %f %s %f %s %f %s %f %s %f %s %f %s %f %s %f', 'delimiter', ' ', 'multipleDelimsAsOne', 1);

    disp(['Identified ', num2str(size(dataFile{1,1},1)), ' Particles for Third-Party Initialisation']);
    disp(' ');
else
    error('No Extracted Particle Data in Specified Directory');
end

particleVars = {'x', 'y', 'z', 'Ux', 'Uy', 'Uz', 'd', 't'};
extractionData = cell(size(dataFile{1,1},1)+1, size(particleVars,2));

for i = 1:size(particleVars,2)
    extractionData{1,i} = particleVars(i);
end

i = 1;
for j = 2:2:size(dataFile,2)

    for k = 1:size(dataFile{1,1},1)
        extractionData{k+1,i} = dataFile{1,j}(k);
    end

    i = i + 1;
end

% Flow Sample Data
if exist([caseFolderExtract, '/postProcessing/planes'], 'dir')
    timeDirs = dir([caseFolderExtract, '/postProcessing/planes']);
else
    error('No Flow Sample Data in Specified Directory');
end

i = 1;
while i <= size(timeDirs,1)
    
    if isnan(str2double(timeDirs(i,1).name))
        timeDirs(i,:) = [];
    else
        i = i + 1;
    end
    
end

disp(['Identified ', num2str(size(timeDirs,1)), ' Inlet Flow Variations']);
disp(' ');
disp(' ');


%% Initialisation

disp('INITIALISATION');
disp('--------------');

% Particle Initialisation
valid = false;
while valid == false
    disp(' ');
    timeOffset = str2double(input('Enter Injection Start Time [s]: ', 's'));

    if isnan(timeOffset)
        disp('    WARNING: Invalid Time Entry');
    elseif length(timeOffset) > 1
        disp('    WARNING: Invalid Time Entry');
    else
        valid = true;
        disp(' ');
    end

end

planeExtraction = planeInput('Extraction');

planeInjection = -1.044 - planeExtraction;
disp(' ');
disp(['Default (Symmetric) Injection Position: [m]', num2str(planeInjection)]);

valid = false;
while valid == false
    disp(' ');
    modify = input('    Modify Default Position? [y/n]: ', 's');
    
    if modify == 'n' | modify == 'N' %#ok<OR2>
        valid = true;
        disp(' ');
    elseif modify == 'y' | modify == 'Y' %#ok<OR2>
        disp(' ');
        planeInjection = planeInput('Injection');
        valid = true;
        disp(' ');
    else
        disp('        WARNING: Invalid Entry');
    end
    
end

planeOffset = planeExtraction - planeInjection;

for i = 2:size(extractionData,1)
    extractionData{i,1} = extractionData{i,1} - planeOffset;
    extractionData{i,8} = extractionData{i,8} - timeOffset;
end

disp('Translating Particles...');

fileID = fopen([caseFolderInject, '/constant/parcelInjectionProperties'], 'w');

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
fprintf(fileID, '%s\n', '(');
fprintf(fileID, '%s\n', '	// (x y z) (Ux Uy Uz) d rho mDot* injectionTime** particleCount');
fprintf(fileID, '%s\n', '	// * Unused');
fprintf(fileID, '%s\n', '	// ** Relative to start of injection');

formatSpec = '	(%g %g %g) (%g %g %g) %g %g %g %g %g\n';
for i = 2:size(extractionData,1)
    fprintf(fileID, formatSpec, extractionData{i,1}, extractionData{i,2}, extractionData{i,3}, extractionData{i,4}, extractionData{i,5}, extractionData{i,6}, extractionData{i,7}, 1000, 1, extractionData{i,8}, 1);
end

fprintf(fileID, '%s\n', ');');
fprintf(fileID, '%s\n', ' ');
fprintf(fileID, '%s\n', '// ************************************************************************* //');

disp('*************************************************************************');
disp('Extracted Particles Successfully Translated and Written to Injection Case');
disp('*************************************************************************');

% Inlet Flow Initialisation
if exist([caseFolderInject, '/constant/boundaryData/inlet'], 'dir') == 0
    mkdir([caseFolderInject, '/constant/boundaryData/inlet']);
else
    rmdir([caseFolderInject, '/constant/boundaryData/inlet'], 's');
    mkdir([caseFolderInject, '/constant/boundaryData/inlet']);
end

disp(' ');
disp('Copying Inlet Flow Variations...');

copyfile([caseFolderExtract, '/postProcessing/planes/', timeDirs(1,1).name, '/thirdPartyExtraction/faceCentres'], [caseFolderInject, '/constant/boundaryData/inlet/points']);

for i = 1:size(timeDirs,1)
    mkdir([caseFolderInject, '/constant/boundaryData/inlet/', timeDirs(i,1).name]);
    copyfile([caseFolderExtract, '/postProcessing/planes/', timeDirs(i,1).name, '/thirdPartyExtraction/scalarField/nut'], [caseFolderInject, '/constant/boundaryData/inlet/', timeDirs(i,1).name, '/nut']);
    copyfile([caseFolderExtract, '/postProcessing/planes/', timeDirs(i,1).name, '/thirdPartyExtraction/scalarField/nuTilda'], [caseFolderInject, '/constant/boundaryData/inlet/', timeDirs(i,1).name, '/nuTilda']);
    copyfile([caseFolderExtract, '/postProcessing/planes/', timeDirs(i,1).name, '/thirdPartyExtraction/vectorField/U'], [caseFolderInject, '/constant/boundaryData/inlet/', timeDirs(i,1).name, '/U']);
end

disp('************************************************************');
disp('Inlet Flow Variations Successfully Written to Injection Case');
disp('************************************************************');


%% Clean-up

% clear


%% Local Functions

function plane = planeInput(type)

    valid = false;
    while valid == false
        disp(' ');
        plane = str2double(input(['Enter X-Position of ', type, ' Plane [m]: '], 's'));

        if isnan(plane)
            disp('    WARNING: Invalid Plane Position');
        elseif length(plane) > 1
            disp('    WARNING: Invalid Plane Position');
        else
            valid = true;
        end

    end
    
end
