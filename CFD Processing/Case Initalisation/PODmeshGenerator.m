%% POD Mesh Generator v1.0

clearvars;
close all;
clc;

disp ('=======================');
disp ('POD MESH GENERATOR v1.0');
disp ('=======================');
disp (' ');


%% Changelog

% v1.0 - Initial Commit


%% Initialisation

% Settings
cellSize = 8e-3;
xRange = [0.31875; 1.504];
yRange = [-0.2585; 0.2585];
zRange = [0; 0.403];

% Generate Mesh
[x, y, z] = meshgrid(min(xRange):cellSize:max(xRange), min(yRange):cellSize:max(yRange), min(zRange):cellSize:max(zRange));
meshPoints = [x(:), y(:), z(:)];

% Delete Previously Generated Mesh
if exist('Output/Files/PODmesh', 'file')
	delete('Output/Files/PODmesh')
end

% Write Mesh
tic;
disp('***********');
disp('  Running  ');
disp(' ');
disp(['    Writing POD Mesh to ~/MATLAB/Output/Files/PODmesh']);

fileID = fopen('Output/Files/PODmesh', 'w');
formatSpec = '\t\t(%g %g %g)\n';

for i = 1:size(meshPoints,1)
	fprintf(fileID, formatSpec, meshPoints(i,1), meshPoints(i,2), meshPoints(i,3));
end

fclose(fileID);

executionTime = toc;

disp(' ');
disp(['    Write Time: ', num2str(executionTime), 's']);
disp(' ');
disp('  Success  ');
disp('***********');


%% Cleaning

clearvars -except cellSize points
disp(' ');