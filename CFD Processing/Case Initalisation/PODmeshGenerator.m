%% POD Mesh Generator v1.0

clearvars;
close all;
clc;

fig = 0; %#ok<*NASGU>
figHold = 0; %#ok<*NASGU>

disp ('=======================');
disp ('POD MESH GENERATOR v1.0');
disp ('=======================');
disp (' ');


%% Changelog

% v1.0 - Initial Commit


%% Initialisation

% Mesh Settings
cellSize = 8e-3;			% m
xRange = [0.31875; 1.504];	% m
yRange = [-0.2585; 0.2585];	% m
zRange = [0.001; 0.403];	% m

% Delete Previously Generated Mesh
if exist(['Output/Files/PODmesh_', num2str(cellSize * 1e3), 'mm'], 'file')
	delete(['Output/Files/PODmesh_', num2str(cellSize * 1e3), 'mm']);
end

% Generate New Mesh
tic;
disp('***********');
disp('  Running  ');
disp(' ');
disp('    Generating POD Mesh');

[x, y, z] = meshgrid(min(xRange):cellSize:max(xRange), min(yRange):cellSize:max(yRange), min(zRange):cellSize:max(zRange));
meshPoints = [x(:), y(:), z(:)];

% Remove Mesh Points Intersecting Vehicle
disp(' ');
disp('    Removing Points Intersecting Vehicle');

[file, path] = uigetfile('~/CAD/CFD Geometries/*.*', 'Select Subject Geometry', 'multiSelect', 'on');

if isa(file, 'cell')

	for i = 1:size(file,2)
		part = file{1,i}(1:end-4);
		geometry.(part) = stlreader([path,file{1,i}]);
		
		X = unique(geometry.(part).vertices, 'rows'); % Unique points in CAD model
		TRI = delaunayTriangulation(X);
		XI = meshPoints;

		intersect = find(~isnan(tsearchn(X, TRI.ConnectivityList, XI)));
		meshPoints(intersect,:) = [];
	end

else
	part = file(1:end-4);
	geometry.(part) = stlreader([path, file]);
	
	X = unique(geometry.(part).vertices, 'rows'); % Unique points in CAD model
	TRI = delaunayTriangulation(X);
	XI = meshPoints;

	intersect = find(~isnan(tsearchn(X, TRI, XI)));
	meshPoints(intersect,:) = [];
end

% Write Mesh
disp(' ');
disp(['    Writing POD Mesh to ~/MATLAB/Output/Files/PODmesh_', num2str(cellSize * 1e3), 'mm']);

fileID = fopen(['Output/Files/PODmesh_', num2str(cellSize * 1e3), 'mm'], 'w');
formatSpec = '                (%g %g %g)\n';

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

% Figure Setup
fig = fig + 1;
figure('name', 'POD Volume Probes');
hold on;
set(figure(fig), 'outerPosition', [25, 25, 800, 800]);

% Plot
parts = fieldnames(geometry);

for i = 1:size(parts,1)
	part = parts{i,1};
	patch(geometry.(part), 'faceColor', [0.5, 0.5, 0.5], 'edgeColor', [0.5, 0.5, 0.5]);
end

scatter3(meshPoints(:,1), meshPoints(:,2), meshPoints(:,3), 5, [0.21176, 0.06667, 0.38824], 'filled')

% Figure Formatting
view([45, 15]);
axis off
set(gca, 'units', 'normalized', 'position', [0.1275, 0.1275, 0.745, 0.745], ...
		 'fontName', 'LM Roman 12', 'fontSize', 12, 'layer', 'top', ...
		 'dataAspectRatio', [1, 1, 1]);
hold off;


%% Cleaning

clearvars -except cellSize meshPoints;
disp(' ');