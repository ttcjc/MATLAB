%% Particle Tracker v1.0

clear variables;
close all;
clc;

disp ('=====================');
disp ('Particle Tracker v1.0');
disp ('=====================');
disp (' ');


%% Changelog

% v1.0 - Initial Commit


%% Case Selection

disp('CASE SELECTION');
disp('--------------');
disp(' ');

disp('Select Extract Case:');
% caseFolder = uigetdir('~/OpenFOAM');
caseFolder = ('/home/cjc96/Mount/Uni/OpenFOAM/ttcjc-7/run/Test_Block');
disp(['Case: ', caseFolder]);
disp(' ');
[timeDirs, deltaT] = timeDirectories(caseFolder);
disp(' ');
disp(' ');

disp('GEOMETRY SELECTION');
disp('------------------');
disp(' ');

disp('Select Case Geometry:');
% [file, path] = uigetfile(caseFolder);
file = 'Test_Block.stl';
path = '/home/cjc96/CAD/CFD Geometries/';
disp(['Geometry: ', path, file]);
geometry = stlread([path, file]);
disp(' ');
disp(' ');


%% Data Acquisition

disp('LAGRANGIAN DATA COLLECTION');
disp('--------------------------');

valid = false;
while ~valid
    disp(' ');
    selection = input('Load Saved Particle Data? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        disp(' ');
        [particleData, nParticles] = lagrangianData(caseFolder, timeDirs);        
        valid = true;
        disp(' ');
        disp(' ');
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        particleData = uigetfile('~/MATLAB/CFD/particleData');
        disp(['    Loading: ', particleData]);
        load(particleData)
        valid = true;
        disp(' ');
        disp(' ');
    else
        disp('    WARNING: Invalid Entry');
    end

end

%% Tracking Method

disp('TRACKING METHOD');
disp('----------------');

valid = false;
while ~valid
    disp(' ');
    disp('Possible Tracking Methods:');
    disp('    A: Self-Contamination Case');
    disp('    B: Third-Party Contamination Case');
    disp(' ');
    selection = input('Select Desired Tracking Method [A/B]: ', 's');
    
    if selection == 'a' | selection == 'A' %#ok<OR2>
        method = 'A';
        valid = true;
    elseif selection == 'b' | selection == 'B' %#ok<OR2>
        method = 'B';
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end
    
end

valid = false;
while ~valid
    disp(' ');
    selection = input('Enable Position Interpolation? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        interp = false;
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        interp = true;          
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end

valid = false;
while ~valid
    disp(' ');
    selection = input('Filter Particle Diameters? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        minD = min(cell2mat(particleData{2,end}(2:end,5))) * 1e6;
        maxD = max(cell2mat(particleData{2,end}(2:end,5))) * 1e6;
        valid = true;
        disp(' ');
        disp(' ');
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        minD = inputD('Minimum');
        maxD = inputD('Maximum');
        valid = true;
        disp(' ');
        disp(' ');
    else
        disp('    WARNING: Invalid Entry');
    end

end


%% Particle Tracking

tic;

particleActive = cell2mat(particleData{2,end}(2:end,1));
particleD = cell2mat(particleData{2,end}(2:end,5)) * 1e6;

trackingTemp = horzcat(particleActive, particleD);
trackingTemp = find(~trackingTemp(:,1) & trackingTemp(:,2) >= minD & trackingTemp(:,2) <= maxD)+1;

contaminantActive = cell2mat(particleData{2,end}(trackingTemp,1));
contaminantID = cell2mat(particleData{2,end}(trackingTemp,2));
contaminantProc = cell2mat(particleData{2,end}(trackingTemp,3));
contaminantPos = cell2mat(particleData{2,end}(trackingTemp,4));
contaminantD = cell2mat(particleData{2,end}(trackingTemp,5)) * 1e6;

trackingData = cell(size(contaminantID,1),size(particleData,2));

for i = 1:size(particleData,2)
    trackingTemp = intersect(cell2mat(particleData{2,i}(2:end,2:3)), horzcat(contaminantID, contaminantProc), 'rows');
    trackingTemp = find(ismember(cell2mat(particleData{2,i}(2:end,2:3)), trackingTemp, 'rows'))+1;

    for j = 1:size(trackingTemp,1)
        trackingData(j,i) = particleData{2,i}(trackingTemp(j,1),4);
    end

end

switch method
    
    case 'A'
        
        disp('TRACKING SELF-CONTAMINATION');
        disp('---------------------------');
        disp(' ');
        disp([num2str(size(contaminantID,1)), ' Particles Deposited on Geometry']);
        disp('    PLOTTING...');

        figure(1);
        hold on;
        set(figure(1), 'outerPosition', [50, 50, 750, 750]);
        set(gca, 'units', 'normalized', 'position', [0.125, 0.125, 0.75, 0.75]);

        title({'Self-Contamination' ; 'Impinging Particle Tracks' ; ' '}, 'interpreter', 'none');
        xlabel('x [m]');
        ylabel('y [m]');
        zlabel('z [m]');

        view([45 15]);
        
        axis('image');
        colormap viridis;
        caxis([minD maxD])
        c = colorbar;
        c.Label.String = 'Particle Diameter [um]';
        faceColor = [0.5 0.5 0.5];
        edgeColor = [0.5 0.5 0.5];

        trisurf(geometry, 'faceColor', faceColor, 'edgeColor', edgeColor);
    
    case 'B'
        
        disp('TRACKING THIRD-PARTY CONTAMINATION');
        disp('---------------------------');
        disp(' ');
        disp([num2str(size(contaminantID,1)), ' Particles Deposited on Geometry']);
        disp('    PLOTTING...');

        figure(1);
        hold on;
        set(figure(1), 'outerPosition', [50, 50, 750, 750]);
        set(gca, 'units', 'normalized', 'position', [0.125, 0.125, 0.75, 0.75]);

        title({'Third-Party Contamination' ; 'Impinging Particle Tracks' ; ' '}, 'interpreter', 'none');
        xlabel('x [m]');
        ylabel('y [m]');
        zlabel('z [m]');

        view([-45 15]);

        axis('image');
        colormap viridis;
        c = colorbar;
        caxis([minD maxD])
        c.Label.String = 'Particle Diameter [um]';
        faceColor = [0.5 0.5 0.5];
        edgeColor = [0.5 0.5 0.5];

        trisurf(geometry, 'faceColor', faceColor, 'edgeColor', edgeColor);
        
end

scatter3(contaminantPos(:,1), contaminantPos(:,2), contaminantPos(:,3), 100, contaminantD(:,1));

for i = 1:size(trackingData,1)
    path = unique(cell2mat(trackingData(i,:)'), 'stable', 'rows');

    if interp && size(path,1) > 1
        path = interparc(round(deltaT/0.001), path(:,1), path(:,2), path(:,3), 'spline');
    end

    posX = path(:,1);
    posY = path(:,2);
    posZ = path(:,3);

    map = colormap;
    index = round(1+(size(map,1)-1)*(contaminantD(i,1)-minD)/(maxD-minD));

    plot3(posX, posY, posZ, 'color', map(index,:));
end

hold off;

executionTime = toc;
disp(' ');
disp(['Tracking Time = ', num2str(executionTime), 's']);

disp(' ');
disp('*****************');
disp('TRACKING COMPLETE');
disp('*****************');


%% Local Functions

function D = inputD(type)

    valid = false;
    while ~valid
        D = str2double(input(['    ', type, ' Diameter of Interest [', char(956), 'm]: '], 's'));

        if isnan(D)
            disp('        WARNING: Invalid Entry');
            disp(' ');
        elseif length(D) > 1
            disp('        WARNING: Invalid Entry');
            disp(' ');
        else
            valid = true;
        end

    end

end
