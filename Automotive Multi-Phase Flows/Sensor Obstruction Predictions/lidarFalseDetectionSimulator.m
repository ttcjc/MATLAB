%% Lidar False Detection Simulator
% ----
% Lorem ipsum


%% Preamble

run preamble;

samplesPerCell = 32; % Used During Numerical Integration

nSteps = 1e5; % Used To Determine 'dT' for the Lidar Range Integral

figSave = false; % Save .fig File(s)

disp('====================================');
disp('Lidar False Detection Simulator v1.0');
disp('====================================');

disp(' ');
disp(' ');


%% Changelog

% v1.0 - Initial Commit


%% Select Region of Interest

disp('Region of Interest');
disp('-------------------');

disp(' ');

disp('Possible Regions of Interest:');
disp('    A: Near Wake');
disp('    B: Mid Wake');
disp('    C: Far Wake (Full-Scale Only)');

valid = false;
while ~valid
    disp(' ');
    
    selection = input('Select Region of Interest [A/B]: ', 's');

    if selection == 'a' | selection == 'A' %#ok<OR2>
        format = 'A';
        
        valid = true;
    elseif selection == 'b' | selection == 'B' %#ok<OR2>
        format = 'B';
        
        valid = true;
    elseif selection == 'c' | selection == 'C' %#ok<OR2>
        format = 'C';
        
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;

disp(' ');
disp(' ');


%% Acquire Volume Field

disp('Volume Field Acquisition');
disp('-------------------------');

valid = false;
while ~valid
    disp(' ');
    [fileName, filePath] = uigetfile([saveLoc, '/Numerical/MATLAB/volumeField/*.mat'], ...
                                      'Select Volumetric Data');
    
    switch format
        
        case 'A'
            
            if contains(filePath, '/nearWake')
                valid = true;
            else
                disp('WARNING: Invalid File Selection');
                clear fileName filePath;
            end
            
        case 'B'
            
            if contains(filePath, '/midWake')
                valid = true;
            else
                disp('WARNING: Invalid File Selection');
                clear fileName filePath;
            end
            
        case 'C'
            
            if contains(filePath, '/farWake')
                valid = true;
            else
                disp('WARNING: Invalid File Selection');
                clear fileName filePath;
            end
            
    end
    
    disp(['Loading ''', fileName, '''...']);
    
    campaignID = load([filePath, fileName], 'campaignID').campaignID;
    caseID = load([filePath, fileName], 'caseID').caseID;
    dataID = load([filePath, fileName], 'dataID').dataID;
    volumeData = load([filePath, fileName], 'volumeData').volumeData;
    cellSize = load([filePath, fileName], 'cellSize').cellSize;
    sampleInt = load([filePath, fileName], 'sampleInt').sampleInt;
    timePrecision = load([filePath, fileName], 'timePrecision').timePrecision;
    dLims = load([filePath, fileName], 'dLims').dLims;
    
    disp('    Success');
end
clear valid;

disp(' ');
disp(' ');


%% Select Relevant Geometry and Define Bounding Box

[geometry, xDims, yDims, zDims, spacePrecision, normLength] = selectGeometry(geoLoc);

disp(' ');
disp(' ');


%% Select Plane of Interest

% Select Plane of Interest
lidarData.PoV.targetPlane = identifyVolumeSlices(volumeData.positionGrid, spacePrecision, false);
planeID = fieldnames(lidarData.PoV.targetPlane); planeID = planeID{1};
lidarData.PoV.targetPlane = lidarData.PoV.targetPlane.(planeID);
clear planeID;

% Extract Planar Position Data
orientation = lidarData.PoV.targetPlane.orientation;

switch orientation
    
    case 'YZ'
        index = find(volumeData.positionGrid(:,1) == lidarData.PoV.targetPlane.position);
        
    case 'XZ'
        index = find(volumeData.positionGrid(:,2) == lidarData.PoV.targetPlane.position);
        
    case 'XY'
        index = find(volumeData.positionGrid(:,3) == lidarData.PoV.targetPlane.position);
        
end

lidarData.positionGrid = volumeData.positionGrid(index,:);

nCells = height(lidarData.positionGrid);

disp(' ');
disp(' ');


%% Select Origin Point

disp('Ray Origin Definition');
disp('----------------------');

% Select Ray Origin Point
lidarData.PoV.originPoint = zeros([1,3]);

valid = false;
while ~valid
    disp(' ')
    disp('Specify Ray Origin Point:')
    
    lidarData.PoV.originPoint(1) = inputPos('X');
    lidarData.PoV.originPoint(2) = inputPos('Y');
    lidarData.PoV.originPoint(3) = inputPos('Z');

    if (lidarData.PoV.originPoint(1) < min(volumeData.positionGrid(:,1)) || ...
        lidarData.PoV.originPoint(1) > max(volumeData.positionGrid(:,1))) || ...
       (lidarData.PoV.originPoint(2) < min(volumeData.positionGrid(:,2)) || ...
        lidarData.PoV.originPoint(2) > max(volumeData.positionGrid(:,1))) || ...
       (lidarData.PoV.originPoint(3) < min(volumeData.positionGrid(:,3)) || ...
        lidarData.PoV.originPoint(3) > max(volumeData.positionGrid(:,1)))
        disp('        WARNING: Origin Point Lies Outside Volume');
        
        continue;
    end
    
    switch orientation
        
        case 'YZ'
            
            if lidarData.PoV.originPoint(1) < lidarData.PoV.targetPlane.position
                disp('        WARNING: Origin Point Must Lie in the Positive X-Direction of the Target Plane');
                
                continue;
            end
            
        case 'XZ'
            
            if lidarData.PoV.originPoint(2) > lidarData.PoV.targetPlane.position
                disp('        WARNING: Origin Point Must Lie in the Negative Y-Direction of the Target Plane');
                
                continue;
            end
            
        case 'XY'
            
            if lidarData.PoV.originPoint(3) > lidarData.PoV.targetPlane.position
                disp('        WARNING: Origin Point Must Lie in the Positive Z-Direction of the Target Plane');
                
                continue;
            end
            
    end
    
    valid = true;    
end
clear valid;

disp(' ');
disp(' ');


%% Select Ray(s) of Interest

% RoI = dsearchn(lidarData.positionGrid(:,[2,3]), [0, 0.76]);
RoI = dsearchn(lidarData.positionGrid(:,[2,3]), [1.945, 0.76]);


%% Calculate Mass Along Ray(s)

disp('Mass Along Ray Calculation');
disp('---------------------------');

disp(' ');

disp('***********');
disp('  RUNNING ');

tic;

maxNumCompThreads(nProc);

evalc('parpool(''threads'');');

%%%%

disp(' ');

disp('    Initialising...');

nTimes = height(volumeData.time);

dL = cellSize.target / samplesPerCell;

% Check if Plane of Interest Intersects Geometry
removeIntersect = false;

switch orientation
    
    case 'YZ'
        
        if lidarData.PoV.targetPlane.position <= xDims(2)
            removeIntersect = true;
        end
        
    case 'XZ'
        
        if lidarData.PoV.targetPlane.position >= yDims(1)
            removeIntersect = true;
        end
        
    case 'XY'
        
        if lidarData.PoV.targetPlane.position <= zDims(2)
            removeIntersect = true;
        end
        
end

% Remove Erroneous Data From Cells Intersecting Geometry
if removeIntersect
    disp('        Removing Erroneous Data From Grid Cells Intersecting Geometry...');

    % Perform Removal
    volumeDataVars = fieldnames(volumeData);
    nonFieldVars = {'positionGrid'; 'time'};
    fieldVars = setdiff(volumeDataVars, nonFieldVars);
    clear volumeDataVars nonFieldVars;
    
    parts = fieldnames(geometry);
    for i = 1:height(parts)
        DT = delaunay(geometry.(parts{i}).vertices);

        index = ~isnan(tsearchn(geoPoints, DT, volumeData.positionGrid));

        for j = 1:height(fields)
            volumeData.(fieldVars{j}).mean(index,:) = NaN;

            for k = 1:nTimes
                volumeData.(fieldVars{j}).inst{k}(index,:) = NaN;
            end

        end
        clear j;

    end
    clear i parts;
    
    clear fieldVars;    
end

% Update Map Boundaries
switch orientation

    case 'YZ'
        xLimsData = lidarData.PoV.targetPlane.position;
        yLimsData = [min(lidarData.positionGrid(:,2)); max(lidarData.positionGrid(:,2))];
        zLimsData = [min(lidarData.positionGrid(:,3)); max(lidarData.positionGrid(:,3))];
        
    case 'XZ'
        xLimsData = [min(lidarData.positionGrid(:,1)); max(lidarData.positionGrid(:,1))];
        yLimsData = lidarData.PoV.targetPlane.position;
        zLimsData = [min(lidarData.positionGrid(:,3)); max(lidarData.positionGrid(:,3))];
        
    case 'XY'
        xLimsData = [min(lidarData.positionGrid(:,1)); max(lidarData.positionGrid(:,1))];
        yLimsData = [min(lidarData.positionGrid(:,2)); max(lidarData.positionGrid(:,2))];
        zLimsData = lidarData.PoV.targetPlane.position;
        
end

disp (' ');

% Reshape Position Data for Improved Interpolation Performance
disp('    Reshaping Position Data for Improved Interpolation Performance...');

gridShape = [height(unique(volumeData.positionGrid(:,1))), ...
             height(unique(volumeData.positionGrid(:,2))), ...
             height(unique(volumeData.positionGrid(:,3)))];

x = reshape(volumeData.positionGrid(:,1), gridShape);
y = reshape(volumeData.positionGrid(:,2), gridShape);
z = reshape(volumeData.positionGrid(:,3), gridShape);

disp(' ');

% Calculate Instantaneous Line of Sight
disp('    Calculating Instantaneous Obstruction...');

lidarData.time = volumeData.time;

% Initialise Progress Bar
wB = waitbar(0, 'Calculating Instantaneous Obstruction', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';
dQ = parallel.pool.DataQueue;
afterEach(dQ, @parforWaitBar);
parforWaitBar(wB, nTimes);

% Perform Calculation
sampleDist = cell(nTimes,1); sampleDist(:) = {cell(height(RoI),1)};
nParticlesDensity = sampleDist;
d20 = sampleDist;
density = sampleDist;

nParticlesVolume = volumeData.nParticles.inst;
d20Volume = volumeData.d20.inst;
densityVolume = volumeData.density.inst;
cellVolume = cellSize.volume;
positionGrid = lidarData.positionGrid(RoI,:);
pointPosition = lidarData.PoV.originPoint;
parfor i = 1:nTimes
    nParticlesField = reshape(full(nParticlesVolume{i}), gridShape) / cellVolume;
    d20Field = reshape(full(d20Volume{i}), gridShape);
    densityField = reshape(full(nParticlesVolume{i}), gridShape);
    
    nParticlesInterp = griddedInterpolant(x, y, z, nParticlesField, 'linear', 'none');
    d20Interp = griddedInterpolant(x, y, z, d20Field, 'linear', 'none');
    densityInterp = griddedInterpolant(x, y, z, densityField, 'linear', 'none');
    
    for j = 1:height(positionGrid)
        dirVec = positionGrid(j,:) - pointPosition;
        distFull = sqrt(dirVec(1)^2 + dirVec(2)^2 + dirVec(3)^2);

        dist = (dL:dL:distFull)';
        samplePoints = pointPosition + (dist * (dirVec / distFull));
        
        nParticlesDensity{i}{j} = nParticlesInterp(samplePoints(:,1), samplePoints(:,2), samplePoints(:,3));
        d20{i}{j} = d20Interp(samplePoints(:,1), samplePoints(:,2), samplePoints(:,3));
        density{i}{j} = densityInterp(samplePoints(:,1), samplePoints(:,2), samplePoints(:,3));
        
        sampleDist{i}{j} = abs(sqrt(samplePoints(:,1).^2 + samplePoints(:,2).^2 + samplePoints(:,3).^2) - ...
                               sqrt(samplePoints(1,1)^2 + samplePoints(1,2)^2 + samplePoints(1,3)^2));
    end
    
    % Remove Unnecessary Data
    nParticlesVolume{i} = [];
    d20Volume{i} = [];
    densityVolume{i} = [];
    
    % Update Waitbar
    send(dQ, []);
end
clear nParticlesVolume d20Volume densityVolume cellVolume positionGrid pointPosition;

delete(wB);

% clear volumeData;

lidarData.nParticlesDensity.inst = nParticlesDensity; clear nParticlesDensity;
lidarData.d20.inst = d20; clear d20;
lidarData.density.inst = density; clear density;

lidarData.sampleDist = sampleDist{1}{1}; clear sampleDist;

%%%%

evalc('delete(gcp(''nocreate''));');

executionTime = toc;

disp(' ');

disp(['    Run Time: ', num2str(executionTime), 's']);

disp(' ');

disp('  SUCCESS  ');
disp('***********');

disp(' ');
disp(' ');


%% Calculate Mie Efficiencies

n_a = 1.00027451; % Refractive Index of Air
n_w = complex(1.32352, 5.15e-7); % Refractive Index of Water

c = 299792458 / n_a; % Speed of Light in Air

k = (tau / 905e-9) * n_a; % Wave Number in Air

D = (1e-7:1e-7:400e-6)'; % Range of Particle Diameters

Q_ext_Range = zeros([height(D),1]);
Q_bck_Range = zeros([height(D),1]);

for i = 1:height(D)
    x = k * (D(i) / 2); % Size Parameter
    
    mieOut = mieScattering(n_w, x);
    
    Q_ext_Range(i) = mieOut(1);
    Q_bck_Range(i) = mieOut(4);
end
clear i x mieOut;

Q_ext_Interp = griddedInterpolant(D, Q_ext_Range, 'nearest', 'nearest');
Q_bck_Interp = griddedInterpolant(D, Q_bck_Range, 'nearest', 'nearest');


%% Solve Lidar Range Equation

alpha0 = 1e-4; % Extinction Coefficient for a Clear Sky (1 / m)
beta0 = 0; % Backscattering Coefficient for a Clear Sky
Gamma_O = 0.2; % Target Reflectivity
A_O = inf; % Target Area
eta_T = 1; % Transmitter Efficiency
eta_R = 1; % Receiver Efficiency
P0 = 1; % Peak Laser Power
t_T = 5e-9; % Laser Pulse Width
axisDisp = 27e-3; % Displacement Between Optical Channels
D_T0 = 4e-3; % Initial Beam Diameter
D_R0 = 40e-3; % Receiver Diameter
A_R0 = (tau * D_R0^2) / 8; % Receiever Area
gamma_T = deg2rad(0.06); % Beam Divergence
gamma_R = deg2rad(0.1); % Divergence of Receiver Channel

R1 = (axisDisp - (D_T0 / 2) - (D_R0 / 2)) / (tan(gamma_T / 2) + tan(gamma_R / 2));
R2 = (axisDisp - (D_R0 / 2) + (D_T0 / 3)) / (tan(gamma_R / 2) - tan(gamma_T / 2));

for i = 1:height(RoI)    
    R_discrete = (0:dL:20)';
%     index_O = find(R_discrete >= (4 * 4.176), 1, 'first');
%     R_O = R_discrete(index_O);
    index_0 = [];
    R_O = 100;
    
    dT = (2 * (max(R_discrete) / c)) / nSteps;
    
%     maxP_R = 0;
    
    for j = 1:4 % [25, 50, 75, 100] % 1:nTimes
        P_R = zeros([height(R_discrete),1]);
        
        % Initialise Spray Data Interpolants
        nParticlesDensityInterp = griddedInterpolant(lidarData.sampleDist, ...
                                                     lidarData.nParticlesDensity.inst{j}{i}, ...
                                                     'linear', 'nearest');
        
        nParticlesDensity = @(R) (nParticlesDensityInterp(R) .* ((R >= 0) & (R <= max(lidarData.sampleDist))));
%         nParticlesDensity = @(R) (0);
        
        d20Interp = griddedInterpolant(lidarData.sampleDist, ...
                                     (lidarData.d20.inst{j}{i} * 1e-6), ...
                                     'linear', 'nearest');
        
        d20 = @(R) (d20Interp(R) .* ((R >= 0) & (R <= max(lidarData.sampleDist))));
%         d20 = @(R) (0);

        densityInterp = griddedInterpolant(lidarData.sampleDist, ...
                                           lidarData.density.inst{j}{i}, ...
                                           'linear', 'nearest');
        
        density = @(R) (densityInterp(R) .* ((R >= 0) & (R <= max(lidarData.sampleDist))));
%         density = @(R) (0);
        
        % Calculate Transmitted Signal
        P_T = @(t) ((P0 .* (sin((tau ./ 2) .* (t ./ t_T)).^2)) .* ((t >= 0) & (t <= t_T)));

        % Calculate Diameter of Transmitted Beam
        D_T = @(R) (D_T0 + (2 .* R .* tan(gamma_T ./ 2)));

        % Calculate Diameter of Receiver Channel
        D_R = @(R) (D_R0 + (2 * R .* tan(gamma_R ./ 2)));
        
        % Calculate Volume of Beam Element
        D1 = @(R) (D_T0 + (2 .* (R - (dL ./ 2)) * tan(gamma_T ./ 2)));
        D2 = @(R) (D_T0 + (2 .* (R + (dL ./ 2)) * tan(gamma_T ./ 2)));
        V_T = @(R) ((tau * dL * (D1(R).^2 + (D1(R) .* D2(R)) + D2(R).^2)) ./ 24);
        
        % Calculate Mass of Spray & Number of Particles in Beam Element
        nParticles = @(R) (nParticlesDensity(R) .* V_T(R));
        mass = @(R) (density(R) .* V_T(R));
        
        % Calculate Local Extinction and Backscattering Coefficients        
        Q_ext = @(R) (Q_ext_Interp(d20(R)));
        Q_bck = @(R) (Q_bck_Interp(d20(R)));

        sigma_ext = @(R) (Q_ext(R) .* ((tau .* d20(R).^2) ./ 8));
        sigma_bck = @(R) (Q_bck(R) .* ((tau .* d20(R).^2) ./ 8));

        alpha = @(R) (alpha0 + (nParticlesDensity(R) .* sigma_ext(R)));
        beta = @(R) (beta0 + (nParticlesDensity(R) .* sigma_bck(R)));
        
        % Calculate Target Impulse Response
        H_O = @(R) ((Gamma_O * (diracDelta((interp1(R_discrete, R_discrete, R, 'nearest') - R_O), 8) / dL)) + beta(R));
%         H_O = @(R) ((Gamma_O .* normc((exp(-((R - R_O) / 1e-12).^2) / (sqrt(tau / 2) * 1e-12)))) + beta(R));

        % Calculate Transmitter/Receiver Overlap Function
        phi_T = @(R) (2 * acos(((D_T(R) ./ 2).^2 - (D_R(R) ./ 2).^2 + axisDisp^2) ./ (2 .* axisDisp .* (D_T(R) ./ 2))));
        phi_R = @(R) (2 * acos(((D_R(R) ./ 2).^2 - (D_T(R) ./ 2).^2 + axisDisp^2) ./ (2 .* axisDisp .* (D_R(R) ./ 2))));

        zeta = @(R) ((((((D_T(R) ./ 2).^2 .* (phi_T(R) - sin(phi_T(R)))) + ...
                     ((D_R(R) ./ 2).^2 .* (phi_R(R) - sin(phi_R(R))))) ./ ...
                    (tau .* (D_T(R) ./ 2).^2)) .* ((R > R1) & (R < R2))) + (1 .* (R >= R2)));

        % Calculate Impulse Response of Optical Channel
        H_C = @(R) (((exp(-sum(alpha(1:R), 'all') .* dL).^2) .* (zeta(R) ./ (tau .* R.^2))));

        % Calculate Impulse Response Function
        H = @(R) (H_C(R) .* H_O(R));

        % Calculate Received Signal
        for k = 1:height(R_discrete)
            t_discrete = (0:dT:((2 * (R_discrete(k) / c)) - dT))';
            
            if length(t_discrete) == 1
                continue;
            end
            
            P_R(k) = trapz(t_discrete, (P_T(t_discrete) .* H(R_discrete(k) - (c * (t_discrete / 2)))));
        end
        clear k;
        
%         % Calculate Received Signal
%         convolution = @(R, t) (P_T(t) .* H(R - (c * (t / 2))));
%         
%         for k = 2:height(R_discrete)
%             timeMax = ((2 * R_discrete(k)) / c); timeMax = timeMax - (timeMax / 1e3);
%             timeR_O = (2 * (R_discrete(k) - R_O)) / c;
%             
%             if timeR_O < 0
%                 P_R(k) = integral(@(t) convolution(R_discrete(k), t), 0, timeMax);
%             elseif timeR_O > 0
%                 P_R(k) = integral(@(t) convolution(R_discrete(k), t), 0, timeR_O) + convolution(R_discrete(k), timeR_O) + integral(@(t) convolution(R_discrete(k), t), (timeR_O + 1e-11), timeMax);
%             else
%                 P_R(k) = integral(@(t) convolution(R_discrete(k), t), 0, timeR_O) + convolution(R_discrete(k), timeR_O);
%             end
%             
%         end
%         clear k;
        
        figTime = num2str(lidarData.time(j), ['%.', num2str(timePrecision), 'f']);
        figTitle = ['{', figTime, ' \it{s}}'];
        
        % Initialise Figure #1
        fig = fig + 1;
        figName = ['Lidar_Signals_Ray_', num2str(RoI), '_T', replace(figTime, '.', '_')];
        set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
                     'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
        pause(0.5);
        hold on;
        set(gca, 'positionConstraint', 'outerPosition', 'plotBoxAspectRatio', [1, 0.75, 0.75], ...
                 'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 22, 'layer', 'top');
             
        % Plot Signals
        plot((R_discrete / normLength), P_R, 'color', graphColours(1), 'lineWidth', 2);
        
        lineHandle = xline((R_O / normLength), 'alpha', 1, ...
                                               'lineStyle', '--', ...
                                               'lineWidth', 2,  ...
                                               'label', 'Target Location', ...
                                               'labelHorizontalAlignment', 'Right', ...
                                               'labelVerticalAlignment', 'Middle');
        lineHandle.Interpreter = 'latex';
        lineHandle.FontSize = 18;
        clear lineHandle;
        
        % Figure Formatting
        title('{-----}', 'interpreter', 'latex');
        subtitle(figTitle);
        axis on;
        box on;
        grid off;
        xlim([0; 5]);
        ylim([0; (1.1 * max(P_R))]);
        tickData = (1:1:4);
        xticks(tickData);
        tickData = [];
        yticks(tickData);
        xtickformat('%.1f');
        xlabel({'{$\ell$}'; '{-----}'}, 'interpreter', 'latex');
        ylabel({'{-----}'; '{Received Signal}'}, 'interpreter', 'latex');
        tightInset = get(gca, 'TightInset');
        set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
                                   (tightInset(2) + 0.00625), ...
                                   (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
                                   (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
        pause(0.5);
        hold off;
        
        % Save Figure
        print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');
        
        
%         % Initialise Figure #2
%         fig = fig + 1;
%         figName = ['Lidar_nParticles_Ray_', num2str(RoI), '_T', replace(figTime, '.', '_')];
%         set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
%                      'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
%         pause(0.5);
%         hold on;
%         set(gca, 'positionConstraint', 'outerPosition', 'plotBoxAspectRatio', [1, 0.75, 0.75], ...
%                  'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 22, 'layer', 'top');
%              
%         % Plot Scattering Particles
%         plot((R_discrete / normLength), nParticles(R_discrete), 'color', graphColours(1), 'lineWidth', 2);
%         
%         % Figure Formatting
%         title('{-----}', 'interpreter', 'latex');
%         subtitle(figTitle);
%         axis on;
%         box on;
%         grid off;
%         xlim([0; 5]);
%         ylim([0; 4.2e3]);
%         tickData = (1:1:4);
%         xticks(tickData);
%         tickData = (0.84e3:0.84e3:3.36e3);
%         yticks(tickData);
%         xtickformat('%.1f');
%         xtickformat('%.1f');
%         axisHandle = gca; axisHandle.YAxis.Exponent = 3; clear axisHandle;
%         xlabel({'{$\ell$}'; '{-----}'}, 'interpreter', 'latex');
%         ylabel({'{-----}'; '{Scattering Particles}'}, 'interpreter', 'latex');
%         tightInset = get(gca, 'TightInset');
%         set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
%                                    (tightInset(2) + 0.00625), ...
%                                    (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
%                                    (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
%         pause(0.5);
%         hold off;
%         
%         % Save Figure
%         print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');
        
        
%         % Initialise Figure #3
%         fig = fig + 1;
%         figName = ['Lidar_d20_Ray_', num2str(RoI), '_T', replace(figTime, '.', '_')];
%         set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
%                      'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
%         pause(0.5);
%         hold on;
%         set(gca, 'positionConstraint', 'outerPosition', 'plotBoxAspectRatio', [1, 0.75, 0.75], ...
%                  'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 22, 'layer', 'top');
%              
%         % Plot Scattering Particles
%         plot((R_discrete / normLength), (d(R_discrete) * 1e6), 'color', graphColours(1), 'lineWidth', 2);
%         
%         % Figure Formatting
%         title('{-----}', 'interpreter', 'latex');
%         subtitle(figTitle);
%         axis on;
%         box on;
%         grid off;
%         xlim([0; 5]);
%         ylim([0; 400]);
%         tickData = (1:1:4);
%         xticks(tickData);
%         tickData = (80:80:320);
%         yticks(tickData);
%         xtickformat('%.1f');
%         xlabel({'{$\ell$}'; '{-----}'}, 'interpreter', 'latex');
%         ylabel({'{-----}'; '{$d_{_{20}}$}'}, 'interpreter', 'latex');
%         tightInset = get(gca, 'TightInset');
%         set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
%                                    (tightInset(2) + 0.00625), ...
%                                    (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
%                                    (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
%         pause(0.5);
%         hold off;
%         
%         % Save Figure
%         print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');


%         % Initialise Figure #5
%         fig = fig + 1;
%         figName = ['Lidar_Mass_Ray_', num2str(RoI), '_T', replace(figTime, '.', '_')];
%         set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
%                      'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
%         pause(0.5);
%         hold on;
%         set(gca, 'positionConstraint', 'outerPosition', 'plotBoxAspectRatio', [1, 0.75, 0.75], ...
%                  'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 22, 'layer', 'top');
%              
%         % Plot Scattering Particles
%         plot((R_discrete / normLength), mass(R_discrete), 'color', graphColours(1), 'lineWidth', 2);
%         
%         % Figure Formatting
%         title('{-----}', 'interpreter', 'latex');
%         subtitle(figTitle);
%         axis on;
%         box on;
%         grid off;
%         xlim([0; 5]);
%         ylim([0; 100e-3]);
%         tickData = (1:1:4);
%         xticks(tickData);
%         tickData = (20e-3:20e-3:80e-3);
%         yticks(tickData);
%         xtickformat('%.1f');
%         xtickformat('%.1f');
%         axisHandle = gca; axisHandle.YAxis.Exponent = -3; clear axisHandle;
%         xlabel({'{$\ell$}'; '{-----}'}, 'interpreter', 'latex');
%         ylabel({'{-----}'; '{Scattering Mass $(kg)$}'}, 'interpreter', 'latex');
%         tightInset = get(gca, 'TightInset');
%         set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
%                                    (tightInset(2) + 0.00625), ...
%                                    (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
%                                    (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
%         pause(0.5);
%         hold off;
%         
%         % Save Figure
%         print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');        
    end
    clear j;
    
end
clear i;


%% Local Functions

function pos = inputPos(orientation)

    pos = str2double(input(['    ', orientation, '-Position [m]: '], 's'));
    
    if isnan(pos) || length(pos) > 1
        disp('        WARNING: Invalid Entry');
        
        pos = -1;
    end
    
end


%         % Test plots
%         fig = fig + 1;
%         figTime = num2str(lidarData.time(j), ['%.', num2str(timePrecision), 'f']);
%         figName = ['Lidar_Detections_Ray_', num2str(RoI), '_T', replace(figTime, '.', '_')];
%         set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
%                      'units', 'pixels', 'outerPosition', [1200, 450, 1500, 1000]);
%         pause(0.5);
%         hold on;
%         set(gca, 'positionConstraint', 'outerPosition', 'plotBoxAspectRatio', [1, 0.75, 0.75], ...
%                  'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 22, 'layer', 'top');
% 
%         tiledlayout(2,2)
%         nexttile(1);
%         hold on;
%         title('H_C');
%         plot(R_discrete, H_C(R_discrete));
%         xlim([0; 20]);
%         hold off;
% 
%         nexttile(2);
%         hold on;
%         title('H_O');
%         plot(R_discrete, H_O(R_discrete));
%         xlim([0; 20]);
%         hold off;
%         
%         nexttile(3);
%         hold on;
%         title('H');
%         plot(R_discrete, H(R_discrete));
%         xlim([0; 20]);
%         hold off;
% 
%         nexttile(4);
%         hold on;
%         title('Received Signal');
%         plot(R_discrete, (P_R / P_O));
%         xline(R_O, 'r--');
%         xlim([0; 20]);
%         hold off;


%         maxP_R = max(maxP_R, max(P_R));
%         
%         % Initialise Figure
%         if j == 1
%             fig = fig + 1;
%             figName = 'Lidar_Signals_Over_Time';
%             set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
%                          'units', 'pixels', 'outerPosition', [50, 50, 795, 880]);
%             pause(0.5);
%             hold on;
%             set(gca, 'positionConstraint', 'outerPosition', 'plotBoxAspectRatio', [1, 0.75, 0.75], ...
%                      'lineWidth', 4, 'fontName', 'LM Mono 12', 'fontSize', 22, 'layer', 'top');
%         end
%              
%         % Plot Signals
%         plot((R_discrete / normLength), P_R, 'color', graphColours(j), 'lineWidth', 2);
%         
%         % Figure Formatting
%         if j == 2
%             lineHandle = xline((R_O / normLength), 'alpha', 1, ...
%                                                    'lineStyle', '--', ...
%                                                    'lineWidth', 2,  ...
%                                                    'label', 'Target Location', ...
%                                                    'labelHorizontalAlignment', 'Right', ...
%                                                    'labelVerticalAlignment', 'Middle');
%             lineHandle.Interpreter = 'latex';
%             lineHandle.FontSize = 18;
%             clear lineHandle;
%             
%             title('{-----}', 'interpreter', 'latex');
%             subtitle('{ }');
%             axis on;
%             box on;
%             grid off;
%             xlim([0; 5]);
%             ylim([0; (1.1 * maxP_R)]);
%             tickData = (1:1:4);
%             xticks(tickData);
%             tickData = [];
%             yticks(tickData);
%             xtickformat('%.1f');
%             xlabel({'{$\ell$}'; '{-----}'}, 'interpreter', 'latex');
%             ylabel({'{-----}'; '{Received Signal}'}, 'interpreter', 'latex');
%             legend({'$10.02\,s$', ...
%                     '$10.04\,s$'}, 'location', 'northEast', ...
%                                    'orientation', 'vertical', ...
%                                    'interpreter', 'latex', ...
%                                    'fontSize', 18, ...
%                                    'box', 'off');
%             tightInset = get(gca, 'TightInset');
%             set(gca, 'innerPosition', [(tightInset(1) + 0.00625), ...
%                                        (tightInset(2) + 0.00625), ...
%                                        (1 - (tightInset(1) + tightInset(3) + 0.0125)), ...
%                                        (1 - (tightInset(2) + tightInset(4) + 0.0125))]);
%             pause(0.5);
%             hold off;
% 
%             % Save Figure
%             print(gcf, [userpath, '/Output/Figures/', figName, '.png'], '-dpng', '-r300');
%         end