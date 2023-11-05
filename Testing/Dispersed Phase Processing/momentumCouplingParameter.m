%% Preamble

clear variables;
close all;
clc;
evalc('delete(gcp(''nocreate''));');

if exist('/mnt/Processing/Data', 'dir')
    saveLocation = '/mnt/Processing/Data';
else
    saveLocation = '~/Data';
end

nProc = maxNumCompThreads - 2; % Number of Processors Used for Parallelisation

fig = 0; % Initialise Figure Tracking
figHold = 0; %#ok<NASGU> % Enable Overwriting of Figures


%% Initialise

evalc('parpool(''threads'');');

load('~/Downloads/Windsor_FS_Geometry.mat');
load('~/Downloads/Windsor_FS_MP_Test.mat');

nTimes = height(LagData.time);

dLims = [20e-6; 400e-6];

xLimsData = 4 * [-0.56075; 1.52725];
yLimsData = 4 * [-0.3915; 0.3915];
zLimsData = 4 * [0; 0.4698];


%% Collate Particles of Interest

index = cell(nTimes,1);

d = LagData.d;
positionCartesian = LagData.positionCartesian;
parfor i = 1:nTimes
    
    if ~isempty(positionCartesian{i})
        index{i} = find(((d{i} * 1e6) >= dLims(1)) & ...
                        ((d{i} * 1e6) <= dLims(2)) & ...
                        (positionCartesian{i}(:,1) >= xLimsData(1)) & ...
                        (positionCartesian{i}(:,1) <= xLimsData(2)) & ...
                        (positionCartesian{i}(:,2) >= yLimsData(1)) & ...
                        (positionCartesian{i}(:,2) <= yLimsData(2)) & ...
                        (positionCartesian{i}(:,3) >= zLimsData(1)) & ...
                        (positionCartesian{i}(:,3) <= zLimsData(2))); %#ok<PFBNS>
    end
    
    % Remove Unnecessary Data
    d{i} = [];
    positionCartesian{i} = [];
end
clear i d positionCartesian;


%% Generate Presentation Grid

cellSize.target = 32e-3;

cellSize.x = (xLimsData(2) - xLimsData(1)) / round(((xLimsData(2) - xLimsData(1)) / cellSize.target));
cellSize.y = (yLimsData(2) - yLimsData(1)) / round(((yLimsData(2) - yLimsData(1)) / cellSize.target));
cellSize.z = (zLimsData(2) - zLimsData(1)) / round(((zLimsData(2) - zLimsData(1)) / cellSize.target));

cellSize.volume = cellSize.x * cellSize.y * cellSize.z;

[x, y, z] = ndgrid(xLimsData(1):(cellSize.x):xLimsData(2), ...
                   yLimsData(1):(cellSize.y):yLimsData(2), ...
                   zLimsData(1):(cellSize.z):zLimsData(2));


%% Assign Particles to Mesh Nodes

totalParticles = cellfun(@height, LagData.positionCartesian);
index = cell(nTimes,1); % Array Position of Closest Mesh Node

positionCartesian = LagData.positionCartesian;
parfor i = 1:nTimes
    
    if totalParticles(i) > 0
        index{i} = zeros(totalParticles(i),3);
        
        for j = 1:totalParticles(i)
            [~, index{i}(j,1)] = min(abs(positionCartesian{i}(j,1) - x(:,1,1))); %#ok<PFBNS>
            [~, index{i}(j,2)] = min(abs(positionCartesian{i}(j,2) - y(1,:,1))); %#ok<PFBNS>
            [~, index{i}(j,3)] = min(abs(positionCartesian{i}(j,3) - z(1,1,:))); %#ok<PFBNS>
        end
    
    end
    
    % Remove Unnecessary Data
    positionCartesian{i} = [];
end
clear i positionCartesian;


%% Calculate Instantaneous Field Variables

volumeData.time = LagData.time;

nParticles = cell(nTimes,1); nParticles(:) = {zeros(size(x), 'single')}; % Number of Particles in Cell
volFraction = nParticles; % Fraction of Cell Volume Occupied by Spray
mass = nParticles; % Total Mass in Cell
d10 = nParticles; % Arithmetic Mean Diameter in Cell
Up = nParticles; % Mean Particle Velocity Magnitude in Cell
Uf = nParticles; % Mean Fluid Velocity Magnitude in Cell

nParticle = LagData.nParticle;
d = LagData.d;
U = LagData.U;
Uslip = LagData.Uslip;
cellVolume = cellSize.volume;
parfor i = 1:nTimes
    
    if totalParticles(i) > 0

        for j = 1:totalParticles(i)      
            nParticles{i}(index{i}(j,1), index{i}(j,2), index{i}(j,3)) = ...
            nParticles{i}(index{i}(j,1), index{i}(j,2), index{i}(j,3)) + ...
            nParticle{i}(j);
            
            mass{i}(index{i}(j,1), index{i}(j,2), index{i}(j,3)) = ...
            mass{i}(index{i}(j,1), index{i}(j,2), index{i}(j,3)) + ...
            (nParticle{i}(j) * ((1 / 12) * tau * (d{i}(j)^3)));
            
            d10{i}(index{i}(j,1), index{i}(j,2), index{i}(j,3)) = ...
            d10{i}(index{i}(j,1), index{i}(j,2), index{i}(j,3)) + ...
            (nParticle{i}(j) * d{i}(j));
            
%             Up{i}(index{i}(j,1), index{i}(j,2), index{i}(j,3)) = ...
%             Up{i}(index{i}(j,1), index{i}(j,2), index{i}(j,3)) + ...
%             (nParticle{i}(j) * sqrt(U{i}(j,1)^2 + U{i}(j,2)^2 + U{i}(j,3)^2));
            
            Up{i}(index{i}(j,1), index{i}(j,2), index{i}(j,3)) = ...
            Up{i}(index{i}(j,1), index{i}(j,2), index{i}(j,3)) + ...
            (nParticle{i}(j) * U{i}(j,1));
            
%             Uf{i}(index{i}(j,1), index{i}(j,2), index{i}(j,3)) = ...
%             Uf{i}(index{i}(j,1), index{i}(j,2), index{i}(j,3)) + ...
%             (nParticle{i}(j) * (sqrt(Uslip{i}(j,1)^2 + Uslip{i}(j,2)^2 + Uslip{i}(j,3)^2) - sqrt(U{i}(j,1)^2 + U{i}(j,2)^2 + U{i}(j,3)^2)));
            
            Uf{i}(index{i}(j,1), index{i}(j,2), index{i}(j,3)) = ...
            Uf{i}(index{i}(j,1), index{i}(j,2), index{i}(j,3)) + ...
            (nParticle{i}(j) * (Uslip{i}(j,1) - U{i}(j,1)));
        end

    end

    % Remove Unnecessary Data
    index{i} = [];
    nParticle{i} = [];
    d{i} = [];

    % Calculate Derived Variables
    volFraction{i} = mass{i} / cellVolume;
    mass{i} = 1000 * mass{i};
    d10{i} = (d10{i} ./ nParticles{i});
    Up{i} = abs((Up{i} ./ nParticles{i}));
    Uf{i} = abs((Uf{i} ./ nParticles{i}));
end
clear i j nParticle d cellVolume;

clear LagData;

volumeData.inst.volFraction = volFraction; clear volFraction;
volumeData.inst.mass = mass; clear mass;
volumeData.inst.d10 = d10; clear d10;
volumeData.inst.Up = Up; clear Up;
volumeData.inst.Uf = Uf; clear Uf;

volumeData.inst.Stk = cell(nTimes,1); volumeData.inst.Stk(:) = {zeros(size(x), 'single')}; % Stokes Number  in Cell
volumeData.inst.PI = volumeData.inst.Stk;

for i = 1:nTimes
    volumeData.inst.Stk{i} = (1000 * volumeData.inst.d10{i}.^2 * 22.22) / (18 * (1.754e-5) * 4.176);
    
    mDotWater = (volumeData.inst.mass{i} / cellSize.x) .* volumeData.inst.Up{i};
    mDotAir = 1.269 * (cellSize.y * cellSize.z) * volumeData.inst.Uf{i};
    
    volumeData.inst.PI{i} = mDotWater ./ (mDotAir .* (1 + volumeData.inst.Stk{i}));
end



%% Present Volume Field

xInit = x;
yInit = y;
zInit = z;
POD = false;

cMap = viridis(3); cMap = cMap(1,:);

figTitle = '-'; % Leave Blank ('-') for Formatting Purposes
    
xLimsPlot = xLimsData; yLimsPlot = yLimsData; zLimsPlot = zLimsData;
    
figHold = fig;

isoValue = 1;

for i = 1:nTimes

    if i ~= 1
        clf(fig);
        fig = figHold;
    end

    fieldData = volumeData.inst.PI{i};
    figName = 'Test';
    figSubtitle = ' ';

    fig = plotVolumeField(xLimsData, yLimsData, zLimsData, xInit, yInit, zInit, POD, fieldData, ...
                          fig, figName, geometry, isoValue, cMap, figTitle, figSubtitle, ...
                          xLimsPlot, yLimsPlot, zLimsPlot, false);
end
clear i;