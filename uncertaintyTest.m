clear variables;
close all;
clc;
evalc('delete(gcp(''nocreate''));');

fig = 0; % Initialise Figure Tracking
figHold = 0; % Enable Overwriting of Figures


%% Load Data

% load('/mnt/Processing/Data/Numerical/MATLAB/planarContaminantMap/Windsor_SB_wW_Upstream_SC/base/T20025_T50000_F400_D1_D120_Norm.mat');
% load('/mnt/Processing/Data/Numerical/MATLAB/planarContaminantMap/Windsor_SB_wW_Upstream_SC/X_020225/T20025_T50000_F400_D1_D120_Norm.mat');
load('~/Downloads/T20025_T50000_F400_D1_D120_Norm.mat');

% Remove Excess Fields
fields = fieldnames(mapData.mean);

for index = 1:height(fields)
    
    if ~matches(fields{index}, 'mass')
        mapData.inst = rmfield(mapData.inst, fields{index});
        mapData.mean = rmfield(mapData.mean, fields{index});
    end
    
end


%% Calculate Settling Time

% Calculate Settling Time
massTotalAll = cellfun(@sum, mapData.inst.mass);
massTotalAllMovMean = movmean(massTotalAll, 0.05*height(mapData.inst.time));


%% Visualise Settling Time

% Figure Setup
fig = fig + 1;

set(figure(fig), 'outerPosition', [25, 25, 1275, 850], 'name', 'Settling Time');
set(gca, 'lineWidth', 2, 'fontName', 'LM Mono 12', ...
         'fontSize', 20, 'layer', 'top');
hold on;

% Plot  
plot(mapData.inst.time, massTotalAll, 'lineWidth', 2, 'color', ([74, 24, 99] / 255));
plot(mapData.inst.time, massTotalAllMovMean, 'lineWidth', 2, 'color', ([230, 0, 126] / 255));

% Figure Formatting
axis on;
box on;
grid off;
xlim([1.9; 5.1]);
ylim('auto');
xticks(2:0.5:5);
yticks([]);
xlabel({' ', '{\bf{Time (\it{s})}}'}, 'fontName', 'LM Roman 12');
ylabel([]);
set(gca, 'outerPosition', [0.05, 0.05, 0.9, 0.9]);
hold off;

pause(2);


%% Initialise Uncertainty Calculation

% Specify Start Time
valid = false;
while ~valid
    disp(' ');
    startTime = inputTime;
    
    if startTime == -1
        continue;
    end
    
    if startTime < mapData.inst.time(1) || startTime >= mapData.inst.time(end)
        disp('    WARNING: Start Time Must Fall Within Data Range');
    else
        valid = true;
    end
    
end
clear valid;

% Remove Prior Time Instances
index = 1;
while index <= height(mapData.inst.time)

    if mapData.inst.time(index) < startTime
        mapData.inst.time(index) = [];
        mapData.inst.mass(index) = [];
    else
        index = index + 1;
    end

end
clear i;

if mod(height(mapData.inst.time),2) ~= 0
    mapData.inst.time(1) = [];
    mapData.inst.mass(1) = [];
end

%% Uncertainty A

% Calculate Divisors
factors = factor(height(mapData.inst.time));
divisors = [1, factors(1)];

for index = factors(2:end)
    divisors = [1;index] * divisors;
    divisors = unique(divisors(:)');
end

divisors = divisors(2:end)';

% Calculate Moving Means
massTotalSteady = cellfun(@sum, mapData.inst.mass);
massTotalSteadyMovMean = cell(height(divisors),1);

for i = 1:height(divisors)
    massTotalSteadyMovMean{i} = movmean(massTotalSteady, divisors(i), 'endPoints', 'discard');
end

sampleTime = divisors * (1 / 400);
uncertaintyPlot = [];

for i = 1:height(divisors)
    subset = [ones(height(massTotalSteadyMovMean{i}),1) * sampleTime(i), massTotalSteadyMovMean{i}];
    uncertaintyPlot = [uncertaintyPlot; subset]; %#ok<AGROW>
end
clear subset;

% Figure Setup
fig = fig + 1;

set(figure(fig), 'outerPosition', [25, 25, 1275, 850], 'name', 'Uncertainty Convergence A');
set(gca, 'lineWidth', 2, 'fontName', 'LM Mono 12', ...
         'fontSize', 20, 'layer', 'top');
hold on;

% Plot  
scatter(uncertaintyPlot(:,1), uncertaintyPlot(:,2), 10, ([74, 24, 99] / 255), 'filled');

% Figure Formatting
axis on;
box on;
grid off;
xlim([-0.1; 2.6]);
ylim('auto');
xticks(0:0.5:2.5);
yticks([]);
xlabel({' ', '{\bf{Sample Length (\it{s})}}'}, 'fontName', 'LM Roman 12');
ylabel([]);
set(gca, 'outerPosition', [0.05, 0.05, 0.9, 0.9]);
hold off;

pause(2);


%% Uncertainty B

% Calculate 95% Confidence Intervals
subsetMean = nan(height(mapData.inst.time),1);
subsetSigma = nan(height(mapData.inst.time),1);
subsetCI = nan(height(mapData.inst.time),2);

for i = 2:height(mapData.inst.time)
    subsetMean(i) = mean(massTotalSteady(1:i));
    subsetSigma(i) = std(massTotalSteady(1:i));

    subsetCI(i,1) = subsetMean(i) - (1.96 * (subsetSigma(i) / sqrt(i)));
    subsetCI(i,2) = subsetMean(i) + (1.96 * (subsetSigma(i) / sqrt(i)));
end

sampleTime = (1:height(mapData.inst.time))' * (1 / 400);

% Figure Setup
fig = fig + 1;

set(figure(fig), 'outerPosition', [25, 25, 1275, 850], 'name', 'Uncertainty Convergence B');
set(gca, 'lineWidth', 2, 'fontName', 'LM Mono 12', ...
         'fontSize', 20, 'layer', 'top');
hold on;

% Plot  
plot(sampleTime, subsetCI(:,1), 'lineWidth', 2, 'color', ([230, 0, 126] / 255));
plot(sampleTime, subsetMean(:,1), 'lineWidth', 2, 'color', ([74, 24, 99] / 255));
plot(sampleTime, subsetCI(:,2), 'lineWidth', 2, 'color', ([230, 0, 126] / 255));

% Figure Formatting
axis on;
box on;
grid off;
xlim([-0.1; 2.6]);
ylim('auto');
xticks(0:0.5:2.5);
yticks([]);
xlabel({' ', '{\bf{Sample Length (\it{s})}}'}, 'fontName', 'LM Roman 12');
ylabel([]);
set(gca, 'outerPosition', [0.05, 0.05, 0.9, 0.9]);
hold off;

pause(2);


%% Local Functions

function time = inputTime

    time = str2double(input('Input Desired Averaging Start Time [s]: ', 's'));
    
    if isnan(time) || length(time) > 1 || time <= 0
        disp('        WARNING: Invalid Entry');
        time = -1;
    end

end