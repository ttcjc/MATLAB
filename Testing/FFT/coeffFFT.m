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
figHold = 0; % Enable Overwriting of Figures


%% Load Data      

% load('~/MATLAB/Testing/FFT/coeffDataQS.mat');
load('~/MATLAB/Testing/FFT/coeffDataFS.mat');


%% Perform FFT

% sampleStart = find(coeffData.time == 0.25) + 1;
sampleStart = find(coeffData.time == 0.75) + 1;

sampleLength = height(coeffData.time(sampleStart:end));

windowSize = sampleLength / 1;

nDFT = 2^nextpow2(sampleLength);

% sampleFreq = 50000;
sampleFreq = 5000;

[PSD, freq] = pwelch(coeffData.Cs(sampleStart:end), round(windowSize), round((windowSize / 2)), nDFT, sampleFreq);

% Sr = (freq * 0.289) / 40;
Sr = (freq * 1.156) / 22.22;


%% Plot Results

% Initialise Figure
fig = fig + 1;
figName = 'FFT';
set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
                 'outerPosition', [25, 25, 650, 650], 'units', 'pixels')
set(gca, 'positionConstraint', 'outerPosition', ...
                 'lineWidth', 2, 'fontName', 'LM Mono 12', 'fontSize', 16, 'layer', 'top', ...
                 'xScale', 'log');
hold on;

% Plot Data
plot(freq, (freq .* PSD), 'lineWidth', 1.5, 'color', ([74, 24, 99] / 255));

% Format Figure
box on;
grid off;
hold off;
xlim('padded');
ylim('padded');
xticks('auto');
yticks('auto');
xtickformat('auto');
ytickformat('auto');
xlabel({'{\bf{Frequency}}'}, 'fontName', 'LM Roman 12');
ylabel({'{\bf{f\cdot\Phi}}'}, 'fontName', 'LM Roman 12');

% Save Figure
pause(2);
% exportgraphics(gcf, [userpath, '/Output/Figures/', figName, '.png'], 'resolution', 600);