clear variables;
close all;
clc;
evalc('delete(gcp(''nocreate''));');

fig = 0; % Initialise Figure Tracking
figHold = 0; % Enable Overwriting of Figures


%% Load Example Data

load('~/Downloads/Uncertainty/coeffData.mat');


%% Calculate Confidence Interval

Cd = coeffData.Cd_Corr;
% Cd = coeffData.Cd_Corr((find(coeffData.time == 0.5) + 1):end);

Cd = Cd(1:2:end);

N = height(Cd);

Cd_mean = nan(N,1);
Cd_sigma = nan(N,1);
Cd_tVal = nan(N,2);
Cd_CI = nan(N,2);

for i = 2:N
    Cd_mean(i) = mean(Cd(1:i));
    Cd_sigma(i) = std(Cd(1:i));
    Cd_tVal(i,:) = tinv([0.025, 0.975], (i - 1));
    
    Cd_CI(i,:) = Cd_mean(i) + Cd_tVal(i,:) * (Cd_sigma(i) / sqrt(i));
end


%% Visualise Confidence Interval

% Figure Setup
fig = fig + 1;

set(figure(fig), 'outerPosition', [25, 25, 1275, 850], 'name', 'Uncertainty');
set(gca, 'lineWidth', 2, 'fontName', 'LM Mono 12', ...
         'fontSize', 20, 'layer', 'top');
hold on;

% Plot  
plot(Cd, 'lineWidth', 2, 'color', ([74, 24, 99] / 255));
plot(Cd_CI, 'lineWidth', 2, 'color', ([230, 0, 126] / 255));

% Figure Formatting
axis on;
box on;
grid off;
xlim([0; N]);
ylim([0.25; 0.45]);
xticks([]);
yticks([]);
xlabel([]);
xlabel({' ', '\bf{Time}'}, 'fontName', 'LM Roman 12');
ylabel({'\bf{C_D}', ' '}, 'fontName', 'LM Roman 12');
set(gca, 'outerPosition', [0.05, 0.05, 0.9, 0.9]);
hold off;

exportgraphics(gcf, '~/Downloads/Uncertainty/Uncertainty.png', 'resolution', 300);
pause(2);
