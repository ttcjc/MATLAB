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


%%

load('/mnt/Processing/Data/Numerical/MATLAB/LagData/Distribution_Test/volume/T12_T12_F5.mat');

diameterBins = single((1:147)' * 1e-6); % um


%%

binnedData.d = interp1(diameterBins, diameterBins, LagData.d{1}, 'nearest');

binnedData.mass = zeros([height(diameterBins), 1]);

for i = 1:height(diameterBins)
    index = find(binnedData.d == diameterBins(i));
    
    for j = 1:height(index)
        binnedData.mass(i) = binnedData.mass(i) + (LagData.nParticle{1}(index(j)) * (1000 * ((tau * (LagData.d{1}(index(j))^3)) / 12)));
    end
    
end

% binnedData.massTotal = 0;
% 
% for i = 1:height(LagData.d{1})
%     binnedData.massTotal = binnedData.massTotal + ...
%                            (LagData.nParticle{1}(i) * (1000 * ((tau * (LagData.d{1}(i)^3)) / 12)));
% end


%%

% Initialise Figure
fig = fig + 1;
figName = 'Particle_Size_Distribution';
figTitle = '-';
figSubtitle = ' ';

set(figure(fig), 'name', figName, 'color', [1, 1, 1], ...
                 'outerPosition', [25, 25, 650, 650], 'units', 'pixels')
set(gca, 'positionConstraint', 'outerPosition', ...
         'lineWidth', 2, 'fontName', 'LM Mono 12', 'fontSize', 16, 'layer', 'top', 'yScale', 'log');
hold on;

% Plot Data
plot((diameterBins * 1e6), binnedData.mass, 'lineWidth', 1.5, 'color', ([74, 24, 99] / 255));

% Format Figure
title(figTitle);
subtitle(figSubtitle);
box on;
grid off;
hold off;
xlim([0; 150]);
ylim([1e-9; 1e-4]);
xticks(25:25:125);
yticks([1e-8, 1e-7, 1e-6, 1e-5]);
% xtickformat('%0.2f');
% ytickformat('%0.2f');
xlabel({'{\bf{Particle Diameter [\mum]}}'; '-'}, 'fontName', 'LM Roman 12');
ylabel({'-'; '{\bf{Mass [kg]}}'}, 'fontName', 'LM Roman 12');

% % Save Figure
% pause(1);
% exportgraphics(gcf, [userpath, '/Output/Figures/', figName, '.png'], 'resolution', 600);
% close(fig);