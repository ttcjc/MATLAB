%% Distribution Discretiser v1.0

clearvars;
close all;
clc;

fig = 0; %#ok<*NASGU>
figHold = 0; %#ok<*NASGU>

disp ('=============================');
disp ('Distribution Discretiser v1.0');
disp ('=============================');
disp (' ');


%% Changelog

% v1.0 - Initial Commit


%% Distribution Comparison

distOrig = load('~/Documents/Engineering/PhD/Data/Numerical/MATLAB/Particle Distributions/Anton_Developed_Spray.mat');
data = load('~/Documents/Engineering/PhD/Data/Numerical/MATLAB/particleData/Lag_Test_Mass_25.mat');
% data = load('~/Documents/Engineering/PhD/Data/Numerical/MATLAB/particleData/Lag_Test_Freq.mat');

% Discretise Diameters
binSize = 1e-6; % m

distOut.diameter = (min(round(data.particleData.d{1,1} / binSize) * binSize):binSize:max(round(data.particleData.d{1,1} / binSize) * binSize))';


% Adjust for Small Particles (D < 1 um)
i = 1;
while distOut.diameter(i) == 0
	distOut.diameter(i) = binSize;
	i = i + 1;
end

% Discretise Population
assignments = discretize(data.particleData.d{1,1}, distOut.diameter);

binTotals = zeros(size(distOut.diameter,1),1);
nParticle = zeros(size(distOut.diameter,1),1);

for i = 1:size(distOut.diameter,1)
	binTotals(i) = sum(assignments == i);
	nParticle(i) = sum(data.particleData.nParticle{1,1}(assignments == i));
end

parcelProb = (binTotals / sum(binTotals)) * 100;
particleProb = (nParticle / sum(nParticle)) * 100;

% Discretise Mass-Weighted Population
mass = zeros(size(distOut.diameter,1),1);

% for i = 1:size(distOut.diameter,1)
% 	mass(i) = binTotals(i) * (1000 * (4 / 3) * pi * (distOut.diameter(i) / 2)^3);
% end
% 
% massProb = (mass / sum(mass)) * 100;

for i = 1:size(distOut.diameter,1)
	mass(i) = nParticle(i) * (1000 * (4 / 3) * pi * (distOut.diameter(i) / 2)^3);
end

massProb = (mass / sum(mass)) * 100;

% Figure Setup
fig = fig + 1;
figure('name', 'Relative Particle Size Distribution');
hold on;
set(figure(fig), 'outerPosition', [25, 25, 800, 800]);

% Plot
plot(distOrig.discrete.diameter, distOrig.discrete.populationByMass, 'r');
plot(distOrig.discrete.diameter, distOrig.discrete.population, 'b');

plot((distOut.diameter * 1e6), massProb, 'g--');
plot((distOut.diameter * 1e6), parcelProb, 'k--');
plot((distOut.diameter * 1e6), particleProb, 'm--');

% Figure Formatting
xlabel({' ', 'Particle Diameter (\it{\mum})'});
ylabel({'Population (%)', ' '});
grid on;
box on;
legend('Original (Mass-Based)', 'Original (Particle-Based)', 'Output (Mass Dist)', 'Output (Parcel Dist)', 'Output (Particle Dist)', 'location', 'northOutside', 'orientation', 'vertical');
legend boxoff;
set(gca, 'units', 'normalized', 'position', [0.1275, 0.1275, 0.745, 0.745], ...
	     'fontName', 'LM Roman 12', 'fontSize', 12, 'layer', 'top');