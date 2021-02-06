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


%% Discretiser

% Import Continuous Spray Distribution
data = importdata('~/Downloads/Anton/Developed Spray Distribution.csv');
dataName = 'Anton_Developed_Spray';

continuous.diameter = data(:,1);
continuous.mass = data(:,2);

% Discretise Diameters
binSize = 1; % um

discrete.diameter = (min(round(continuous.diameter / binSize) * binSize):binSize:max(round(continuous.diameter / binSize) * binSize))';

% Adjust for Small Particles (D < 1 um)
i = 1;
while discrete.diameter(i) == 0
	discrete.diameter(i) = binSize;
	i = i + 1;
end

% Produce Equal Spacing
interpSize = 0.25;

interpA = (0:interpSize:max(discrete.diameter) + 0.5)';
[interpB, index, ~] = unique(continuous.diameter);
interpC = continuous.mass(index);
interpD = interp1(interpB, interpC, interpA, 'spline');

% Discretise Population
discrete.mass = zeros(size(discrete.diameter,1),1);

for i = 1:max(discrete.diameter)
	index = find(interpA == i);
	a = interpD(index - (0.5 / interpSize));
	b = interpD(index + (0.5 / interpSize));
	discrete.mass(i) = ((a + b) / 2) * 1;
end

% Figure Setup
fig = fig + 1;
figure('name', 'Absolute Particle Size Distribution');
hold on;
set(figure(fig), 'outerPosition', [25, 25, 800, 800]);

% Plot
plot(continuous.diameter, continuous.mass, 'color', [0.21176, 0.06667, 0.38824]);
plot(interpA, interpD, 'color', [0.71765, 0.00000, 0.38431]);
plot(discrete.diameter, discrete.mass, 'color', [0.94902 0.41569 0.21961]);

% Figure Formatting
xlabel({' ', 'Particle Diameter (\it{\mum})'});
ylabel({'Mass (\it{kg})', ' '});
grid on;
box on;
legend('Continuous Data', 'Continuous Interpolation', 'Discretised Data', 'location', 'northOutside', 'orientation', 'vertical');
legend boxoff;
set(gca, 'units', 'normalized', 'position', [0.1275, 0.1275, 0.745, 0.745], ...
	     'fontName', 'LM Roman 12', 'fontSize', 12, 'layer', 'top');
hold off;

discrete.populationByMass = (discrete.mass ./ sum(discrete.mass)) * 100;

discrete.frequency = zeros(size(discrete.diameter,1),1);

for i = 1:size(discrete.diameter,1)
	discrete.frequency(i) = discrete.mass(i) / (1000 * ((4 / 3) * pi * ((discrete.diameter(i) / 2) / 1e6)^3));
end

discrete.population = (discrete.frequency ./ sum(discrete.frequency)) * 100;

% Figure Setup
fig = fig + 1;
figure('name', 'Relative Particle Size Distribution');
hold on;
set(figure(fig), 'outerPosition', [25, 25, 800, 800]);

% Plot
plot(discrete.diameter, discrete.populationByMass, 'color', [0.21176, 0.06667, 0.38824]);
plot(discrete.diameter, discrete.population, 'color', [0.71765, 0.00000, 0.38431]);

% Figure Formatting
xlabel({' ', 'Particle Diameter (\it{\mum})'});
ylabel({'Population (%)', ' '});
grid on;
box on;
legend('Mass-Based', 'Particle-Based', 'location', 'northOutside', 'orientation', 'vertical');
legend boxoff;
set(gca, 'units', 'normalized', 'position', [0.1275, 0.1275, 0.745, 0.745], ...
	     'fontName', 'LM Roman 12', 'fontSize', 12, 'layer', 'top');

% Save Data 
save(['~/Documents/Engineering/PhD/Data/Numerical/MATLAB/Particle Distributions/', dataName, '.mat'], 'continuous', 'discrete');


%% Cleaning

clearvars -except continuous discrete
disp(' ');