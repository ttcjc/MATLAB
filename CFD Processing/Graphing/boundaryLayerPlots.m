%% Boundary Layer Probing v1.0

clearvars;
close all;
clc;

fig = 0;
figHold = 0; %#ok<*NASGU>

disp ('===========================');
disp ('Boundary Layer Probing v1.0');
disp ('===========================');
disp (' ');


%% Changelog

% v1.0 - Initial Commit


%% Case Selection

disp('CASE SELECTION');
disp('--------------');
disp(' ');

disp('Select Case Folder:');
% caseFolder = uigetdir('~/OpenFOAM');
caseFolder = ('~/Mount/Uni/OpenFOAM/ttcjc-7/results/BL_Test/BL_Test_9.408_SA');
disp(['Case: ', caseFolder]);
disp(' ');
[timeDirs, deltaT] = timeDirectories(caseFolder);
disp(' ');
disp(' ');

if ~exist([caseFolder, '/postProcessing/probes'], 'dir')
    error('No Probe Data Found for Target Case');
end


%% Location A (x = -0.537 l)

disp('Location A (x = -0.537 l)');
disp('-------------------------');
disp(' ');

% Figure Setup
fig = fig + 1;
figure('name', 'Location A');
hold on;
set(figure(fig), 'outerPosition', [1945, 25, 750, 750]);

% Experimental Data (Luckhurst)
probeData = csvread('~/Documents/Engineering/PhD/Literature/Loughborough Wind Tunnel/Boundary Layer Measurements/Luckhurst_Exp_Location_A.csv');
u = probeData(:,1) * 40;
z = probeData(:,2) / (1e3 * 1.044);

U = max(u);
index = find(u == U, 1, 'first');
u = u(1:index);
z = z(1:index);

scatter(u/U, z, 10, [0.21176 0.06667 0.38824]);

% Numerical Data (Luckhurst)
probeData = csvread('~/Documents/Engineering/PhD/Literature/Loughborough Wind Tunnel/Boundary Layer Measurements/Luckhurst_CFD_Location_A.csv');
u = probeData(:,1) * 40;
z = probeData(:,2) / (1e3 * 1.044);

U = max(u);
index = find(u == U, 1, 'first');
u = u(1:index);
z = z(1:index);

LuckhurstA = table(z, u, u/U, 'variableNames', {'z (l)', 'u [m/s]', 'u/U'});
delta = interp1(u/U, z, 0.99);
deltaStar = trapz(z, 1 - (u / U));
theta = trapz(z, (u / U) .* ( 1 -(u / U)));

plot(u/U, z, 'color', [0.71765 0.00000 0.38431]);

disp('Luckhurst (Numerical)');
disp(['    Boundary Layer Thickness (', char(948), ') ~ ', num2str(delta), ' (l)']);
disp(['    Boundary Layer Displacement Thickness (', char(948), '*) ~ ', num2str(deltaStar), ' (l)']);
disp(['    Boundary Layer Momentum Thickness (', char(952), ') ~ ', num2str(theta), ' (l)']);
disp(' ');

% Numerical Data (Crickmore)
fileID = fopen([caseFolder, '/postProcessing/probes/', timeDirs(end,1).name, '/boundaryLayerLocationA_U.xy']);
probeData = textscan(fileID, '%f %f %f %f', 'headerLines', 1, 'delimiter', ' ', 'MultipleDelimsAsOne', 1);
u = vertcat(0, probeData{1,2});
z = vertcat(0, probeData{1,1} / 1.044);

[~,index] = unique(u);
z = z(index);
u = u(index);

U = max(u);
index = find(u == U, 1, 'first');
u = u(1:index);
z = z(1:index);

CrickmoreA = table(z, u, u/U, 'variableNames', {'z (l)', 'u [m/s]', 'u/U'});
delta = interp1(u/U, z, 0.99);
deltaStar = trapz(z, 1 - (u / U));
theta = trapz(z, (u / U) .* ( 1 -(u / U)));

plot(u/U, z, 'color', [0.94902 0.41569 0.21961]);

disp('Crickmore (Numerical)');
disp(['    Boundary Layer Thickness (', char(948), ') ~ ', num2str(delta), ' (l)']);
disp(['    Boundary Layer Displacement Thickness (', char(948), '*) ~ ', num2str(deltaStar), ' (l)']);
disp(['    Boundary Layer Momentum Thickness (', char(952), ') ~ ', num2str(theta), ' (l)']);
disp(' ');
disp(' ');

% Figure Formatting
xlim([0, 1.1]);
ylim([0, 0.08]);
xticks(0.1:0.1:1);
yticks(0.005:0.01:0.075);
xtickformat('%.1f');
ytickformat('%.3f');
xlabel({' ', '(\it{U})'});
ylabel({'z (\it{l})', ' '});
grid on;
box on;
legend('Experimental (Luckhurst)', 'Numerical (Luckhurst)', 'Working Section = 9.011 \it{l}', 'location', 'northOutside', 'orientation', 'vertical');
legend boxoff;
set(gca, 'units', 'normalized', 'position', [0.125, 0.125, 0.75, 0.75], ...
         'fontName', 'LM Roman 12', 'fontSize', 11, 'layer', 'top');
hold off;

print(fig, '~/MATLAB/Output/Figures/Boundary_Probe_A', '-dpng', '-r300');

clearvars -except fig figHold caseFolder timeDirs;


%% Location B (X = 0)

disp('Location B (x = 0 l)');
disp('--------------------');
disp(' ');

% Figure Setup
fig = fig + 1;
figure('name', 'Location B');
hold on;
set(figure(fig), 'outerPosition', [1945, 100, 750, 750]);

% Experimental Data (Luckhurst)
probeData = csvread('~/Documents/Engineering/PhD/Literature/Loughborough Wind Tunnel/Boundary Layer Measurements/Luckhurst_Exp_Location_B.csv');
u = probeData(:,1) * 40;
z = probeData(:,2) / (1e3 * 1.044);

U = max(u);
index = find(u == U, 1, 'first');
u = u(1:index);
z = z(1:index);

scatter(u/U, z, 10, [0.21176 0.06667 0.38824]);

% Numerical Data (Luckhurst)
probeData = csvread('~/Documents/Engineering/PhD/Literature/Loughborough Wind Tunnel/Boundary Layer Measurements/Luckhurst_CFD_Location_B.csv');
u = probeData(:,1) * 40;
z = probeData(:,2);

U = max(u);
index = find(u == U, 1, 'first');
u = u(1:index);
z = z(1:index) / (1e3 * 1.044);

LuckhurstB = table(z, u, u/U, 'variableNames', {'z (l)', 'u [m/s]', 'u/U'});
delta = interp1(u/U, z, 0.99);
deltaStar = trapz(z, 1 - (u / U));
theta = trapz(z, (u / U) .* ( 1 -(u / U)));

plot(u/U, z, 'color', [0.71765 0.00000 0.38431]);

disp('Luckhurst (Numerical)');
disp(['    Boundary Layer Thickness (', char(948), ') ~ ', num2str(delta), ' (l)']);
disp(['    Boundary Layer Displacement Thickness (', char(948), '*) ~ ', num2str(deltaStar), ' (l)']);
disp(['    Boundary Layer Momentum Thickness (', char(952), ') ~ ', num2str(theta), ' (l)']);
disp(' ');

% Numerical Data (Crickmore)
fileID = fopen([caseFolder, '/postProcessing/probes/', timeDirs(end,1).name, '/boundaryLayerLocationB_U.xy']);
probeData = textscan(fileID, '%f %f %f %f', 'headerLines', 1, 'delimiter', ' ', 'MultipleDelimsAsOne', 1);
u = vertcat(0, probeData{1,2});
z = vertcat(0, probeData{1,1} /1.044);

[~,index] = unique(u);
z = z(index);
u = u(index);

U = max(u);
index = find(u == U, 1, 'first');
u = u(1:index);
z = z(1:index);

CrickmoreB = table(z, u, u/U, 'variableNames', {'z (l)', 'u [m/s]', 'u/U'});
delta = interp1(u/U, z, 0.99);
deltaStar = trapz(z, 1 - (u / U));
theta = trapz(z, (u / U) .* ( 1 -(u / U)));

plot(u/U, z, 'color', [0.94902 0.41569 0.21961]);

disp('Crickmore (Numerical)');
disp(['    Boundary Layer Thickness (', char(948), ') ~ ', num2str(delta), ' (l)']);
disp(['    Boundary Layer Displacement Thickness (', char(948), '*) ~ ', num2str(deltaStar), ' (l)']);
disp(['    Boundary Layer Momentum Thickness (', char(952), ') ~ ', num2str(theta), ' (l)']);
disp(' ');
disp(' ');

% Figure Formatting
xlim([0, 1.1]);
ylim([0, 0.08]);
xticks(0.1:0.1:1);
yticks(0.005:0.01:0.075);
xtickformat('%.1f');
ytickformat('%.3f');
xlabel({' ', '(\it{U})'});
ylabel({'z (\it{l})', ' '});
grid on;
box on;
legend('Experimental (Luckhurst)', 'Numerical (Luckhurst)', 'Working Section = 9.011 \it{l}', 'location', 'northOutside', 'orientation', 'vertical');
legend boxoff;
set(gca, 'units', 'normalized', 'position', [0.125, 0.125, 0.75, 0.75], ...
         'fontName', 'LM Roman 12', 'fontSize', 11, 'layer', 'top');
hold off;

print(fig, '~/MATLAB/Output/Figures/Boundary_Probe_B', '-dpng', '-r300');

clearvars;