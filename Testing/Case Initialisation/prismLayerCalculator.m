%% Prism Layer Calculator v1.1

clearvars;
close all;
clc;

fig = 0; %#ok<*NASGU>
figHold = 0; %#ok<*NASGU>

disp ('===========================');
disp ('Prism Layer Calculator v1.1');
disp ('===========================');
disp (' ');


%% Changelog

% v1.0 - Initial Commit
% v1.1 - Added Basic Error Catching


%% Data Input

variables = {
             'Density'
             'Freestream Velocity [m/s]'
             'Characteristic Length [m]'
             'Dynamic Viscosity [Pa*s]'
            };

valid = false;
while ~valid
    userInput = inputdlg(variables, 'Input Flow Properties', ...
                         [1, 50], {'1.269', '40', '9.44', '1.754e-5'});

    values = str2double(userInput);

    if any(isnan(values))
        disp('WARNING: Invalid Entry (Numeric Values Expected)');
        disp(' ');
    elseif size(values,1) ~= 4
        disp('WARNING: Incomplete Entry');
        disp(' ');
    else
        valid = true;
    end
    
end

rho = values(1);
U = values(2);
L = values(3);
mu = values(4);

disp('Flow Properties');
disp('---------------');
disp(['Density:                       ' num2str(rho) ' [kg/m^3]']);
disp(['Freestream Velocity:           ' num2str(U) ' [m/s]']);
disp(['Characteristic Length:         ' num2str(L) ' [m]']);
disp(['Dynamic Viscosity:             ' num2str(mu) ' [Pa*s]']);
disp(' ');
disp(' ');


%% Prism Layer Property Calculation

variables = {
             'Target y+'
             'Base Cell Size [m]'
             'Expansion Ratio'
            };

valid = false;
while ~valid
    userInput = inputdlg(variables, 'Input Desired Prism Layer Parameters', ...
                         [1, 50], {'30', '0.008', '1.3'});

    values = str2double(userInput);

    if any(isnan(values))
        disp('WARNING: Invalid Entry (Numeric Values Expected)');
        disp(' ');
    elseif size(values,1) ~= 3
        disp('WARNING: Incomplete Entry');
        disp(' ');
    else
        valid = true;
    end
    
end

yPlus = values(1);
base = values(2);
ratio = values(3);

disp('Prism Layer Input');
disp('-----------------');
disp(['Target y+:                     ' num2str(yPlus)]);
disp(['Base Cell Size:                ' num2str(base) ' [m]']);
disp(['Expansion Ratio:               ' num2str(ratio)]);
disp(' ');
disp(' ');


finalTar = base;
Re = (rho * U * L) / mu;
Cf = (2 * log10(Re) - 0.65)^(-2.3);
Tw = Cf * 0.5 * rho * U^2;
uFriction = sqrt(Tw / rho);

firstOptimal = 2 * ((yPlus * mu) / (rho * uFriction));
firstOptimalRel = firstOptimal / base;

cellSize(1) = firstOptimal;
n = 1;
while cellSize(end) < finalTar
    cellSize(n,1) = cellSize(1) * (ratio^(n - 1)); %#ok<SAGROW>
    n = n + 1;
end
cellSize(end) = [];

disp('Prism Layer Output');
disp('------------------');
disp(['Optimal First Cell Thickness {Absolute}: ' num2str(firstOptimal) ' [m]']);
disp(['Optimal First Cell Thickness {Relative}: ' num2str(firstOptimal / base) ' x Base}']);
disp(' ');
disp(['Final Cell Thickness {Absolute}:         ' num2str(cellSize(end)) ' [m]']);
disp(['Final Cell Thickness {Relative}:         ' num2str(cellSize(end) / base) ' x Base}']);
disp(' ');
disp(['Total Prism Thickness {Absolute}:        ' num2str(sum(cellSize(1:end))) ' [m]']);
disp(['Total Prism Thickness {Relative}:        ' num2str(sum(cellSize(1:end)) / base) ' x Base}']);
disp(' ');
disp(['Prism to Core Transition Ratio:          ' num2str(base/cellSize(end-1))]);
disp(' ');
disp(['Required Prism Layers:                   ' num2str(size(cellSize,1))]);


%% Cleaning

clearvars -except base firstOptimal firstOptimalRel cellSize;
disp(' ');