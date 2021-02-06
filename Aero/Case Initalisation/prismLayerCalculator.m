%% Prism Layer Calculator v2.0

clear variables;
close all;
clc;

disp ('===========================');
disp ('Prism Layer Calculator v2.0');
disp ('===========================');
disp (' ');


%% Changelog

% v1.0 - Initial Commit
% v2.0 - Layout Restructure


%% Flow Properties

data = {'Density', 'Freestream Velocity [m/s]', ...
        'Characteristic Length [m]', 'Dynamic Viscosity [Pa*s]'};
dlg_title = 'Flow Properties';
default = {'1.269', '40', '9.408', '1.754e-5'};
resize = 'on';
values = inputdlg(data, dlg_title, [1,60], default, resize);
input = str2double(values);

rho = input(1);
U = input(2);
L = input(3);
mu = input(4);

disp('Flow Properties');
disp('---------------');
disp(['Density:                       ' num2str(rho) ' [kg/m^3]']);
disp(['Freestream Velocity:           ' num2str(U) ' [m/s]']);
disp(['Characteristic Length:         ' num2str(L) ' [m]']);
disp(['Dynamic Viscosity:             ' num2str(mu) ' [Pa*s]']);
disp(' ');
disp(' ');


%% Prism Layer Calculator

data = {'Target y+', 'Base Cell Size [m]', 'Expansion Ratio', ...
        'Target Final Cell Thickness {Relative to Base}'};
dlg_title = 'Prism Layer Calculator';
default = {'30', '0.008', '1.3', '1'};
resize = 'on';
values = inputdlg(data, dlg_title, [1,60], default, resize);
input = str2double(values);

yPlus = input(1);
base = input(2);
ratio = input(3);
finalTar = input(4)*base;

disp('Prism Layer Input');
disp('-----------------');
disp(['Target y+:                     ' num2str(yPlus)]);
disp(['Base Cell Size:                ' num2str(base) ' [m]']);
disp(['Expansion Ratio:               ' num2str(ratio)]);
disp(' ');
disp(' ');

Re = (rho*U*L)/mu;
Cf = (2*log10(Re)-0.65)^(-2.3);
Tw = Cf*0.5*rho*U^2;
uFriction = sqrt(Tw/rho);
firstOptimal = 2*((yPlus*mu)/(rho*uFriction));

cellSize(1) = firstOptimal;
n = 1;
while cellSize(end) < finalTar
    cellSize(n,1) = cellSize(1)*(ratio^(n-1));
    n = n+1;
end
cellSize(end) = [];

disp('Prism Layer Output');
disp('------------------');
disp(['Optimal First Cell Thickness {Absolute}: ' num2str(firstOptimal) ' [m]']);
disp(['Optimal First Cell Thickness {Relative}: ' num2str(firstOptimal/base) ' x Base}']);
disp(' ');
disp(['Final Cell Thickness {Absolute}:         ' num2str(cellSize(end-1)) ' [m]']);
disp(['Final Cell Thickness {Relative}:         ' num2str(cellSize(end-1)/base) ' x Base}']);
disp(' ');
disp(['Total Prism Thickness {Absolute}:        ' num2str(sum(cellSize(1:end-1))) ' [m]']);
disp(['Total Prism Thickness {Relative}:        ' num2str(sum(cellSize(1:end-1))/base) ' x Base}']);
disp(' ');
disp(['Prism to Core Transition Ratio:          ' num2str(base/cellSize(end-1))]);
disp(' ');
disp(['Required Prism Layers:                   ' num2str(size(cellSize,1)-1)]);


%% Clean-up

clear data default dlg_title input n resize values;