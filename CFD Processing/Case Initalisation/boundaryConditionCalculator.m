%% Turbulence Boundary Condition Calculator v1.0

clearvars;
close all;
clc;

disp ('=============================================');
disp ('Turbulence Boundary Condition Calculator v1.0');
disp ('=============================================');
disp (' ');


%% Changelog

% v1.0 - Initial Commit


%% Data Input

data = {'Characteristic Length [m]', 'Turbulence Intensity [%]', ...
        'Freestream Velocity [m/s]'};
dlg_title = 'Turbulence Parameters';
default = {'1.044', '0.2', '40'};
resize = 'on';
values = inputdlg(data, dlg_title, [1,60], default, resize);
input = str2double(values);

L = input(1);
I = input(2)/100;
U = input(3);

disp(['Characteristic Length:         ' num2str(L) ' [m]']);
disp(['Turbulence Intensity:          ' num2str(I*100) ' [%]']);
disp(['Freestream Velocity:           ' num2str(U) ' [m/s]']);
disp(' ')


%% Prism Layer Calculator

Lt = 0.07 * L;

kappa = 1.5 * (U^2) * (I^2);
epsilon = 0.09 * ((kappa^(1.5)) / Lt);
omega = 0.09^(-0.25) * (sqrt(kappa) / Lt);
nuTilda = 0.09 * ((kappa^2) / epsilon);

disp([char(954), ' {Freestream} = ', num2str(kappa), ' [m^2 s^-2]']);
disp([char(949), ' {Freestream} = ', num2str(epsilon), ' [m^2 s^-1]']);
disp([char(969), ' {Freestream} = ', num2str(omega), ' [s^-1]']);
disp([char(957), ' {Freestream} = ', num2str(nuTilda), ' [m^2 s^-1]']);


%% Cleaning

clearvars -except kappa epsilon omega nuTilda;
disp(' ');