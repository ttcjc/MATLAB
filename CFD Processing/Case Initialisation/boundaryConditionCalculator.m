%% Turbulence Boundary Condition Calculator v1.1

clearvars;
close all;
clc;

fig = 0; %#ok<*NASGU>
figHold = 0; %#ok<*NASGU>

disp ('=============================================');
disp ('Turbulence Boundary Condition Calculator v1.1');
disp ('=============================================');
disp (' ');


%% Changelog

% v1.0 - Initial Commit
% v1.1 - Added Basic Error Catching


%% Data Input

variables = {
             'Characteristic Length [m]'
             'Turbulence Intensity [%]'
             'Freestream Velocity [m/s]'
            };

valid = false;
while ~valid
    userInput = inputdlg(variables, 'Input Turbulence Parameters', ...
                         [1, 50], {'1.044', '0.2', '40'});
    
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

L = values(1);
I = values(2) / 100;
U = values(3);

disp(['Characteristic Length:         ' num2str(L) ' [m]']);
disp(['Turbulence Intensity:          ' num2str(I * 100) ' [%]']);
disp(['Freestream Velocity:           ' num2str(U) ' [m/s]']);
disp(' ')


%% Boundary Condition Calculation

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