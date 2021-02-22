%% Aerodynamic Coefficient Plotter v1.0

clearvars;
close all;
clc;

fig = 0; %#ok<*NASGU>
figHold = 0; %#ok<*NASGU>

disp ('=======================================');
disp ('%% Aerodynamic Coefficient Plotter v1.0');
disp ('=======================================');
disp (' ');


%% Changelog

% v1.0 - Initial Commit


%% Case Selection

disp('CASE SELECTION');
disp('--------------');
disp(' ');

disp('Select Case Folder:');
% caseFolder = uigetdir('~/OpenFOAM');
caseFolder = ('~/Mount/Uni/OpenFOAM/ttcjc-7/results/Windsor_Square_wW_SC_Study');

disp(['Case: ', caseFolder]);
disp(' ');


%% Data Acquisition

if exist([caseFolder, '/postProcessing/forceCoeffs'], 'dir')
	coeffDirs = dir([caseFolder, '/postProcessing/forceCoeffs']);
	
	i = 1;
	while i <= size(coeffDirs,1)

		if isnan(str2double(coeffDirs(i,1).name))
			coeffDirs(i,:) = [];
		else
			i = i + 1;
		end

	end
	
	if isempty(coeffDirs)
		error('No Coefficient Data Found for Target Case');
	end
	
else
	error('No Coefficient Data Found for Target Case');
end

for i = 1:size(coeffDirs,1)
	coeffData.directory{i,1} = coeffDirs(i,1).name;
	
	fileID = fopen([caseFolder, '/postProcessing/forceCoeffs/', coeffData.directory{i,1}, '/forceCoeffs.dat']);
	content = textscan(fileID, '%f %f %f %f %f %f', 'headerLines', 9, 'delimiter', ' ', 'multipleDelimsAsOne', 1);
	
	coeffData.time{i,1} = content{1,1};
% 	coeffData.Cm{i,1} = content{1,2};
	coeffData.Cd{i,1} = content{1,3};
	coeffData.Cl{i,1} = content{1,4};
% 	coeffData.Cl_f{i,1} = content{1,5};
% 	coeffData.Cl_r{i,1} = content{1,6};
	
	fclose(fileID);
end

coeffData.Cd_Exp = 0.3889;
coeffData.Cl_Exp = 0.1435;


%% Blockage Correction

% Experimental
At = ((2 * (0.96 + (1.695 * tan(atan(0.01 / 3.6))))) * 1.32) - (4 * 0.01125); 
Am = (0.289 * 0.389) + (2 * (0.05 * 0.055));
U = 39.56;
rho = 1.269;

Ucorr = (U * At) / (At - Am);

Fd = coeffData.Cd_Exp * (0.5 * rho * U^2 * Am);
Fl = coeffData.Cl_Exp * (0.5 * rho * U^2 * Am);

coeffData.Cd_Exp_Corr = Fd / (0.5 * rho * Ucorr^2 * Am);
coeffData.Cl_Exp_Corr = Fl / (0.5 * rho * Ucorr^2 * Am);

% Numerical
At = (2 * (0.96 + (4.704 * tan(atan(0.02613 / 9.408)))) * 1.32);
Am = (0.289 * 0.389) + (2 * (0.046 * 0.055));
U = 40;
rho = 1.269;

Ucorr = (U * At) / (At - Am);

for i = 1:size(coeffData.directory,1)
	Fd = coeffData.Cd{i,1} * (0.5 * rho * U^2 * Am);
	Fl = coeffData.Cl{i,1} * (0.5 * rho * U^2 * Am);

	coeffData.Cd_Corr{i,1} = Fd / (0.5 * rho * Ucorr^2 * Am);
	coeffData.Cl_Corr{i,1} = Fl / (0.5 * rho * Ucorr^2 * Am);
end


%% Plot

% Figure Setup
fig = fig + 1;
figure('name', 'Force Coefficients');
hold on;
set(figure(fig), 'outerPosition', [25, 25, 800, 800]);

% Plot
for i = 1:size(coeffData.directory,1)
	plot(coeffData.time{i,1}, coeffData.Cd_Corr{i,1}, 'color', [0.21176, 0.06667, 0.38824]);
	plot(coeffData.time{i,1}, coeffData.Cl_Corr{i,1}, 'color', [0.71765, 0.00000, 0.38431]);
end

% Figure Formatting
xlim([0, 3])
ylim([-0.8, 0.8]);
xticks(0.5:0.5:2.5);
yticks(-0.6:0.2:0.6);
xtickformat('%.2f');
ytickformat('%.2f');
xlabel({' ', 'Time (\it{s})'});
ylabel({'Force Coefficient', ' '});
grid on;
box on;
legend('C_D', 'C_L', 'location', 'northOutside', 'orientation', 'vertical');
legend boxoff;
set(gca, 'units', 'normalized', 'position', [0.1275, 0.1275, 0.745, 0.745], ...
	     'fontName', 'LM Roman 12', 'fontSize', 12, 'layer', 'top');
hold off;

namePos = max(strfind(caseFolder, '/'));
savefig(fig, ['~/MATLAB/Output/Figures/', caseFolder(namePos(end):end), '_Force_Coeffs_']);
print(fig, ['~/MATLAB/Output/Figures/', caseFolder(namePos(end):end), '_Force_Coeffs_'], '-dpng', '-r300');


%% Time-Averaged Values

Cd_Corr = coeffData.Cd_Corr{1,1};
Cl_Corr = coeffData.Cl_Corr{1,1};

if size(coeffData.directory,1) > 1
	
	for i = 2:size(coeffData.directory,1)
		Cd_Corr = vertcat(Cd_Corr, coeffData.Cd_Corr{i,1}); %#ok<AGROW>
		Cl_Corr = vertcat(Cl_Corr, coeffData.Cl_Corr{i,1}); %#ok<AGROW>
	end
	
end

coeffData.Cd_Corr = Cd_Corr;
coeffData.Cl_Corr = Cl_Corr;

coeffData.Cd_Corr_Mean = round(mean(coeffData.Cd_Corr((end / 3):end)),4);
coeffData.Cl_Corr_Mean = round(mean(coeffData.Cl_Corr((end / 3):end)),4);

coeffData.Cd_Error = ((coeffData.Cd_Exp_Corr - coeffData.Cd_Corr_Mean) / coeffData.Cd_Exp_Corr) * 100;
coeffData.Cl_Error = ((coeffData.Cl_Exp_Corr - coeffData.Cl_Corr_Mean) / coeffData.Cl_Exp_Corr) * 100;

disp(['Time-Averaged Numerical Drag Coefficient = ', num2str(coeffData.Cd_Corr_Mean)]);
disp(['Time-Averaged Numerical Lift Coefficient = ', num2str(coeffData.Cl_Corr_Mean)]);
disp(' ');
disp(['Experimental Drag Coefficient = ', num2str(coeffData.Cd_Exp_Corr)]);
disp(['Experimental Lift Coefficient = ', num2str(coeffData.Cl_Exp_Corr)]);
disp(' ');
disp(['Drag Coefficient Error = ', num2str(coeffData.Cd_Error), '%']);
disp(['Lift Coefficient Error = ', num2str(coeffData.Cl_Error), '%']);


%% Cleaning

clearvars -except coeffData;
disp(' ');