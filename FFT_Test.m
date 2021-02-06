%% Fast Fourier Transforms v1.0

clearvars;
close all;
clc;

fig = 0; %#ok<*NASGU>
figHold = 0; %#ok<*NASGU>

disp ('============================');
disp ('Fast Fourier Transforms v1.0');
disp ('============================');
disp (' ');


%% Changelog

% v1.0 - Initial Commit


%% Data Read

caseFolder = ('/home/cjc96/Mount/Athena/OpenFOAM/ttcjc-7/run/Windsor_Square_Sampling');

% Force Coefficients
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
	content = textscan(fileID, '%f %f %f %f %f %f', 'headerLines', 9, 'delimiter', ' ', 'multipleDelimsAsOne', 1, 'collectOutput', 1);
	
	coeffData.time{i,1} = content{1,1}(:,1);
	coeffData.Cd{i,1} = content{1,1}(:,3);
	coeffData.Cl{i,1} = content{1,1}(:,4);
	
	fclose(fileID);
end

% Velocity Probes
if exist([caseFolder, '/postProcessing/probesWake'], 'dir')
	probeDirs = dir([caseFolder, '/postProcessing/probesWake']);
	
	i = 1;
	while i <= size(probeDirs,1)

		if isnan(str2double(probeDirs(i,1).name))
			probeDirs(i,:) = [];
		else
			i = i + 1;
		end

	end
	
	if isempty(probeDirs)
		error('No Probe Data Found for Target Case');
	end
	
else
	error('No Probe Data Found for Target Case');
end

for i = 1:size(probeDirs,1)
	probeData.directory{i,1} = probeDirs(i,1).name;
	
	fileID = fopen([caseFolder, '/postProcessing/probesWake/', probeData.directory{i,1}, '/U']);
	content = textscan(fileID, '%f (%f %f %f) (%f %f %f) (%f %f %f) (%f %f %f) (%f %f %f)', 'headerLines', 7, 'delimiter', ' ', 'multipleDelimsAsOne', 1, 'collectOutput', 1);
	
	probeData.time{i,1} = content{1,1}(:,1);
	probeData.probeA{i,1} = content{1,1}(:,2:4);
	probeData.probeB{i,1} = content{1,1}(:,5:7);
	probeData.probeC{i,1} = content{1,1}(:,8:10);
	probeData.probeD{i,1} = content{1,1}(:,11:13);
	probeData.probeE{i,1} = content{1,1}(:,14:16);
	
	probeData.probeA{i,1}(:,4) = sqrt(probeData.probeA{i,1}(:,1).^2 + probeData.probeA{i,1}(:,2).^2 + probeData.probeA{i,1}(:,3).^2);
	probeData.probeB{i,1}(:,4) = sqrt(probeData.probeB{i,1}(:,1).^2 + probeData.probeB{i,1}(:,2).^2 + probeData.probeB{i,1}(:,3).^2);
	probeData.probeC{i,1}(:,4) = sqrt(probeData.probeC{i,1}(:,1).^2 + probeData.probeC{i,1}(:,2).^2 + probeData.probeC{i,1}(:,3).^2);
	probeData.probeD{i,1}(:,4) = sqrt(probeData.probeD{i,1}(:,1).^2 + probeData.probeD{i,1}(:,2).^2 + probeData.probeD{i,1}(:,3).^2);
	probeData.probeE{i,1}(:,4) = sqrt(probeData.probeE{i,1}(:,1).^2 + probeData.probeE{i,1}(:,2).^2 + probeData.probeE{i,1}(:,3).^2);
	
	fclose(fileID);
end


%% Plot Against Time

% Figure Setup
fig = fig + 1;
figure('name', 'Force Coefficients');
hold on;
set(figure(fig), 'outerPosition', [25, 25, 800, 800]);

% Plot
for i = 1:size(coeffData.directory,1)
	plot(coeffData.time{i,1}, coeffData.Cd{i,1});
	plot(coeffData.time{i,1}, coeffData.Cl{i,1});
end

% Figure Formatting
xlim([0.45, 1.05])
xlabel({' ', 'Time (\it{s})'});
ylabel({'Force Coefficient', ' '});
grid on;
box on;
set(gca, 'units', 'normalized', 'position', [0.1275, 0.1275, 0.745, 0.745], ...
	     'fontName', 'LM Roman 12', 'fontSize', 12, 'layer', 'top');
hold off;

% Figure Setup
fig = fig + 1;
figure('name', 'Velocity Probes');
hold on;
set(figure(fig), 'outerPosition', [25, 25, 800, 800]);

% Plot
for i = 1:size(probeData.directory,1)
% 	plot(probeData.time{i,1}, probeData.probeA{i,1}(:,4));
% 	plot(probeData.time{i,1}, probeData.probeB{i,1}(:,4));
% 	plot(probeData.time{i,1}, probeData.probeC{i,1}(:,4));
% 	plot(probeData.time{i,1}, probeData.probeD{i,1}(:,4));
% 	plot(probeData.time{i,1}, probeData.probeE{i,1}(:,4));
end

% Figure Formatting
xlim([0.45, 1.05])
xlabel({' ', 'Time (\it{s})'});
ylabel({'U (\it{m/s})', ' '});
grid on;
box on;
set(gca, 'units', 'normalized', 'position', [0.1275, 0.1275, 0.745, 0.745], ...
	     'fontName', 'LM Roman 12', 'fontSize', 12, 'layer', 'top');
hold off;


%% Fast Fourier Transform

Cd_Mean = mean(coeffData.Cd{1,1}(length(coeffData.Cd{1,1})/3:end));
Cd_Prime = coeffData.Cd{1,1} - Cd_Mean;
Cl_Mean = mean(coeffData.Cl{1,1}(length(coeffData.Cl{1,1})/3:end));
Cl_Prime = coeffData.Cl{1,1} - Cl_Mean;

t = coeffData.time{1,1};
Fs = 1/2e-5;

xDr = Cd_Prime(length(Cd_Prime)/3:end);
nDr = pow2(nextpow2(length(xDr)));
xLi = Cl_Prime(length(Cl_Prime)/3:end);;
nLi = pow2(nextpow2(length(xLi)));

yDr = fft(xDr,nDr);
fDr = (0:nDr-1)*Fs/nDr;
yLi = fft(xLi,nLi);
fLi = (0:nLi-1)*Fs/nLi;

powerDr = abs(yDr).^2/nDr;
powerLi = abs(yLi).^2/nLi;

% Figure Setup
fig = fig + 1;
figure('name', 'Force Coefficient FFT');
hold on;
set(figure(fig), 'outerPosition', [25, 25, 800, 800]);

% Plot
plot(fDr(1:floor(nDr/2)),powerDr(1:floor(nDr/2)));
plot(fLi(1:floor(nLi/2)),powerLi(1:floor(nLi/2)));

% Figure Formatting
xlabel({' ', 'Frequency (\it{Hz})'});
ylabel({'Power', ' '});
grid on;
box on;
set(gca, 'Xscale', 'log');
set(gca, 'units', 'normalized', 'position', [0.1275, 0.1275, 0.745, 0.745], ...
	     'fontName', 'LM Roman 12', 'fontSize', 12, 'layer', 'top');
hold off;

U_Mean_A = mean(probeData.probeA{i,1}(:,4));
U_Prime_A = probeData.probeA{i,1}(:,4) - U_Mean_A;
U_Mean_B = mean(probeData.probeB{i,1}(:,4));
U_Prime_B = probeData.probeB{i,1}(:,4) - U_Mean_B;
U_Mean_C = mean(probeData.probeC{i,1}(:,4));
U_Prime_C = probeData.probeC{i,1}(:,4) - U_Mean_C;
U_Mean_D = mean(probeData.probeD{i,1}(:,4));
U_Prime_D = probeData.probeD{i,1}(:,4) - U_Mean_D;
U_Mean_E = mean(probeData.probeE{i,1}(:,4));
U_Prime_E = probeData.probeE{i,1}(:,4) - U_Mean_E;

t = probeData.time{1,1};
Fs = 1/2e-5;

xA = U_Prime_A;
nA = pow2(nextpow2(length(xA)));
xB = U_Prime_B;
nB = pow2(nextpow2(length(xB)));
xC = U_Prime_C;
nC = pow2(nextpow2(length(xC)));
xD = U_Prime_D;
nD = pow2(nextpow2(length(xD)));
xE = U_Prime_E;
nE = pow2(nextpow2(length(xE)));

yA = fft(xA,nA);
fA = (0:nA-1)*Fs/nA;
yB = fft(xB,nB);
fB = (0:nB-1)*Fs/nB;
yC = fft(xC,nC);
fC = (0:nC-1)*Fs/nC;
yD = fft(xD,nD);
fD = (0:nD-1)*Fs/nD;
yE = fft(xE,nE);
fE = (0:nE-1)*Fs/nE;

powerA = abs(yA).^2/nA;
powerB = abs(yB).^2/nB;
powerC = abs(yC).^2/nC;
powerD = abs(yD).^2/nD;
powerE = abs(yE).^2/nE;

% Figure Setup
fig = fig + 1;
figure('name', 'Velocity Probe FFT');
hold on;
set(figure(fig), 'outerPosition', [25, 25, 800, 800]);

% Plot
% plot(fA(1:floor(nA/2)),powerA(1:floor(nA/2)));
% plot(fB(1:floor(nB/2)),powerB(1:floor(nB/2)));
% plot(fC(1:floor(nC/2)),powerC(1:floor(nC/2)));
% plot(fD(1:floor(nD/2)),powerD(1:floor(nD/2)));
% plot(fE(1:floor(nE/2)),powerE(1:floor(nE/2)));

% Figure Formatting
xlabel({' ', 'Frequency (\it{Hz})'});
ylabel({'Power', ' '});
grid on;
box on;
set(gca, 'Xscale', 'log');
set(gca, 'units', 'normalized', 'position', [0.1275, 0.1275, 0.745, 0.745], ...
	     'fontName', 'LM Roman 12', 'fontSize', 12, 'layer', 'top');
hold off;