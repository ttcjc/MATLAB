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


%% Case Initialisation

% caseFolder = uigetdir('~/OpenFOAM', 'Select Case');
caseFolder = ('~/Mount/Uni/OpenFOAM/ttcjc-7/results/Windsor_Square_wW_SC_Study');

disp(['Case: ', caseFolder]);
disp(' ');
disp(' ');

% Confirm Valid Case
if exist([caseFolder, '/0'], 'dir')
	timeDirs = dir(caseFolder);
else
	error('Invalid Case Directory');
end

% Identify Sample Data
if exist([caseFolder, '/postProcessing/forceCoeffs'], 'dir')
	forceCoeffs = true;
else
	forceCoeffs = false;
end

if exist([caseFolder, '/postProcessing/probesWake'], 'dir')
	velocityProbes = true;
else
	velocityProbes = false;
end

if ~forceCoeffs && ~velocityProbes
	error('No Suitable Sample Data Found for Target Case');
end
	
	
%% Force Coefficients
	
if forceCoeffs
	
	% Read Data
	coeffDirs = dir([caseFolder, '/postProcessing/forceCoeffs']);
		
	i = 1;
	while i <= size(coeffDirs,1)

		if isnan(str2double(coeffDirs(i,1).name))
			coeffDirs(i,:) = [];
		else
		i = i + 1;
		end

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
	
	% FFT
	if size(coeffData.directory,1) == 1
		FFT.forceCoeffs.time = coeffData.time{1,1};
		FFT.forceCoeffs.Cd = coeffData.Cd{1,1};
		FFT.forceCoeffs.Cl = coeffData.Cl{1,1};
	else
		
		for i = 1:size(coeffData.directory,1)-1
			FFT.forceCoeffs.time = vertcat(coeffData.time{i,1}, coeffData.time{i+1,1});
			FFT.forceCoeffs.Cd = vertcat(coeffData.Cd{i,1}, coeffData.Cd{i+1,1});
			FFT.forceCoeffs.Cl = vertcat(coeffData.Cl{i,1}, coeffData.Cl{i+1,1});
		end
		
	end
	
	FFT.forceCoeffs.time = FFT.forceCoeffs.time(end/3:end);
	FFT.forceCoeffs.Cd = FFT.forceCoeffs.Cd(end/3:end);
	FFT.forceCoeffs.Cl = FFT.forceCoeffs.Cl(end/3:end);
	
	FFT.forceCoeffs.Cd_Mean = mean(FFT.forceCoeffs.Cd);
	FFT.forceCoeffs.Cl_Mean = mean(FFT.forceCoeffs.Cl);
	FFT.forceCoeffs.Cd_Prime = FFT.forceCoeffs.Cd - FFT.forceCoeffs.Cd_Mean;
	FFT.forceCoeffs.Cl_Prime = FFT.forceCoeffs.Cl - FFT.forceCoeffs.Cl_Mean;
		
	FFT.forceCoeffs.time = FFT.forceCoeffs.time - min(FFT.forceCoeffs.time);
	FFT.forceCoeffs.freq_Samping = 1 / 2e-5;
	FFT.forceCoeffs.freq = linspace(0, 1, fix(size(FFT.forceCoeffs.time,1) / 2) + 1) * (FFT.forceCoeffs.freq_Samping / 2);
	
	FFT.forceCoeffs.Cd_FFT = fft(FFT.forceCoeffs.Cd_Prime) / size(FFT.forceCoeffs.time,1);
	FFT.forceCoeffs.Cl_FFT = fft(FFT.forceCoeffs.Cl_Prime) / size(FFT.forceCoeffs.time,1);
	
	
	% Figure Setup
	fig = fig + 1;
	figure('name', 'Force Coefficients');
	set(figure(fig), 'outerPosition', [25, 25, 800, 800]);
	
	% Plot in Time Domain
	subplot(2,1,1);
	hold on;
	
	plot(FFT.forceCoeffs.time, FFT.forceCoeffs.Cd, 'color', [0.21176, 0.06667, 0.38824]);
% 	plot(FFT.forceCoeffs.time, FFT.forceCoeffs.Cl, 'color', [0.71765, 0.00000, 0.38431]);
	
	% Figure Formatting
	xlim([min(FFT.forceCoeffs.time), max(FFT.forceCoeffs.time)])
	xlabel({' ', 'Time (\it{s})'});
	ylabel({'Force Coefficient', ' '});
	grid on;
	box on;
	hold off;
	
	% Plot in Frequency Domain
	subplot(2,1,2);
	hold on;
	
	index = 1:length(FFT.forceCoeffs.freq);
	plot(FFT.forceCoeffs.freq, abs(FFT.forceCoeffs.Cd_FFT(index)), 'color', [0.21176, 0.06667, 0.38824]);
% 	plot(FFT.forceCoeffs.freq, abs(FFT.forceCoeffs.Cl_FFT(index)), 'color', [0.71765, 0.00000, 0.38431]);
	
	% Figure Formatting
	xlabel({' ', 'Frequency (\it{Hz})'});
	ylabel({'Magnitude', ' '});
	grid on;
	box on;
	hold off;
	
end


%% Velocity Probes

if velocityProbes
	% Read Data
	probeDirs = dir([caseFolder, '/postProcessing/probesWake']);
	
	i = 1;
	while i <= size(probeDirs,1)

		if isnan(str2double(probeDirs(i,1).name))
			probeDirs(i,:) = [];
		else
			i = i + 1;
		end

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
	
	% FFT
	if size(probeData.directory,1) == 1
		FFT.velocityProbes.time = probeData.time{1,1};
		FFT.velocityProbes.probeA = probeData.probeA{1,1}(:,4);
		FFT.velocityProbes.probeB = probeData.probeB{1,1}(:,4);
		FFT.velocityProbes.probeC = probeData.probeC{1,1}(:,4);
		FFT.velocityProbes.probeD = probeData.probeD{1,1}(:,4);
		FFT.velocityProbes.probeE = probeData.probeE{1,1}(:,4);
	else
		
		for i = 1:size(probeData.directory,1)-1
			FFT.velocityProbes.time = vertcat(probeData.time{i,1}, probeData.time{i+1,1});
			FFT.velocityProbes.probeA = vertcat(probeData.probeA{i,1}, probeData.probeA{i+1,1});
			FFT.velocityProbes.probeB = vertcat(probeData.probeB{i,1}, probeData.probeB{i+1,1});
			FFT.velocityProbes.probeC = vertcat(probeData.probeC{i,1}, probeData.probeC{i+1,1});
			FFT.velocityProbes.probeD = vertcat(probeData.probeD{i,1}, probeData.probeD{i+1,1});
			FFT.velocityProbes.probeE = vertcat(probeData.probeE{i,1}, probeData.probeE{i+1,1});
		end
		
	end
	
	FFT.velocityProbes.probeA_Mean = mean(FFT.velocityProbes.probeA);
	FFT.velocityProbes.probeB_Mean = mean(FFT.velocityProbes.probeB);
	FFT.velocityProbes.probeC_Mean = mean(FFT.velocityProbes.probeC);
	FFT.velocityProbes.probeD_Mean = mean(FFT.velocityProbes.probeD);
	FFT.velocityProbes.probeE_Mean = mean(FFT.velocityProbes.probeE);
	FFT.velocityProbes.probeA_Prime = FFT.velocityProbes.probeA - FFT.velocityProbes.probeA_Mean;
	FFT.velocityProbes.probeB_Prime = FFT.velocityProbes.probeB - FFT.velocityProbes.probeB_Mean;
	FFT.velocityProbes.probeC_Prime = FFT.velocityProbes.probeC - FFT.velocityProbes.probeC_Mean;
	FFT.velocityProbes.probeD_Prime = FFT.velocityProbes.probeD - FFT.velocityProbes.probeD_Mean;
	FFT.velocityProbes.probeE_Prime = FFT.velocityProbes.probeE - FFT.velocityProbes.probeE_Mean;
	
	FFT.velocityProbes.time = FFT.velocityProbes.time - min(FFT.velocityProbes.time);
	FFT.velocityProbes.freq_Samping = 1 / 2e-5;
	FFT.velocityProbes.freq = linspace(0, 1, fix(size(FFT.velocityProbes.time,1) / 2) + 1) * (FFT.velocityProbes.freq_Samping / 2);
	
	FFT.velocityProbes.probeA_FFT = fft(FFT.velocityProbes.probeA_Prime) / size(FFT.velocityProbes.time,1);
	FFT.velocityProbes.probeB_FFT = fft(FFT.velocityProbes.probeB_Prime) / size(FFT.velocityProbes.time,1);
	FFT.velocityProbes.probeC_FFT = fft(FFT.velocityProbes.probeC_Prime) / size(FFT.velocityProbes.time,1);
	FFT.velocityProbes.probeD_FFT = fft(FFT.velocityProbes.probeD_Prime) / size(FFT.velocityProbes.time,1);
	FFT.velocityProbes.probeE_FFT = fft(FFT.velocityProbes.probeE_Prime) / size(FFT.velocityProbes.time,1);
	
	% Figure Setup
	fig = fig + 1;
	figure('name', 'Velocity Probes');
	set(figure(fig), 'outerPosition', [25, 25, 800, 800]);
	
	% Plot in Time Domain
	subplot(2,1,1);
	hold on;
	
	plot(FFT.velocityProbes.time, FFT.velocityProbes.probeA);
% 	plot(FFT.velocityProbes.time, FFT.velocityProbes.probeB);
% 	plot(FFT.velocityProbes.time, FFT.velocityProbes.probeC);
% 	plot(FFT.velocityProbes.time, FFT.velocityProbes.probeD);
% 	plot(FFT.velocityProbes.time, FFT.velocityProbes.probeE);
	
	% Figure Formatting
	xlim([min(FFT.velocityProbes.time), max(FFT.velocityProbes.time)])
	xlabel({' ', 'Time (\it{s})'});
	ylabel({'|U| (\it{m/s})', ' '});
	grid on;
	box on;
	hold off;
	
	% Plot in Frequency Domain
	subplot(2,1,2);
	hold on;
	
	index = 1:length(FFT.velocityProbes.freq);
	plot(FFT.velocityProbes.freq, abs(FFT.velocityProbes.probeA_FFT(index)));
% 	plot(FFT.velocityProbes.freq, abs(FFT.velocityProbes.probeB_FFT(index)));
% 	plot(FFT.velocityProbes.freq, abs(FFT.velocityProbes.probeC_FFT(index)));
% 	plot(FFT.velocityProbes.freq, abs(FFT.velocityProbes.probeD_FFT(index)));
% 	plot(FFT.velocityProbes.freq, abs(FFT.velocityProbes.probeE_FFT(index)));
	
	% Figure Formatting
	xlabel({' ', 'Frequency (\it{Hz})'});
	ylabel({'Magnitude', ' '});
	grid on;
	box on;
	hold off;	
end


%% Cleaning

clearvars -except FFT;
disp(' ');