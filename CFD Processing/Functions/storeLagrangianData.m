%% Lagrangian Data Storage v1.0


%% Changelog

% v1.0 - Initial Commit


%% Main Function

function [] = storeLagrangianData(caseFolder)

	clc;

	% Confirm Support
	if ~contains(caseFolder, ["Lag_Test", "Test_Block", "Windsor_Square"])
		error('Unsupported Case');
	end
	
	% Identify Time Directories
	[timeDirs, ~] = timeDirectories(caseFolder);
	
	disp(' ');
	disp(' ');
	
	%  Store Lagrangian Data
	lagrangianData(caseFolder, timeDirs);

end