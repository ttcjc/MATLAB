%% POD Temporal Analysis Tool v1.0

clear variables;
close all;
clc;
evalc('delete(gcp(''nocreate''));');

nModes = 8; % Analyse First ‘nModes’ Energetic Modes [Even Integer]

nProc = maxNumCompThreads - 2; % Number of Processors Used for Parallel Collation

fig = 0; % Initialise Figure Tracking
figHold = 0; % Enable Overwriting of Figures

disp('===============================');
disp('POD Temporal Analysis Tool v1.0');
disp('===============================');

disp(' ');
disp(' ');


%% Changelog

% v1.0 - Initial Commit


%% Select Analysis Type

disp('Analysis Type');
disp('--------------');

disp(' ');

disp('Possible Analysis Types:');
disp('    A: Single-Field Fourier Transform');
disp('    B: Dual-Field Temporal Comparison');

valid = false;
while ~valid
    disp(' ');
    selection = input('Select Analysis Type [A/B]: ', 's');

    if selection == 'a' | selection == 'A' %#ok<OR2>
        format = 'A';
        valid = true;
    elseif selection == 'b' | selection == 'B' %#ok<OR2>
        format = 'B';
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;

disp(' ');
disp(' ');


%% Load POD Data

disp('POD Data Acquisition');
disp('---------------------');

valid = false;
while ~valid
    disp(' ');
    [fileName, filePath] = uigetfile('/mnt/Processing/Data/Numerical/MATLAB/*.mat', ...
                                      'Select POD Data');

    if contains(filePath, 'POD')
        disp(['Loading ''', fileName, '''...']);
        temporalData.A.dataID = load([filePath, fileName], 'dataID').dataID;
        temporalData.A.PODdata = load([filePath, fileName], 'PODdata').PODdata;
        disp('    Success');
        
        valid = true;
    else
        disp('WARNING: Invalid File Selection');
        clear fileName filePath;
    end
    
end
clear valid;

if contains(filePath, 'planarContaminantPOD')
    namePos = strfind(filePath, '/');
    fieldName = filePath((namePos(end - 1) + 1):(namePos(end) - 1));
elseif contains(filePath, 'planarPressurePOD')
    fieldName = 'pressure';
elseif contains(filePath, 'planarVelocityPOD')
    fieldName = 'velocity';
else
    error('Unknown POD Format');
end

temporalData = renameStructField(temporalData, 'A', fieldName);

switch format
    
    case 'B'
        
        valid = false;
        while ~valid
            disp(' ');
            [fileName, filePath] = uigetfile('/mnt/Processing/Data/Numerical/MATLAB/*.mat', ...
                                             'Select POD Data');
                                         
            if contains(filePath, 'POD')
                disp(['Loading ''', fileName, '''...']);
                temporalData.B.dataID = load([filePath, fileName], 'dataID').dataID;
                temporalData.B.PODdata = load([filePath, fileName], 'PODdata').PODdata;
                disp('    Success');
                
                valid = true;
            else
                disp('WARNING: Invalid File Selection');
                clear fileName filePath;
            end
            
        end
        clear valid;
        
        if contains(filePath, 'planarContaminantPOD')
            namePos = strfind(filePath, '/');
            fieldName = filePath((namePos(end - 1) + 1):(namePos(end) - 1));
        elseif contains(filePath, 'planarPressurePOD')
            fieldName = 'pressure';
        elseif contains(filePath, 'planarVelocityPOD')
            fieldName = 'velocity';
        else
            error('Unknown POD Format');
        end
        
        PODfields = fieldnames(temporalData);
        
        if strcmp(fieldName, PODfields{1})
            temporalData = renameStructField(temporalData, fieldName, [fieldName, 'A']);
            temporalData = renameStructField(temporalData, 'B', [fieldName, 'B']);
        else
            temporalData = renameStructField(temporalData, 'B', fieldName);
        end

end

PODfields = fieldnames(temporalData);

if mod(nModes, 2)
    nModes = nModes + 1;
end

disp(' ');
disp(' ');


%% Perform Fourier Transfom

for i = 1:height(PODfields)
    temporalData.(PODfields{i}).PODdata.time = temporalData.(PODfields{i}).PODdata.time - ...
                                            min(temporalData.(PODfields{i}).PODdata.time);
    n = pow2(nextpow2(height(temporalData.(PODfields{i}).PODdata.time)));
    Fs = round(1 / (temporalData.(PODfields{i}).PODdata.time(2) - temporalData.(PODfields{i}).PODdata.time(1)));
    
    % Figure Setup
    fig = fig + 1;
    figName = [PODfields{i}, '_POD_Mode_Temporal_Analysis'];
    set(figure(fig), 'outerPosition', [(25 + 875 * (fig - 1)), 25, 850, 850], 'name', figName);
    set(gca, 'lineWidth', 2, 'fontName', 'LM Mono 12', ...
             'fontSize', 20, 'layer', 'top');
    tiledlayout((nModes / 2),2);
    
    for j = 1:nModes
        x = zeros(n,1);
        x(1:height(temporalData.(PODfields{i}).PODdata.A_coeff(:,j))) = temporalData.(PODfields{i}).PODdata.A_coeff(:,j);
        [PSD, freq] = pwelch(x, 128, 64, n, Fs);
        Sr = (freq * 0.289) / 40;
        
        % Plot
        nexttile;
        hold on;
        plot(Sr(1:floor(n / 2)), rescale(PSD(1:floor(n / 2))), ...
             'lineWidth', 1.5, 'color', ([74, 24, 99] / 255));
        
        % Figure Formatting
        set(gca, 'xScale', 'log', 'yScale', 'log');
        title(['M', num2str(j)]);
        axis on;
        box on;
        grid off;
        xlim([1e-3; 1]);
        ylim([0; 10]);
        tickData = [1e-3, 1e-2, 1e-1, 1];
        xticks(tickData);
        tickData = [];
        yticks(tickData);
        xlabel('{\bf{Sr}}_{{\it{H}}}', 'fontName', 'LM Roman 12');
        ylabel('{\bf{PSD}}', 'fontName', 'LM Roman 12');
        hold off;
    
    end
    
end


%% Perform Mode Comparison

r = nan(nModes,1);
rMode = r;

for i = 1:nModes
    A1 = rescale(temporalData.(PODfields{1}).PODdata.A_coeff(:,i), -1, 1);
    
    for j = 1:nModes
        A2 = rescale(temporalData.(PODfields{2}).PODdata.A_coeff(:,j), -1, 1);
        rTemp = corr(A1, A2);
        
        if max(abs(rTemp), abs(r(i))) == abs(rTemp)
            r(i) = rTemp;
            rMode(i) = j;
        end
        
    end
    clear rTemp;
    
end

% Figure Setup
fig = fig + 1;
figName = 'POD_Mode_Temporal_Correlation';
set(figure(fig), 'outerPosition', [1945, 25, 850, 850], 'name', figName);
set(gca, 'lineWidth', 2, 'fontName', 'LM Mono 12', ...
         'fontSize', 20, 'layer', 'top');
tiledlayout((nModes / 2),2);

for i = 1:nModes
    A1 = rescale(temporalData.(PODfields{1}).PODdata.A_coeff(:,i), -1, 1);
    A2 = rescale(temporalData.(PODfields{2}).PODdata.A_coeff(:,rMode(i)), -1, 1);
    
    % Plot
    nexttile;
    hold on;
    scatter(A1, A2, 10, ([74, 24, 99] / 255), 'filled');

    % Figure Formatting
    title(['r = ', num2str(r(i))]);
    axis on;
    box on;
    grid off;
    xlim([-1.2; 1.2]);
    ylim([-1.2; 1.2]);
    tickData = [];
    xticks(tickData);
    tickData = [];
    yticks(tickData);
    xlabel(['{\bf{', PODfields{1}, '}}_{{M', num2str(i), '}}'], 'fontName', 'LM Roman 12');
    ylabel(['{\bf{', PODfields{2}, '}}_{{M', num2str(rMode(i)), '}}'], 'fontName', 'LM Roman 12');
    hold off;
end