run preamble;

% caseFolder = '/home/lunet/ttcjc/OpenFOAM/ttcjc-7/results/Windsor_Upstream_2023/Windsor_SB_wW_Upstream_SC';
caseFolder = '/home/lunet/ttcjc/OpenFOAM/ttcjc-7/results/Windsor_fullScale/Windsor_SB_fullScale_multiPhase_20deg';


%%

tic;

evalc('parpool(''threads'');');

%%%%

% Identify Distributed Directories
dataDirs = dir([caseFolder, '/LagrangianExtractionPlane']);

i = 1;
while i <= height(dataDirs)

    if isnan(str2double(dataDirs(i).name))
        dataDirs(i,:) = [];
    else
        i = i + 1;
    end

end
clear i;

dataFiles = dir([caseFolder, '/LagrangianExtractionPlane/', dataDirs(1).name, '/LagrangianExtractionPlaneData_*']);

for i = 1:height(dataFiles)

    % Concatenate Distributed Files
    if i ~= 1
        disp(' ');
    end

    disp(['    Concatenating ', dataFiles(i).name, ' Files...']);
    
    contentFloat = [];

    for j = 1:height(dataDirs)
        fileID = fopen([caseFolder, '/LagrangianExtractionPlane/', dataDirs(j).name, '/', dataFiles(i).name]);
        contentRaw = textscan(fileID, '%f32 %u32 %u32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32', 'headerLines', 0, 'delimiter', '\n');

        contentFloat = [contentFloat; cell2mat(contentRaw([1,(4:end)]))]; %#ok<AGROW>

        fclose(fileID);
    end
    clear j;

    % Identify Plane Position
    if contains(dataFiles(i).name, '-')
        plane = ['X_N', dataFiles(i).name((strfind(dataFiles(i).name, '_') + 2):(strfind(dataFiles(i).name, '.') - 1))...
                 '_', dataFiles(i).name((strfind(dataFiles(i).name, '.') + 1):end)];
    else
        plane = ['X_P', dataFiles(i).name((strfind(dataFiles(i).name, '_') + 1):(strfind(dataFiles(i).name, '.') - 1)), ...
                 '_', dataFiles(i).name((strfind(dataFiles(i).name, '.') + 1):end)];
    end

    LagData.(plane).time = unique(contentFloat(:,1));
    
    nTimes = height(LagData.(plane).time);

    % Collate Particle Data
    disp(['        Collating ''', plane, ''' Data']);

    % Initialise Progress Bar
    wB = waitbar(0, ['Collating ''', plane, ''' Data'], 'name', 'Progress');
    wB.Children.Title.Interpreter = 'none';
    dQ = parallel.pool.DataQueue;
    afterEach(dQ, @parforWaitBar);
    parforWaitBar(wB, nTimes);

    % Perform Collation
    d = cell(nTimes, 1);
    nParticle = d;
    Up = d;
    
    uniqueTimes = LagData.(plane).time;
    content_time = contentFloat(:,1);
    content_d = contentFloat(:,2);
    content_nParticle = contentFloat(:,3);
    content_Up = contentFloat(:,[7,8,9]);
    parfor j = 1:nTimes
        index = content_time == uniqueTimes(j);
        
        d{j} = content_d(index); %#ok<PFBNS>
        nParticle{j} = content_nParticle(index); %#ok<PFBNS>
        Up{j} = content_Up(index,:); %#ok<PFBNS>
        
        % Update Waitbar
        send(dQ, []);
    end
    clear uniqueTimes content_time content_d content_nParticle content_Up;

    delete(wB);
    
    LagData.(plane).d = d; clear d;
    LagData.(plane).nParticle = nParticle; clear nParticle;
    LagData.(plane).Up = Up; clear Up;
end
clear i;

save('~/MATLAB/Testing/Dispersed Phase Processing/LagDataPlaneRawFS_20deg.mat', 'LagData', '-v7.3', '-noCompression');

%%%%

evalc('delete(gcp(''nocreate''));');

executionTime = toc;

disp(' ');

disp(['Run Time: ', num2str(executionTime), 's']);
