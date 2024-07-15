run preamble;

% caseFolder = '/home/lunet/ttcjc/OpenFOAM/ttcjc-7/results/Windsor_Upstream_2023/Windsor_SB_wW_Upstream_SC';
caseFolder = '/home/lunet/ttcjc/OpenFOAM/ttcjc-7/results/Windsor_fullScale/Windsor_SB_fullScale_multiPhase_coupled';


%%

tic;

evalc('parpool(''threads'');');

%%%%

% Identify Distributed Directories
dataDirs = dir([caseFolder, '/LagrangianSurfaceContamination']);

i = 1;
while i <= height(dataDirs)

    if isnan(str2double(dataDirs(i).name))
        dataDirs(i,:) = [];
    else
        i = i + 1;
    end

end
clear i;

% Concatenate Distributed Files
disp('    Concatenating ''LagrangianSurfaceContaminationData'' Files...');

contentInt = [];
contentFloat = [];

for i = 1:height(dataDirs)
    dataFile = dir([caseFolder, '/LagrangianSurfaceContamination/', dataDirs(i).name, '/LagrangianSurfaceContaminationData']);

    fileID = fopen([caseFolder, '/LagrangianSurfaceContamination/', dataDirs(i).name, '/', dataFile.name]);
    contentRaw = textscan(fileID, '%f32 %u32 %u32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32', 'headerLines', 0, 'delimiter', '\n');

    contentInt = [contentInt; cell2mat(contentRaw(2:3))]; %#ok<AGROW>
    contentFloat = [contentFloat; cell2mat(contentRaw([1,(4:end)]))]; %#ok<AGROW>

    fclose(fileID);
end
clear i;

LagData.time = unique(contentFloat(:,1));

nTimes = height(LagData.time);

% Collate Particle Data
disp('        Collating Surface Data...');

% Initialise Progress Bar
wB = waitbar(0, 'Collating Surface Data', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';
dQ = parallel.pool.DataQueue;
afterEach(dQ, @parforWaitBar);
parforWaitBar(wB, nTimes);

% % Identify Base Data
% index = (round(contentFloat(:,4), 3) == single(1.933) & contentFloat(:,6) > single(0.182));
% contentInt = contentInt(index,:);
% contentFloat = contentFloat(index,:);

% Identify Underbody Data
index = (round(contentFloat(:,6), 3) == single(0.182));
contentInt = contentInt(index,:);
contentFloat = contentFloat(index,:);

% Perform Collation
d = cell(nTimes, 1);
nParticle = d;
Up = d;

uniqueTimes = LagData.time;
content_time = contentFloat(:,1);
content_d = contentFloat(:,2);
content_nParticle = contentFloat(:,3);
content_pos = contentFloat(:,[7,8,9]);
content_Up = contentFloat(:,[7,8,9]);
parfor i = 1:nTimes
    index = (content_time == uniqueTimes(i));

    d{i} = content_d(index); %#ok<PFBNS>
    nParticle{i} = content_nParticle(index); %#ok<PFBNS>
    Up{i} = content_Up(index,:); %#ok<PFBNS>

    % Update Waitbar
    send(dQ, []);
end
clear uniqueTimes content_time content_d content_nParticle content_Up;

delete(wB);

LagData.d = d; clear d;
LagData.nParticle = nParticle; clear nParticle;
LagData.Up = Up; clear Up;

save('~/MATLAB/Testing/Dispersed Phase Processing/LagDataSoilingRaw_Underbody_FS_Coupled.mat', 'LagData', '-v7.3', '-noCompression');

%%%%

evalc('delete(gcp(''nocreate''));');

executionTime = toc;

disp(' ');

disp(['Run Time: ', num2str(executionTime), 's']);