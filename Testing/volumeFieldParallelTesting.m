clc;
evalc('delete(gcp(''nocreate''));');


%% Collate Particles of Interest

evalc('parpool(''threads'');');

tic;
index = cell(nTimes,1);

d = LagData.d;
positionCartesian = LagData.positionCartesian;
parfor i = 1:nTimes
    
    if ~isempty(positionCartesian{i})
        index{i} = find(((d{i} * 1e6) >= dLims(1)) & ...
                        ((d{i} * 1e6) <= dLims(2)) & ...
                        (positionCartesian{i}(:,1) >= xLimsData(1)) & ...
                        (positionCartesian{i}(:,1) <= xLimsData(2)) & ...
                        (positionCartesian{i}(:,2) >= yLimsData(1)) & ...
                        (positionCartesian{i}(:,2) <= yLimsData(2)) & ...
                        (positionCartesian{i}(:,3) >= zLimsData(1)) & ...
                        (positionCartesian{i}(:,3) <= zLimsData(2)));
    end

end
clear i d positionCartesian;
toc;

evalc('delete(gcp(''nocreate''));');


disp(' ');
disp(' ');


tic;
index = cell(nTimes,1);

for i = 1:nTimes
    
    if ~isempty(LagData.positionCartesian{i})
        index{i} = find(((LagData.d{i} * 1e6) >= dLims(1)) & ...
                        ((LagData.d{i} * 1e6) <= dLims(2)) & ...
                        (LagData.positionCartesian{i}(:,1) >= xLimsData(1)) & ...
                        (LagData.positionCartesian{i}(:,1) <= xLimsData(2)) & ...
                        (LagData.positionCartesian{i}(:,2) >= yLimsData(1)) & ...
                        (LagData.positionCartesian{i}(:,2) <= yLimsData(2)) & ...
                        (LagData.positionCartesian{i}(:,3) >= zLimsData(1)) & ...
                        (LagData.positionCartesian{i}(:,3) <= zLimsData(2)));
    end

end
clear i;
toc;


%% Assign Particles to Volume Nodes

evalc('parpool(''threads'');');

tic;
totalParticles = cellfun(@height, LagData.positionCartesian);
index = cell(nTimes,1); % Array Position of Closest Mesh Node

positionCartesian = LagData.positionCartesian;
parfor i = 1:nTimes
    
    if totalParticles(i) > 0
        index{i} = zeros(height(positionCartesian{i}),3);
        
        for j = 1:totalParticles(i)
            [~, index{i}(j,1)] = min(abs(positionCartesian{i}(j,1) - x(:,1,1))); %#ok<PFBNS>
            [~, index{i}(j,2)] = min(abs(positionCartesian{i}(j,2) - y(1,:,1))); %#ok<PFBNS>
            [~, index{i}(j,3)] = min(abs(positionCartesian{i}(j,3) - z(1,1,:))); %#ok<PFBNS>

        end
    
    end

end
clear positionCartesian;
toc;

evalc('delete(gcp(''nocreate''));');


disp(' ');
disp(' ');


tic;
totalParticles = cellfun(@height, LagData.positionCartesian);
index = cell(nTimes,1); % Array Position of Closest Mesh Node

for i = 1:nTimes
    
    if totalParticles(i) > 0
        index{i} = zeros(height(LagData.positionCartesian{i}),3);
        
        for j = 1:totalParticles(i)
            [~, index{i}(j,1)] = min(abs(LagData.positionCartesian{i}(j,1) - x(:,1,1)));
            [~, index{i}(j,2)] = min(abs(LagData.positionCartesian{i}(j,2) - y(1,:,1)));
            [~, index{i}(j,3)] = min(abs(LagData.positionCartesian{i}(j,3) - z(1,1,:)));

        end
    
    end

end
toc;