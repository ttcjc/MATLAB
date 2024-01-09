run preamble;

nCellsA = 1000;
nValuesA = 2e6;

nCellsB = 900;
nValuesB = 1.8e6;

nCellsC = 1100;
nValuesC = 2.2e6;


%% Make Some Big Cell Arrays and Extract a Random Value From One

testVarA = cell(nCellsA,1);

for i = 1:nCellsA
    testVarA{i} = rand([nValuesA,1]);
end
clear i;

testVarB = cell(nCellsB,1);

for i = 1:nCellsB
    testVarB{i} = rand([nValuesB,1]);
end
clear i;

testVarC = cell(nCellsC,1);

for i = 1:nCellsC
    testVarC{i} = rand([nValuesC,1]);
end
clear i;

randVarB1 = randperm(nCellsB, 1); randVarB2 = randperm(nValuesB, 1);
testVarB = testVarB{randVarB1}(randVarB2);

randVarA1 = randperm(nCellsA, 1); randVarA2 = randperm(nValuesA, 1);
testVarA = testVarA{randVarA1}(randVarA2);

randVarC1 = randperm(nCellsC, 1); randVarC2 = randperm(nValuesC, 1);
testVarC = testVarC{randVarC1}(randVarC2);


%% Make Some Big Cell Arrays and Extract a Random Value From Each Using a Function Instead

[testVarD, testVarE, testVarF] = bigVars(nCellsA, nValuesA, nCellsB, nValuesB, nCellsC, nValuesC);


%% Local Functions

function [outA, outB, outC] = bigVars(inA1, inA2, inB1, inB2, inC1, inC2)
    
    outA = cell(inA1,1);
    
    for i = 1:inA1
        outA{i} = rand([inA2,1]);
    end
    clear i;
    
    outB = cell(inB1,1);
    
    for i = 1:inB1
        outB{i} = rand([inB2,1]);
    end
    clear i;
    
    outC = cell(inC1,1);
    
    for i = 1:inC1
        outC{i} = rand([inC2,1]);
    end
    clear i;
    
    randVarB1 = randperm(inB1, 1); randVarB2 = randperm(inB2, 1);
    outB = outB{randVarB1}(randVarB2);
    
    randVarA1 = randperm(inA1, 1); randVarA2 = randperm(inA2, 1);
    outA = outA{randVarA1}(randVarA2);
    
    randVarC1 = randperm(inC1, 1); randVarC2 = randperm(inC2, 1);
    outC = outC{randVarC1}(randVarC2);
    
end
