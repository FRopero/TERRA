% This function solves the TSP for the UGV in the TERRA algorithm.

% INPUTS
% V3: list of vertices
% ugv_tsp: UGV's TSP Genetic Algoritm Parameters

% OUTPUT
% UGV_path: LoL

function [rte,minDist,UGV_path] = tsp_ga_ugv(V3, ugv_tsp)
% Load Instance
xy = V3';
last_depot = 1;
popSize     = ugv_tsp.popSize;
tournaments = ugv_tsp.tournaments;
mutOper     = ugv_tsp.mutOper;
mutRate     = ugv_tsp.mutRate;
crossOper   = ugv_tsp.crossOper;
eliteP      = ugv_tsp.eliteP;

% Set Fixed Seval and compute the remaining parameters
seval   = 30000;
numIter = seval/popSize;
elite   = popSize*eliteP / 100;

% Initialize default configuration of the genetic algorithm
dmat        = [];
minDist     = 0;
elitesMut   = mutRate * popSize;

if isempty(dmat)
    nPoints = size(xy,1);
    a = meshgrid(1:nPoints);
    dmat = reshape(sqrt(sum((xy(a,:)-xy(a',:)).^2,2)),nPoints,nPoints);
end

% Verify Inputs
[N,~] = size(xy);
[nr,nc] = size(dmat);
if N ~= nr || N ~= nc
    error('Invalid XY or DMAT inputs!')
end
n = N;

% Sanity Checks
popSize     = 4*ceil(popSize/4);
numIter     = max(1,round(real(numIter(1))));

% Initialize the Population
pop = zeros(popSize,n);
for k = 1:popSize
    pop(k,:) = randperm(n);
end

% Run the GA
globalMin = Inf;
totalDist = zeros(1,popSize);
distHistory = zeros(3,numIter);
tmpPop = zeros(1,n);
for iter = 1:numIter
    newPop = [];
    
    % Evaluate Each Population Member (Calculate Total Distance)
    for p = 1:popSize
        if (last_depot == 1)
            d = dmat(pop(p,n),pop(p,1)); % Closed Path
        else
            d = 0; % Open Path
        end
        for k = 2:n
            d = d + dmat(pop(p,k-1),pop(p,k));
        end
        totalDist(p) = d;
    end
    
    % Find the Best Route in the Population
    [minDist,index] = min(totalDist);
    distHistory(1,iter) = minDist;
    if minDist < globalMin
        globalMin = minDist;
        optRoute = pop(index,:);
    end
    distHistory(2,iter) = mean(totalDist);
    [maxDist, ~] = max(totalDist);
    distHistory(3,iter) = maxDist;
    
    %Only for demonstrating to a reviewer an issue
    if (iter>40)
        elitesMut = 0;
    end
    
    % Mutation Operator for the Genetic Algorithm
    for p = 1:elitesMut
        randIdx = ceil(1 + (popSize-1)*rand(1,tournaments)); % interval (a,b) with the formula r = a + (b-a).*rand(N,1).
        rtes = pop(randIdx(1:tournaments),:);
        dists = totalDist(randIdx(1:tournaments));
        [~,idx] = min(dists);
        best= rtes(idx,:);
        routeInsertionPoints = sort(ceil(n*rand(1,2)));
        I = routeInsertionPoints(1);
        J = routeInsertionPoints(2);
        % Mutate the Best
        tmpPop(1,:) = best;
        switch mutOper
            case 1 % Flip
                tmpPop(1,I:J) = tmpPop(1,J:-1:I);
            case 2 % Swap
                tmpPop(1,[I J]) = tmpPop(1,[J I]);
            case 3 % Slide
                tmpPop(1,I:J) = tmpPop(1,[I+1:J I]);
            otherwise % Do Nothing
        end
        newPop = [newPop; tmpPop];

    end
    
    if (elitesMut>0)
        %Mecanismo para mantener el tamaño de la población
        [r,~] = size(newPop);
        restPop = popSize - r;
    else
        restPop = popSize;
    end
    
    %Elitism Selection Mechanism
    if (elite>0)
        restPop = restPop - elite;
        tmpTotalDist = totalDist;
        p=0;
        while (p<elite)
            [~,index] = min(tmpTotalDist);
            tmpTotalDist(index) = +Inf;
            newPop = [newPop; pop(index,:)];
            p = p+1;
        end
    end
    
    tmpPop = zeros(1,n);
    while (restPop>0)
        randomOrder = randperm(popSize);
        firstParent = zeros(1,n);
        secondParent = zeros(1,n);
            randIdx = ceil(tournaments + (popSize-tournaments)*rand(1,1)); % interval (a,b) with the formula r = a + (b-a).*rand(N,1).
            rtes = pop(randomOrder(randIdx-tournaments+1:randIdx),:);
            dists = totalDist(randomOrder(randIdx-tournaments+1:randIdx));
            [~,idx] = min(dists);
            firstParent = rtes(idx,:);
            dists(1,idx) = +Inf;
            [~,idx] = min(dists);
            secondParent = rtes(idx,:);
            
        switch crossOper
            case 1 %Order Crossover (OX)
                tmpPop = OX(firstParent,secondParent);
                restPop = restPop - 1;
            case 2 %Cycle Crossover (CX)
                tmpPop = CX(firstParent,secondParent);
                restPop = restPop - 2;
                if (restPop<0)
                    aux = tmpPop(1,:);
                    tmpPop = [];
                    tmpPop = aux;
                end
            case 3 %Order Base Crossover (OBX)
                tmpPop = OBX(firstParent,secondParent);
                restPop = restPop - 1;
                
            otherwise % Do Nothing
        end
        newPop = [newPop; tmpPop];
    end
    
    pop = newPop;
end

%Translate the path to a Home to Home path i.e., 1,..,..,..,1
idx = find(optRoute==1);
len = length(optRoute);
cycles = len - idx + 1;
rte = circshift(optRoute,cycles);
if last_depot == 1
    %rte = optRoute([1:n 1]); %Closed path
    rte = [rte last_depot]; %Closed path
end

UGV_path = [V3(1,rte);V3(2,rte)];
UGV_path = UGV_path';

%fprintf('Result for irace=%g\n', minDist);
end


function child = OX(p1,p2)

    [~,c1] = size(p1);
    [~,c2] = size(p2);
    child(1,1:c1) = -1;
    
    %First
    randIdx1 = sort(ceil(1 + (c1-1)*rand(1,2)));
    child(1,randIdx1(1,1):randIdx1(1,2)) = p1(1,randIdx1(1,1):randIdx1(1,2));
    
    %Second
    for p = 1:c2
        k = find(child(1,:)==p2(1,p), 1);
        if (isempty(k))
            pos = 0;
            t = 1;
            while (pos==0)
                if (child(1,t)==-1)
                    pos = t;
                end
                t = t + 1;
            end
            child(1,pos) = p2(1,p);
        end
    end
    
    
end

function childs = CX(parent1,parent2)
    totalCycles = [];
    childs = [];
    [~,c] = size(parent1);
    completed = false;
    %Compute N cycles, till parent1 = []
    while (~completed)
        [cycleN, parent1, parent2] = cycling(parent1,parent2);
        
        totalCycles = [totalCycles; cycleN];

        count = 0;
        for h = 1:c
            if (parent1(1,h)==-1)
                count = count + 1;
            end
        end
        
        if count==c
            completed = true;
        end
    end
    
    %Filling the new chromosomes
    [r,n] = size(totalCycles);
    numCycles = r/2;
    
    c = 0;
    for i=1:numCycles
        p = mod(i,2);
        if p==1 
            for j=1:n
                if (totalCycles(i+c,j)~=-1)
                   childs(1,j) = totalCycles(i+c,j);
                end
                if (totalCycles(i+c+1,j)~=-1)
                   childs(2,j) = totalCycles(i+c+1,j);
                end
            end
        else 
            for j=1:n
                if (totalCycles(i+c,j)~=-1)
                   childs(2,j) = totalCycles(i+c,j);
                end
                if (totalCycles(i+c+1,j)~=-1)
                   childs(1,j) = totalCycles(i+c+1,j);
                end
            end
        end
        c = c + 1;
    end
    
end

function child = OBX(parent1,parent2)
    [~,c] = size(parent1);
    child(1,1:c) = -1;
    
    %Step 1: Random indexs
    set = ceil(1 + (c-1)*rand(1,1));
    
    selectedPositions = randperm(set);
    selectedinP1 = parent1(selectedPositions);
    
    %Step 2: Copy parent2 to child
    for i=1:c
        if (isempty(find(selectedinP1==parent2(1,i), 1)))
            child(1,i) = parent2(1,i);
        end
    end
    
    %Step 3: Copy parent1 to child
     for i=1:c
        if (~isempty(find(selectedinP1==parent1(1,i), 1)))
            for j=1:c
                if (child(1,j)==-1)
                    child(1,j) = parent1(1,i);
                    break;
                end
            end
        end
    end
    
end

function [cycle, newParent1, newParent2] = cycling(parent1,parent2)
   
    [~,c] = size(parent1);
    cycle(1:2,1:c) = -1;
    
    %Find first element of the new parent
    r = 1;
    while (parent1(1,r)==-1)
        r = r + 1;
    end
    cycle = find_cycle(cycle,parent1(1,r),parent1(1,r),parent1,parent2);
    
    %Update new parents
    newParent1 = [];
    for p = 1:c
        idx = find(cycle(1,:)==parent1(1,p));
        if isempty(idx) %
           newParent1 = [newParent1 parent1(1,p)];
        else
           newParent1 = [newParent1 -1];
        end
    end
    newParent2 = [];
    for p = 1:c
        idx = find(cycle(2,:)==parent2(1,p));
        if isempty(idx) %
           newParent2 = [newParent2 parent2(1,p)];
        else
           newParent2 = [newParent2 -1];
        end
    end
    
end

function [cycle] = find_cycle(cycle,currValue,initValue,parent1,parent2)
    
    idx = find(parent1==currValue);
    cycle(1,idx) = currValue;
    
    currValue = parent2(1,idx);
    cycle(2,idx) = currValue;
    
    t = length(cycle);
    
    if (t==1)
        initValue = cycle(1,idx);
    elseif (t>1)
        if(initValue~=parent2(1,idx)) %Cycle NOT found!
           cycle = find_cycle(cycle,currValue,initValue,parent1,parent2);
        end
    end

end

