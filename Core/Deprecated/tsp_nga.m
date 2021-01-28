function [ rte, minDist, optGO] = tsp_nga( paths_x, paths_y, last_depot)

% Initialize default configuration of the genetic algorithm
popSize     = 50;
numIter     = 1e4;
minDist     = 0;

% Sanity Checks
popSize     = 4*ceil(popSize/4);
numIter     = max(1,round(real(numIter(1))));

% Initialize the Population
[n,~] = size(paths_x');
pop = zeros(popSize,n);
m_go = zeros(popSize,n);

for k = 1:popSize
    pop(k,:) = randperm(n);
    p = 1;
    q = 1;
    %Ensure the home waypoint constraint
    for i = 1:n
        if pop(k,i)==1
            p = i;
        end
    end
    pop(k,p) = pop(k,1);
    pop(k,1) = 1;
    
    %Ensure the latest UGV depot waypoint constraint
    if last_depot ~= 1
        for i = 1:n
            if pop(k,i)==last_depot
                q = i;
            end
        end
        pop(k,q) = pop(k,n);
        pop(k,n) = last_depot;
    end
end

% Run the GA
globalMin = Inf;
totalDist = zeros(1,popSize);
tmpPop = zeros(4,n);
newPop = zeros(popSize,n);

for iter = 1:numIter
    % Evaluate Each Population Member.
    for p = 1:popSize
         [dist, g_o] = nga(pop(p,:),last_depot,paths_x,paths_y);
         totalDist(p) = dist;
         m_go(p,:) = g_o;
    end
    
    % Find the Best Route in the Population
    [minDist,index] = min(totalDist);
    if minDist < globalMin
        globalMin = minDist;
        optRoute = pop(index,:);
        optGO = m_go(index,:);
    end
    
    % Genetic Algorithm Operators
    randomOrder = randperm(popSize);
    for p = 4:4:popSize
        rtes = pop(randomOrder(p-3:p),:);
        dists = totalDist(randomOrder(p-3:p));
        [~,idx] = min(dists);
        bestOf4Route = rtes(idx,:);
        routeInsertionPoints = sort(ceil((n-2)*rand(1,2))+1); %+1 altera las posiciones de los cromosomas a mutar
        I = routeInsertionPoints(1);
        J = routeInsertionPoints(2);
        for k = 1:4 % Mutate the Best to get Three New Routes
            tmpPop(k,:) = bestOf4Route;
            switch k
                case 2 % Flip
                    tmpPop(k,I:J) = tmpPop(k,J:-1:I);
                case 3 % Swap
                    tmpPop(k,[I J]) = tmpPop(k,[J I]);
                case 4 % Slide
                    tmpPop(k,I:J) = tmpPop(k,[I+1:J I]);
                otherwise % Do Nothing
            end
        end
        newPop(p-3:p,:) = tmpPop;
    end
    pop = newPop;
    
end

if last_depot == 1
    rte = [optRoute last_depot]; %Closed path
    optGO = [optGO optGO(1)];
else
    rte = optRoute; %Open path
end

end

