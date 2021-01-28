function [minDist, g_o] = nga(chr, last_depot, paths_x, paths_y)

% Initialize default configuration of the genetic algorithm
popSize     = 50;
numIter     = 1e2;
minDist     = 0;
x           = paths_x';
y           = paths_y';

% Sanity Checks
popSize     = 4*ceil(popSize/4);
numIter     = max(1,round(real(numIter(1))));

[n,c] = size(x);

% Initialize the Population
pop = zeros(popSize,n);
g_o = zeros(1,n);
for k = 1:n
    pop(1,k) = 1;
end

for k = 1:popSize
    pop(k,:) = randi([1 c],1,n);
end

% Run the GA
globalMin = Inf;
totalDist = zeros(1,popSize);
tmpPop = zeros(4,n);
newPop = zeros(popSize,n);

for iter = 1:numIter
    % Evaluate Each Population Member. Calculate Total Distance following
    % the positions of the current chromosome
    for p = 1:popSize
        if (last_depot == 1)
            A = [x(chr(1),pop(p,1));y(chr(1),pop(p,1))];
            B = [x(chr(n),pop(p,n));y(chr(n),pop(p,n))];
            d = dmat(A,B); % Closed Path
        else
            d = 0; % Open Path
        end
        for k = 2:n
            A = [x(chr(k-1),pop(p,k-1));y(chr(k-1),pop(p,k-1))];
            B = [x(chr(k),pop(p,k));y(chr(k),pop(p,k))];
            d = d + dmat(A,B);
            %d = d + dmat(pop(p,k-1),pop(p,k));
        end
        totalDist(p) = d;
    end
    
    % Find the Best Route in the Population
    [minDist,index] = min(totalDist);
    if minDist < globalMin
        globalMin = minDist;
        g_o = pop(index,:);
    end
    
    % Genetic Algorithm Operators
    randomOrder = randperm(popSize);
    for p = 4:4:popSize
        rtes = pop(randomOrder(p-3:p),:);
        dists = totalDist(randomOrder(p-3:p));
        [~,idx] = min(dists);
        bestOf4Route = rtes(idx,:);
        routeInsertionPoints = sort(ceil(rand(1,2))); %+1 altera las posiciones de los cromosomas a mutar
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

minDist = globalMin;

end

function [ds] = dmat(A, B)
    ds = sqrt(((A(1,1)-B(1,1))^2) + ((A(2,1)-B(2,1))^2));
end

