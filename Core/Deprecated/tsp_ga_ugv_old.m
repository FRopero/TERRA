function [ rte, minDist, xy] = tsp_ga_ugv(x,y, path_x, path_y, showResults, type, saveResults, fullname)

% Initialize default configuration of the genetic algorithm
%xy          = 10*rand(50,2);
dmat        = [];
popSize     = 200;
numIter     = 100;
showProg    = false;
showResult  = showResults;
showWaitbar = false;
minDist     = 0;
xy = [];
for i=1:length(path_x)
    xy(i,1) = path_x(i);
end
for i=1:length(path_y)
    xy(i,2) = path_y(i);
end

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
showProg    = 0;
showResult  = false;
showWaitbar = logical(showWaitbar(1));

% Initialize the Population
pop = zeros(popSize,n);
%pop(1,:) = (1:n);

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
    last_depot = 1;
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
distHistory = zeros(1,numIter);
tmpPop = zeros(4,n);
newPop = zeros(popSize,n);
if showProg
    %fig2 = figure('Name','TSP_GA | Current Best Solution','Numbertitle','off','Position', [2550, 300, 470, 350]);
    fig2 = figure('Name','TSP_GA | Current Best Solution','Numbertitle','off');
    hAx = gca;
end
if showWaitbar
    hWait = waitbar(0,'Searching for near-optimal solution ...');
end

for iter = 1:numIter
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
    distHistory(iter) = minDist;
    if minDist < globalMin
        globalMin = minDist;
        optRoute = pop(index,:);
        if showProg
            % Plot the Best Route
            set(0,'CurrentFigure',fig2);
            rte = optRoute([1:n 1]);
            if dims > 2, plot3(hAx,xy(rte,1),xy(rte,2),xy(rte,3),'r.-');
            else plot(hAx,xy(rte,1),xy(rte,2),'r.-'); end
            title(hAx,sprintf('Total Distance = %1.4f, Iteration = %d',minDist,iter));
            drawnow;
        end
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
    
    % Update the waitbar
    if showWaitbar && ~mod(iter,ceil(numIter/325))
        waitbar(iter/numIter,hWait);
    end
    
end
if showWaitbar
    close(hWait);
end

if showResult
    % Plots the GA Results
    figure('Name','TSP_GA | Results','Numbertitle','off');
    subplot(2,2,1);
    pclr = ~get(0,'DefaultAxesColor');
    if dims > 2, plot3(xy(:,1),xy(:,2),xy(:,3),'.','Color',pclr);
    else plot(xy(:,1),xy(:,2),'.','Color',pclr); end
    title('City Locations');
    subplot(2,2,2);
    imagesc(dmat(optRoute,optRoute));
    title('Distance Matrix');
    subplot(2,2,3);
    rte = optRoute([1:n 1]);
    if dims > 2, plot3(xy(rte,1),xy(rte,2),xy(rte,3),'r.-');
    else plot(xy(rte,1),xy(rte,2),'r.-'); end
    title(sprintf('Total Distance = %1.4f',minDist));
    subplot(2,2,4);
    plot(distHistory,'b','LineWidth',2);
    title('Best Solution History');
    set(gca,'XLim',[0 numIter+1],'YLim',[0 1.1*max([1 distHistory])]);
end

%Plots results on the main figure
%set(0,'CurrentFigure',fig1)

%subplot(3,2,5);
if last_depot == 1
    %rte = optRoute([1:n 1]); %Closed path
    rte = [optRoute last_depot]; %Closed path
else
    rte = optRoute; %Open path
end

% if dims > 2, plot3(xy(rte,1),xy(rte,2),xy(rte,3),'r.-');
% else plot(xy(rte,1),xy(rte,2),'black-',home_x,home_y,'black*',sol_x,sol_y,'black+',avg_x,avg_y,'r+'); end
%     title(sprintf('SP with a Genetic Algorithm for TSP [f = %1.4f]',minDist));

% START Draw %
if (showResults)
    figure;
    plot(xy(rte,1),xy(rte,2),'black-+',path_x(1,1),path_y(1,1),'black*',x,y,'blue.');
    l = length(type);
    if (10~=l)
        title(sprintf('[STEP 4] - (%s) for the UGV-TSP [fugv = %1.3f km]',type,minDist));
    else
        title(sprintf('[STEP 3&4] - (%s) for the UGV-TSP [fugv = %1.3f km]',type,minDist));
    end
    axis equal;
    if (saveResults)
        tp = char(type);
        dir = strcat(fullname,'\','ugv_tsp_',tp,'.fig');
        saveas(gcf,dir,'fig');
    end
end
% END Draw %

end

