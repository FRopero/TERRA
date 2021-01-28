%This function is an modified search algorithm of the real A*. The goal of
%this function is to find out which is the PATH WITH LESS RTH (Return To Home) OPERATIONS for the UAV
%starting in the vertex_node{path_x(1), path_y(1)} and finishing in the
%same node, visiting only one time each waypoint. Thus, it is a
%modification of the TSP, with the restriction of the UAV's energy.

%This algorithm is going to create a graph non-directed to achieve its
%goal. Each node of the graph is a matlab struct like this:

% struct {
%   struct parent
%   coordinates x, y
%   t = estimated time of the whole path
%   g = distance of a single trip
%   h = remaining nodes
%   f = t + h
%   Ec = total energy consumed to this point
% }

function [ rte, dist, time, stops] = search_uav_operations( subpath, uav_data, cfgParams)

path = round(subpath,4);
home = [path(1,1);path(2,1)];
stops = 0;

%N = nº total of nodes
[~,N] = size(path);
h = N - 1;

%Start Node: g = 0 : t = 0 : h = N-1 : f = t + h : Ec = 0
vertex_node = struct('parent',-1,'x',path(1,1),'y',path(2,1),'t',0,'g',0,'h',0,'f',0,'r',N-1);
vertex_node.h = computeH_LKH2(path,vertex_node,uav_data, cfgParams);
vertex_node.f = vertex_node.h;

sol_found = false;
expanded_nodes = 0;
dist = 0;
time = 0;
rte = [];

% The ordered list of currently discovered nodes that are not evaluated yet.
% Initially, only the start node is known.
OPEN = vertex_node;
CLOSE = [];

while (~isempty(OPEN))
    
    [current_node, OPEN] = find_lowest(OPEN);
    
    %Delete current_node from OPEN
    OPEN(1).f = Inf; %OPEN(1) = current_node
    [OPEN] = deleteNode(OPEN);
    
    %%Goal Node is HOME with h=0
    if (isSamePoint([current_node.x;current_node.y],home) && current_node.r == 0)
        sol_found = true;
        break
    end
    
    [SUCS] = expand_graph(current_node, path, vertex_node, uav_data, cfgParams); %List of successors
  
    [~,n_sucs] = size(SUCS);
    for i=1:n_sucs
        [in_closed] = is_visited(CLOSE, SUCS(i)); %Optimiza Distancia, quitarlo
        if (0==in_closed) %Optimiza Distancia, quitarlo
            [in_open] = is_visited(OPEN, SUCS(i));
            if (0==in_open)
                SUCS(i).t = +Inf;
                SUCS(i).parent = '';
                suc = SUCS(i);
            else
                suc = OPEN(in_open);
            end
            
            %Get the suc replicant in OPEN
            [OPEN] = UpdateVertex(OPEN,current_node,suc,uav_data,vertex_node);
            
        end
    end
    
    expanded_nodes = expanded_nodes + n_sucs;
    %[sanityCheck] = SanityCheck(OPEN);
    %Put current_node in list of nodes currently visited CLOSED
    CLOSE = [CLOSE current_node];
end

if (sol_found)
    [rte, dist] = reconstruct_path(current_node, rte, dist, path);
    [~,c] = size(rte);
    for i=2:c-1
        if (rte(i)==1)
            stops = stops + 1;
        end
    end
    time = current_node.t;
else
    disp('UAV Subpath solution NOT found!');
end

end

function [OPEN] = UpdateVertex(OPEN,current_node,suc,uav_data,vertex_node)
    d = EUC_2D_Distance([suc.x;suc.y],[current_node.x;current_node.y]);
    t = EUC_2D_Time(d,uav_data);
    
    if (isSamePoint([current_node.x;current_node.y],[vertex_node.x;vertex_node.y])) 
        t = uav_data.To + t ;
    else
        if (isSamePoint([suc.x;suc.y],[vertex_node.x;vertex_node.y]))  
            tc = uav_data.Tt - t;
            t = t + uav_data.Tl + tc;
        end
    end
    
    if (current_node.t + t < suc.t)
       suc.t = current_node.t + t;
       suc.parent = current_node;
       if (is_visited(OPEN, suc)>0)
           [OPEN] = remove(OPEN,suc);
       end
       suc.f = suc.t + suc.h;
       OPEN = [OPEN suc];
    end
end

function [OPEN] = remove(OPEN,node)
    list = OPEN;
    idx = 0;
    if (~isempty(list))
        list = struct2table(list);
        list = [list(:,2:3) list(:,6) list(:,8)];
        node = struct2table(node);
        node = [node(:,2:3) node(:,6) node(:,8)];
        [r,~] = size(list);
        for i=1:r
            if (list{i,:}==node{1,:})
                idx = i;
                break;
            end
        end
        OPEN(idx).f = Inf;
        [OPEN] = deleteNode(OPEN);
    end
end

function [nodeIdx] = is_visited(list, node)
    nodeIdx = 0;
    if (~isempty(list))
        list = struct2table(list);
        list = [list(:,2:3) list(:,6) list(:,8)];
        node = struct2table(node);
        node = [node(:,2:3) node(:,6) node(:,8)];
        [r,~] = size(list);
        for i=1:r
            if (list{i,:}==node{1,:})
                nodeIdx = i;
                break;
            end
        end
    end
end

function [in_list] = is_visited_Working(list, node)
    in_list = 0;
    if (~isempty(list))
        list = struct2table(list);
        list = [list(:,2:3) list(:,6)];
        node = struct2table(node);
        node = [node(:,2:3) node(:,6)];
        in_list = sum(ismember(list,node,'rows')) >=1;
    end
end

function [sol] = isSamePath(nodeA,nodeB)
    if ((class(nodeA)=='struct') & (class(nodeB)=='struct'))
        if (isSamePoint([nodeA.x;nodeA.y],[nodeB.x;nodeB.y]))
            sol = isSamePath(nodeA.parent,nodeB.parent);
        else
            sol = false;
        end
    else
        if ((class(nodeA)~='struct') & (class(nodeB)~='struct'))
            if (isSamePoint([nodeA.x;nodeA.y],[nodeB.x;nodeB.y]))
                sol = true;
            else
                sol = false;
            end
        else
            sol = false;
        end
    end
end

function [OPEN] = deleteNode(OPEN) 

L = [];
[~,c] = size(OPEN);
for i=1:c 
    if (OPEN(i).f ~= Inf)
        L = [L OPEN(i)];
    end
end
OPEN = L;

end

function [rte, d] = reconstruct_path(current_node, rte, d, path)

[~,c] = size(path);

for i=1:c
    if (isSamePoint([current_node.x;current_node.y],path(:,i)))
        rte = [i rte];
        if (class(current_node.parent)=='struct')
            d = d + EUC_2D_Distance([current_node.x;current_node.y],[current_node.parent.x;current_node.parent.y]);
            [rte, d] = reconstruct_path(current_node.parent, rte, d, path);
            break;
        end
    end
end

end

function [node, L] = find_lowest(OPEN)
%This function order the OPEN list from the node with the lowest F value to
%the node with the highest F value

[~,co] = size(OPEN);
L = [];
lowest = 1;
[~,cl] = size(L);

while(cl<co)
    for j=1:co
        if (OPEN(j).f~=Inf)
            node = OPEN(j);
            lowest = j;
            break
        end
    end
    if (class(node)=='struct')
        for i=1:co
            if (OPEN(i).f~=Inf)
                if (OPEN(i).f < node.f)
                    node = OPEN(i); %Catch the lowest node
                    lowest = i;
                end
            end
        end
        
        OPEN(lowest).f = Inf;
        L = [L node]; %Return value
        node = Inf;
        [~,cl] = size(L);
    end
end

node = L(1); %Return value

end

function [L] = expand_graph(current_node, path, vertex_node, uav_data, cfgParams)
[~,c] = size(path);
L = [];
for i=1:c
    suc = [path(1,i);path(2,i)];
    if (~isSamePoint([current_node.x;current_node.y],suc))% UAV cannot flight to the target that he is already in it!
        if (isSamePoint([vertex_node.x;vertex_node.y],suc))
            is_rel = false; %The UAV always can to return to home
        else
            [is_rel] = is_inpath(current_node, suc); %Only expand the suc if it is not already in the path
        end
        if (~is_rel)
            d_toSuc = EUC_2D_Distance([current_node.x;current_node.y],suc);
            t_toSuc = EUC_2D_Time(d_toSuc,uav_data);
            if (isSamePoint([current_node.x;current_node.y],[vertex_node.x;vertex_node.y]))
                g = uav_data.To + t_toSuc ; %+ uav_data.Tg
                t_toSuc = current_node.t + uav_data.To + t_toSuc ; %+ uav_data.Tg
            else
                g = current_node.g + t_toSuc;% + uav_data.Tg;
                if (isSamePoint(suc,[vertex_node.x;vertex_node.y]))
                    tc = uav_data.Tt - g;
                    t_toSuc = current_node.t + t_toSuc + uav_data.Tl + tc;
                else
                    t_toSuc = current_node.t + t_toSuc ; % + uav_data.Tg
                end
            end
            d_tohome = EUC_2D_Distance([vertex_node.x;vertex_node.y],suc);
            t_tohome = EUC_2D_Time(d_tohome,uav_data);
            if (uav_data.Tt >= floor(g + t_tohome + uav_data.Tl)) %The UAV has enough energy to go this successor node and return to home?
                suc = struct('parent',current_node,'x',suc(1,1),'y',suc(2,1),'t',t_toSuc,'g',g,'h',0,'f',0,'r',0);
                if (isSamePoint([suc.x;suc.y],[vertex_node.x;vertex_node.y]))
                    r = current_node.r;
                else
                    r = current_node.r - 1;
                end
                suc.r = r;
                suc.h = computeH_LKH2(path,suc,uav_data,cfgParams);
                suc.f = t_toSuc+suc.h;
                L = [L suc];
            end
        end
    end
end

end

function [is_rel] = is_inpath(current, suc)
%This function search if a node has been already visited in the path
if (class(current)=='struct')
    if (isSamePoint([current.x;current.y],suc))
        is_rel = true;
    else
        [is_rel] = is_inpath(current.parent,suc);
    end
else
    %The last, it will be current.parent == -1, it is not a struct, is a
    %double
    is_rel = false;
end

end

function [bool] = isSamePoint(A,B)
bool = A(1,1) == B(1,1) && A(2,1) == B(2,1);
end

function [h] = computeH(path,suc,uav_data)
[~, c] = size(path);
rnodes = [];
h = 0;
for i=2:c
    if (~is_inpath(suc,path(:,i)))
        rnodes = [rnodes path(:,i)];
    end
end
if (~isempty(rnodes))
    [~, c] = size(rnodes);
    for i=2:c
        h = h + EUC_2D_Time(EUC_2D_Distance(path(:,1),path(:,i)),uav_data);
    end
    h = (h/uav_data.Tt)*(uav_data.To+uav_data.Tl + 250); %Good Performance!
    h = round(h,0);
end
end

function [h] = computeH_LKH1(path,suc,uav_data)

%1. Estimate the ideal Tf with the Clarke-Wright heuristic (Convex-Hull) or LKH ¿Empezando en current node y terminando en base?
%2. Estimate the Travel Costs Tch to Home for recharging (mean(EveryRemNodetoHome)*2)
%3. Compute the optimal Number of Charging Stops So = (Tf/uav_data.Tt) for the ideal Tf
%4. Estimate the total charging time required for the ideal Tf
%   idx = kmeans(rnodes,So). Each cluster has the current node as starting point to compute the wasted battery
%   Tc = Tc_cluster1 + Tc_cluster2 ....
%5. h = Tf+[(Tch*So)*Tc]
end

function [h] = computeH_LKH2(path,suc,uav_data,cfgParams)

%1. Estimate the number of charging stops So = ceil(Sum(EveryRemNodetoHome)/uav_data.Tt);
[~, n] = size(path);
h = 0;
rnodes = [suc.x;suc.y];
costToHome = 0;
for i=2:n
    if (~is_inpath(suc,path(:,i)))
        rnodes = [rnodes path(:,i)];
        costToHome = costToHome + EUC_2D_Time(EUC_2D_Distance(path(:,1),path(:,i)),uav_data);
    end
end
So = int64(ceil(costToHome/uav_data.Tt)); %Pesimista (ceil) / Optimista (floor)

if (So>0)
    %2. Order nodes with LKH solver (Greedy TSP)
    [~,n] = size(rnodes);
    cost_matrix = [];
    for i=1:n
        for j=1:n
            cost_matrix(i,j) = EUC_2D_Time(EUC_2D_Distance(rnodes(:,i),rnodes(:,j)),uav_data);
        end
    end
    if (n>2)

        [TSPsolution,~] = LKH_TSP(cost_matrix,struct('CostMatrixMulFactor',1,'user_comment',''),'tsp_solution',cfgParams.LKHdir,cfgParams.LKHdir);
    else
        TSPsolution = [1 2 1];
        %TSPcost = sum(cost_matrix(1,:));
    end
    
    ordered_nodes = rnodes(:,TSPsolution(1,1:n));
    %3. Cluster nodes in 'So' clusters -> idx = kmeans(rnodes,So),, without considering 'suc'
    idx = kmeans(ordered_nodes(:,2:n)',So);
    
    %4. Compute Tf of each cluster from Current_Node to Home
    cluster_matrix = zeros(2,So);
    cluster_matrix(2,:) = ordered_nodes(1,1).* ones(So,1);
    cluster_matrix(3,:) = ordered_nodes(2,1).* ones(So,1);
    for i=2:n %¿n+1?
        cluster_matrix(1,idx(i-1)) = cluster_matrix(1,idx(i-1)) + EUC_2D_Time(EUC_2D_Distance(cluster_matrix(2:3,idx(i-1)),ordered_nodes(:,i)),uav_data);
        cluster_matrix(2:3,idx(i-1)) = ordered_nodes(:,i);
    end
    %Sum the cost to return to home
    Tch = mean(costToHome)*2; %Estimated travel cost to home
    for j=1:So
        %Tf = cluster_matrix(1,j);
        cluster_matrix(1,j) = cluster_matrix(1,j) + EUC_2D_Time(EUC_2D_Distance(cluster_matrix(2:3,j),path(:,1)),uav_data);
        Se = cluster_matrix(1,j)/uav_data.Tt; %Estimated stops in cluster j
        %h(j) = (Se*Tch)+cluster_matrix(1,j)+(Se*uav_data.Tt); %Estimated Cost for cluster j
        h(j) = (Se*Tch)+(2*cluster_matrix(1,j));
    end
    %6. Final H
    h = sum(h);
end

end

function [h] = computeH_LKH3(path,suc,uav_data)

%1. Estimate the number of charging stops So = ceil(Sum(EveryRemNodetoHome)/uav_data.Tt);
[~, n] = size(path);
h = 0;
rnodes = [suc.x;suc.y];
costToHome = 0;
for i=2:n
    if (~is_inpath(suc,path(:,i)))
        rnodes = [rnodes path(:,i)];
        costToHome = costToHome + EUC_2D_Time(EUC_2D_Distance(path(:,1),path(:,i)),uav_data);
    end
end
So = int64(floor(costToHome/uav_data.Tt)); %Pesimista (ceil) / Optimista (floor)

if (So>0)
    %2. Cluster nodes in 'So' clusters -> idx = kmeans(rnodes,So),, without considering 'suc'
    [~, n] = size(rnodes);
    idx = kmeans(rnodes(:,2:n)',So);
    [r,~] = size(idx);
    %3. Compute cost matrix of each cluster
    %Tf(1,:) = EUC_2D_Time(EUC_2D_Distance([suc.x;suc.y],path(:,1)),uav_data)*ones(So,1);%Init to cost Return to Home
    Tf(1,:) = zeros(1,So);
    for k=1:So
        kcluster = [suc.x;suc.y];
        for h=2:r+1
            if (idx(h-1)==k)
                kcluster = [kcluster rnodes(:,h)];
            end
        end
        [~,n] = size(kcluster);
        
        cost_matrix = [];
        for i=1:n
            for j=1:n
                cost_matrix(i,j) = EUC_2D_Time(EUC_2D_Distance(kcluster(:,i),kcluster(:,j)),uav_data);
            end
        end
        %4. Compute TSP of each cluster
        if (n>2)
            LKHdir = 'C:\\Users\\ferna\\OneDrive\\Universidad\\Doctorado\\Matlab\\TERRA\\EstimatingHeuristic\\';
            [TSPsolution,TSPcost] = LKH_TSP(cost_matrix,struct('CostMatrixMulFactor',1,'user_comment',''),'tsp_solution',LKHdir,LKHdir);
        else
            TSPsolution = [1 2 1];
            TSPcost = sum(cost_matrix(1,:));
        end
        %5. Compute Tf = TSPcost
        Tf(So) = Tf(So) + TSPcost - EUC_2D_Time(EUC_2D_Distance(kcluster(:,TSPsolution(n)),kcluster(:,1)),uav_data); %Remove last path to the path (we want to reach home)
        Tf(So) = Tf(So) + EUC_2D_Time(EUC_2D_Distance(kcluster(:,TSPsolution(n)),path(:,1)),uav_data); %Add Return to Home
    end
    
    %6. Sum the cost to return to home
    Tch = mean(costToHome)*2; %Estimated travel cost to home
    for j=1:So
        Se = Tf(j)/uav_data.Tt; %Estimated stops in cluster j
        h(j) = (Se*Tch)+(2*Tf(j));
    end
    
    %7. Final H
    h = sum(h);
    
end
end

function [d] = EUC_2D_Distance(last,next)
    d = sqrt(((last(1,1) - next(1,1))^2) + ((last(2,1) - next(2,1))^2));
    d = round(d,4);
end

function [t] = EUC_2D_Time(d,uav_data)
    t = round(d/uav_data.Vuav,0);
end
