%function [ rte, dist, stops] = search_uav_3Dpath( subpath, uav_data)
%This function is an modified search algorithm of the real A* algorithm. The goal of
%this function is to find out which is the SHORTEST PATH for the UAV
%starting in the start_node{path_x(1), path_y(1)} and finishing in the
%same node, visiting only one time each waypoint. Thus, it is a
%modification of the TSP, with the restriction of the UAV's energy.
%NOTE: Every node have the same altitude and it does not compute obstacles
%avoidance

%This algorithm is going to create a graph non-directed to achieve its
%goal. Each node of the graph is a matlab struct like this:

% struct {
%   struct parent
%   coordinates x, y
%   g = dE (dE = Euclidean distance accumulated from Home) 
%   k = dE accumulated from the last stop at Home
%   h = nº of remaining visited nodes
%   f = g + h
% }

function [ rte, dist, stops] = search_uav_path( subpath, uav_data)

Radius = uav_data.R;
path = subpath;
home = [path(1,1);path(2,1)];
stops = 0;
%N = nº total of nodes
[~,N] = size(path);
h = N - 1;
%Start Node: g = 0 : k = 0 : h = N-1 : f = g + h

start_node = struct('parent',-1,'x',path(1,1),'y',path(2,1),'g',0,'k',0,'h',h,'f',h);

sol_found = false;
dist = 0;
rte = [];

% The ordered list of currently discovered nodes that are not evaluated yet.
% Initially, only the start node is known.
OPEN = start_node;
CLOSED = [];

while (~isempty(OPEN))
    
    [current_node, OPEN] = find_lowest(OPEN);
    
    %Delete current_node from OPEN
    OPEN(1).f = Inf; %OPEN(1) = current_node
    [OPEN] = deleteNode(OPEN);
    
    %%Goal Node is HOME with h=0
    if ((current_node.x == home(1,1)) && (current_node.y == home(2,1)) && (current_node.h == 0))
        %disp('UAV Subpath solution found!');
        sol_found = true;
        break
    end
    
    [SUCS] = expand_graph(current_node, path, Radius, start_node); %List of successors
    %Check if those successors has been visited yet
    [~,n_sucs] = size(SUCS);
    for i=1:n_sucs
        [in_closed] = is_visited(CLOSED, SUCS(i));
        if (~in_closed)
            [in_open] = is_visited(OPEN, SUCS(i));
            if (~in_open)
                SUCS(i).g = +Inf;
                SUCS(i).parent = '';
            end
            [OPEN] = UpdateVertex(OPEN,current_node,SUCS(i));
        end
    end
    
    %Put current_node in list of nodes currently visited CLOSED
    CLOSED = [CLOSED current_node];
end

if (sol_found)
    [rte, dist] = reconstruct_path(current_node, rte, dist, path);
    [~,c] = size(rte);
    for i=2:c-1
        if (rte(i)==1)
            stops = stops + 1;
        end
    end
else
    disp('UAV Subpath solution NOT found!');
end

end

function [OPEN] = UpdateVertex(OPEN,current_node,suc)
   d = sqrt(((suc.x-current_node.x)^2) + ((suc.y-current_node.y)^2));
   if (current_node.g + d < suc.g)
       suc.g = current_node.g + d;
       suc.parent = current_node;
       [in_open] = is_visited(OPEN, suc);
       if (in_open)
           [OPEN] = remove(OPEN,suc);
       end
       suc.f = suc.g + suc.h;
       OPEN = [OPEN suc];
   end
end

function [OPEN] = remove(OPEN,suc)
    list = OPEN;
    if (~isempty(list))
        list = struct2table(list);
        list = [list(:,2:3) list(:,6)];
        suc = struct2table(suc);
        suc = [suc(:,2:3) suc(:,6)];
        [~, Locb] = ismember(list,suc,'rows');
        idx = find(Locb,1);
        OPEN(idx).f = Inf;
        [OPEN] = deleteNode(OPEN);
        
    end
end

function [in_list] = is_visited(list, node)
    in_list = 0;
    if (~isempty(list))
        list = struct2table(list);
        list = [list(:,2:3) list(:,6)];
        node = struct2table(node);
        node = [node(:,2:3) node(:,6)];
        in_list = sum(ismember(list,node,'rows')) >=1;
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

function [rte, dist] = reconstruct_path(current_node, rte, dist, path)

[~,c] = size(path);

    for i=1:c
        if (current_node.x == path(1,i) && current_node.y == path(2,i))
            rte = [i rte];
            if (class(current_node.parent)=='struct')
                dist = dist + sqrt(((current_node.x - current_node.parent.x)^2) + ((current_node.y - current_node.parent.y)^2));
                [rte, dist] = reconstruct_path(current_node.parent, rte, dist, path);
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

function [L] = expand_graph(current_node, path, Radius, start_node)

[~,c] = size(path);
L = [];
home = [start_node.x;start_node.y];
for i=1:c
    suc = [path(1,i);path(2,i)];
    %If current node is not the root of the tree, and the successor it is
    %not this same point
    if (floor(current_node.x) ~= floor(suc(1,1)) || floor(current_node.y) ~= floor(suc(2,1)))
    %if (current_node.x ~= suc(1,1) || current_node.y ~= suc(2,1))
        %Check if node(i) is a relative of the current node
        if (suc(1,1)==start_node.x && suc(2,1)==start_node.y)
            is_rel = false;
        else
            [is_rel] = is_relative(current_node, suc);
        end
        if (~is_rel)
            if (current_node.x == start_node.x && current_node.y == start_node.y)
                k = sqrt(((path(1,i)-current_node.x)^2) + ((path(2,i)-current_node.y)^2));
            else
                k = current_node.k + sqrt(((path(1,i)-current_node.x)^2) + ((path(2,i)-current_node.y)^2));
            end
            g = current_node.g + sqrt(((path(1,i)-current_node.x)^2) + ((path(2,i)-current_node.y)^2));
            d_tohome = sqrt(((path(1,i)-start_node.x)^2) + ((path(2,i)-start_node.y)^2));
            %The UAV has enough energy to go this successor node and return to home
            if (Radius*2 >= floor(k + d_tohome))
                if (suc==home)
                    h = current_node.h;
                else
                    h = current_node.h - 1;
                end
                %No relative, so, we add the successor to the OPEN List
                suc = struct('parent',current_node,'x',path(1,i),'y',path(2,i),'g',g,'k',k,'h',h,'f',g+h);
                L = [L suc];
            end
        end
    end
end

end

function [is_rel] = is_relative(current, suc)
%This function search if a node has been already visited in the path
if (class(current)=='struct')
    if (current.x == suc(1,1) && current.y == suc(2,1))
        is_rel = true;
    else
        [is_rel] = is_relative(current.parent,suc);
    end
else
    %The last, it will be current.parent == -1
    is_rel = false;
end

end
