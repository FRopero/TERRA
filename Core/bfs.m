function [current_node] = bfs(home,path,current_node,uav_data)
%h = 0;
OPEN = []; % FIFO Queue (TRY TO ORDER BY COST 't' IN THE SAME LVL)
CLOSE = []; %Set of visited nodes
%t = accumulated cost
%g = accumulated cost per trip
%r = remaining nodes to visit
current_node = struct('parent',current_node.parent,'x',current_node.x,'y',current_node.y,'t',0,'g',current_node.g,'r',current_node.r);

OPEN = [current_node OPEN];

while (~isempty(OPEN))
    
    %[current_node, OPEN] = find_lowest(OPEN);
    current_node = OPEN(1);
    [OPEN] = remove(OPEN,current_node);
    
    if (isSamePoint([current_node.x;current_node.y],home) && current_node.r == 0)
        %h = current_node.t;
        break
    end
    
    [SUCS] = expand_graph(current_node, path, home, uav_data); %List of successors
    
    [~,n_sucs] = size(SUCS);
    for i=1:n_sucs
        [in_closed] = is_visited(CLOSE, SUCS(i)); 
        if (0==in_closed) 
            [in_open] = is_visited(OPEN, SUCS(i));
            if (0==in_open)
                OPEN = [SUCS(i) OPEN];
            end
        end
    end
    
    CLOSE = [CLOSE current_node];
end

end


function [SUCS] = expand_graph(current_node, path, home, uav_data)

[~,c] = size(path);
SUCS = [];
for i=1:c
    suc = [path(1,i);path(2,i)];
    if (~isSamePoint([current_node.x;current_node.y],suc))% UAV cannot flight to the target that he is already in it!
        if (isSamePoint(home,suc)) 
            is_rel = false; %The UAV always can to return to home
        else
            [is_rel] = is_inpath(current_node, suc); %Only expand the suc if it is not already in the path
        end
        if (~is_rel)
            d_toSuc = EUC_2D_Distance([current_node.x;current_node.y],suc);
            t_toSuc = EUC_2D_Time(d_toSuc,uav_data); 
            if (isSamePoint([current_node.x;current_node.y],home)) 
                g = uav_data.To + t_toSuc ; %+ uav_data.Tg
                t_toSuc = current_node.t + uav_data.To + t_toSuc ; %+ uav_data.Tg 
            else
                g = current_node.g + t_toSuc + uav_data.Tg;
                if (isSamePoint(suc,home))  
                    tc = uav_data.Tt - g;
                    t_toSuc = current_node.t + t_toSuc + uav_data.Tl + tc;
                else
                    t_toSuc = current_node.t + t_toSuc ; % + uav_data.Tg
                end
            end
            d_tohome = EUC_2D_Distance(home,suc);
            t_tohome = EUC_2D_Time(d_tohome,uav_data); 
            
            if (uav_data.Tt >= floor(g + t_tohome + uav_data.Tl)) %The UAV has enough energy to go this successor node and return to home?
                suc = struct('parent',current_node,'x',suc(1,1),'y',suc(2,1),'t',t_toSuc,'g',g,'r',0);
                if (isSamePoint([suc.x;suc.y],home))
                    suc.r = current_node.r;
                else
                    suc.r = current_node.r - 1;
                end
                SUCS = [SUCS suc];
            end
        end
    end
end

end

function [nodeIdx] = is_visited(list, node)
    nodeIdx = 0;
    if (~isempty(list))
        list = struct2table(list);
        list = [list(:,2:4) list(:,6)];
        node = struct2table(node);
        node = [node(:,2:4) node(:,6)];
        [r,~] = size(list);
        for i=1:r
            if (list{i,:}==node{1,:})
                nodeIdx = i;
                break;
            end
        end
    end
end

function [bool] = isSamePoint(A,B)
bool = A(1,1) == B(1,1) && A(2,1) == B(2,1);
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

function [OPEN] = remove(OPEN,suc)
list = OPEN;
if (~isempty(list))
    list = struct2table(list);
    list = [list(:,2:4) list(:,6)];
    suc = struct2table(suc);
    suc = [suc(:,2:4) suc(:,6)];
    [~, Locb] = ismember(list,suc,'rows');
    idx = find(Locb,1);
    OPEN(idx).f = Inf;
    [OPEN] = deleteNode(OPEN);
    
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

function [d] = EUC_2D_Distance(last,next)
    d = sqrt(((last(1,1) - next(1,1))^2) + ((last(2,1) - next(2,1))^2));
    d = round(d,4);
end

function [t] = EUC_2D_Time(d,uav_data)
    t = round(d/uav_data.Vuav,0);
end