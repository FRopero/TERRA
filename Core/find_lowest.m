function [node, L] = find_lowest(OPEN)
%This function order the OPEN list from the node with the lowest F value to
%the node with the longest F value

[~,co] = size(OPEN);
L = [];
[~,cl] = size(L);
lowest = 1;

while(cl<co)
    for j=1:co
        if (OPEN(j).f~=Inf)
            node = OPEN(j);
            lowest = j;
            break
        end
    end
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
    [~,cl] = size(L);
    node = Inf;
end
    
   node = L(1); %Return value

end