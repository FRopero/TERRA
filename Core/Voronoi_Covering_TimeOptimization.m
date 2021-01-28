function [ V1, covered_Tps, figV] = Voronoi_Covering_TimeOptimization( problem_params, uav_data, cfgParams)
figV = [];

%Algorithm variables
V1 = [];
uncovered_Tps = problem_params.T;
covered_Tps = []; 

trusted_Vertices = [];
used_vertices = [];
artificial_Vertices = problem_params.Home;

duplicates = 0;

if (cfgParams.printResults)
    vis = 'on';
else
    vis = 'off';
end

% figV = [];
% %First iteration plot
% figInit = figure('visible',vis);
% figInit.Name = 'Voronoi_FirstIteration';
% if (cfgParams.saveResults)
%     figV = [figV figInit];
% end
% [vx,vy] = voronoi(wp_nc(1,:),wp_nc(2,:));
% hold on
% plot(wp_nc(1,:),wp_nc(2,:),'blue.',vx,vy,'black-',vx,vy,'green+');
% [~,c]= size(vx);
% theta = linspace(0,2*pi);
% for i=1:c
%     c_x(i,:) = problem_params.R*cos(theta) + vx(i);
%     c_y(i,:) = problem_params.R*sin(theta) + vy(i);
%     plot(c_x(i,:),c_y(i,:),'r:');
% end
% hold off
% axis equal
% title('First Voronoi Iteration');


while ( ~isempty(uncovered_Tps) )
    
    % [STEP 1] - Calculate the Voronoi's Diagram over the WPm not covered
     nearest_vertices = [];
     
    if (~isempty(trusted_Vertices)) %For n iterations
        %Search and catch the nearest vertice to each wp_nc and add to
        %near_v
        [~,c_wp] = size(uncovered_Tps);
        [~,c_vx] = size(trusted_Vertices);
        
        for i=1:c_wp
            idx = 1;
            tmin = EUC_2D_Time(EUC_2D_Distance(uncovered_Tps(:,i),trusted_Vertices(:,1)),uav_data);
            for j=2:c_vx
                t = EUC_2D_Time(EUC_2D_Distance(uncovered_Tps(:,i),trusted_Vertices(:,j)),uav_data);
                if (t < tmin)
                    tmin = t;
                    idx = j;
                end
            end
            nearest_vertices = [nearest_vertices trusted_Vertices(:,idx)];
        end
        
        %Delete duplicates
        nearest_vertices = unique(nearest_vertices.','rows').';  
        
        %Detect duplicates and so, infinite loop in voronoi diagram
        duplicates = isDuplicates(nearest_vertices,used_vertices);
        used_vertices = [used_vertices nearest_vertices];
    end

    %Add these vertices to the next voronoi iteration
    VOR_Vertices = [uncovered_Tps nearest_vertices];
    
    %We can't do voronoi if we have 2 points (one wp and its nearest vertice)
    %or Only 2 wp and 1 vértice, and the distance between the WPs and their common nearest
    %vertice, is less than the R. This situation cause an infinite loop.
    if ((length(VOR_Vertices(1,:))==2)||(duplicates))
        %m is the intermediate point between the WPs. Each iteration,
        %this vertice is nearest to the first wp
        m_x = (VOR_Vertices(1,1) + VOR_Vertices(1,2)) / 2;
        m_y = (VOR_Vertices(2,1) + VOR_Vertices(2,2)) / 2;
        d = EUC_2D_Distance(VOR_Vertices(:,1),[m_x;m_y]);
        %R = (uav_data.Tt - uav_data.To - uav_data.Tl)/2; %m
        %ttrip = EUC_2D_Time(d,uav_data);
        R = problem_params.R; %m
        %Rtrip = uav_data.Vuav*(uav_data.To + EUC_2D_Time(d,uav_data)*2 + uav_data.Tl); %m
        while (d > R)
            m_x = (m_x + VOR_Vertices(1,1)) / 2;
            m_y = (m_y + VOR_Vertices(2,1)) / 2;
            d = EUC_2D_Distance(VOR_Vertices(:,1),[m_x;m_y]);
            %ttrip = EUC_2D_Time(d,uav_data);
        end
        artificial_Vertices = [artificial_Vertices [m_x;m_y]];
    else
        %Create the voronoi diagram for the n iteration
        [vx,vy] = voronoi(VOR_Vertices(1,:),VOR_Vertices(2,:));
        
        %Delete duplicates vertices
        trusted_Vertices = [vx(1,:); vy(1,:)];
        trusted_Vertices = unique(trusted_Vertices.','rows').'; 
    end
    
    %Add artificial vertices created in the last iteration
    trusted_Vertices  = [artificial_Vertices trusted_Vertices];
    
    %Compute the new set of uncovered waypoints
    [wx,wy] = size(uncovered_Tps);
    [~,ct] = size(trusted_Vertices);
    for i=1:wy
        R = problem_params.R; %m
        %ttrip_max = uav_data.Tt;
        best_v = -1;
        tp_covered = -1;
        for j=1:ct
            if (uncovered_Tps(1,i)~=-1) 
                %Time Cost from base = To+(2*Tf)+Tl
                d = EUC_2D_Distance(uncovered_Tps(:,i),trusted_Vertices(:,j));
                %ttrip = uav_data.To + (2*EUC_2D_Time(d,uav_data)) + uav_data.Tl;
                %Rtrip =  uav_data.Vuav*(uav_data.To + EUC_2D_Time(d,uav_data)*2 + uav_data.Tl) ; %m
                if (d < R && best_v ~= 1) %%Si el home ya lo ha cogido, es la mejor opción
                    R = d;
                    best_v = j;
                    tp_covered = i;
                end
            end
        end
        if (best_v ~= -1)% The Waypoint uncovered_Tps(x,i) was covered by unless one vertice!!
            covered_Tps = [covered_Tps uncovered_Tps(:,tp_covered)];
            uncovered_Tps(:,tp_covered) = [-1;-1];
            V1 = [V1 trusted_Vertices(:,best_v)];
        end
    end
    
    del = 0; %Number of wp covered this iteration
    for i=1:wy
        if (uncovered_Tps(1,i)==-1)
            del = del + 1;
        end
    end
    
    %Update uncovered_Tps Table for the next iteration
    if (del>0)
        tmp_wp_nc = uncovered_Tps;
        uncovered_Tps = zeros(wx,wy-del);
        for i=1:wx
            t_y = 0;
            for j=1:wy
                if (tmp_wp_nc(i,j)~=-1)
                    t_y = t_y + 1;
                    uncovered_Tps(i,t_y) = tmp_wp_nc(i,j);
                end
            end
        end
    end
    
end

%Plot Results
figLast = figure('visible',vis);
figLast.Name = 'Voronoi_Solution';
if (cfgParams.saveResults)
    figV = [figV figLast];
end
hold on
plot(problem_params.T(1,:),problem_params.T(2,:),'blue.',V1(1,:),V1(2,:),'green+');
[~,c]= size(V1);
theta = linspace(0,2*pi);
for i=1:c
    c_x(i,:) = problem_params.R*cos(theta) + V1(1,i);
    c_y(i,:) = problem_params.R*sin(theta) + V1(2,i);
    plot(c_x(i,:),c_y(i,:),'r:');
end
hold off
axis equal
title('Voronoi Solution');

end

% Search duplicates between two vectors
function [bool] = isDuplicates(lA,lB)
    bool = false;
    %Count duplicates and so, infinite loop in voronoi diagram
    [~,lv] = size(lA);
    [~,lu] = size(lB);

    for i=1:lv
        for j=1:lu
            if (lA(:,i)==lB(:,j))
                bool = true;
                break
            end
        end
    end
end

function [d] = EUC_2D_Distance(last,next)
    d = sqrt(((last(1,1) - next(1,1))^2) + ((last(2,1) - next(2,1))^2));
    d = round(d,4);
end

function [t] = EUC_2D_Time(d,uav_data)
    t = round(d/uav_data.Vuav,0);
end

