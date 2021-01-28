function [uav_path1, uav_path2, distance, time, stops] = uav_compute_path(wp_c, scp_table, SolL, Vsol, cfgParams, ugv_path, uav_data)
%uav_compute_path(xyS, routeS, home_x, home_y, x, y, wp_c, scp_table, SolL, v_op, showResults, ugv_path, uav_data, type, saveResults, fullname)
uav_path1 = [];
uav_path2 = [];
stops = 0;

distance = 0;
time = 0;
wp_c_aux = wp_c;
[r,c] = size(scp_table);

for i=1:c
    for k=1:length(SolL)
        if (i==SolL(k))
            subpath = [Vsol(1,i);Vsol(2,i)]; %Catch the vertice of the solution
            for j=1:r
                if (scp_table(j,i)==1 && wp_c_aux(1,j) ~= -1) %Waypoint j
                    subpath = [subpath [wp_c(1,j);wp_c(2,j)]];
                    wp_c_aux(1,j) = -1;
                end
            end
            
            % DISTANCE BASED SEARCHING
            %[rteP, dis, st] = search_uav_path( subpath, uav_data); 
            
            % TIME BASED SEARCHING
            [rteT, dis, t, st] = search_uav_operations(subpath, uav_data, cfgParams); 
            time = time + t;
            
            distance = distance + dis;
            stops = stops + st;
            %UAV Path for Distance-Based Searching or Time-Based Searching
            [ru,~] = size(ugv_path);
            for z=1:ru
                if (ugv_path(z,:)==Vsol(:,i)')
                    xy_uav = [subpath(1,:)',subpath(2,:)'];
                    % DISTANCE BASED SEARCHING
                    %uav_path1 = [ uav_path1 struct('Coordinates',[xy_uav(rteP,1) xy_uav(rteP,2)])];
                    % TIME BASED SEARCHING
                    uav_path2 = [ uav_path2 struct('Coordinates',[xy_uav(rteT,1) xy_uav(rteT,2)])];
                end
            end
            
            
        end
    end
end

%%% START Drawing
if (cfgParams.printResults)
    %figure(fig2);
    %xy_uav = [subpath(1,:)',subpath(2,:)'];
    %plot(xy_uav(rteP,1),xy_uav(rteP,2),'blue:');
    %figure(fig3);
    %xy_uav = [subpath(1,:)',subpath(2,:)'];
    %plot(xy_uav(rteP,1),xy_uav(rteP,2),'red:');
%               figure;
%               [~,c] = size(rte);
%               plot(1:c,rte-1,'black-+');
%               title(sprintf('UAV SubPath ID[%d] - Fc = %0.4f\n',t,minD));
%               t = t + 1;
end
%%% STOP Drawing

end

