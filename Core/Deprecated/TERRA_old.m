function [data_sol, path_sol] = TERRA(showResults, home, xy, last_depot, uav_data, ugv_data, saveResults, fullname)
% This function executes the whole solution designed and developed to solve
% the Fuel Constrained UAV Refueling Problem - Moving Refueling Station
% (FCURP-MRS). In fact, the problem presented here is a version of the
% FCURP-MRS, with some differences. The aim of this problem is to visit all
% the waypoints with the UAV, supporting it with UGV as a recharging
% station for the UAV battery.

% INPUTS %
% showResults :Boolean to show the figures of a single result.
% home        :Coordinates of the home location of the both robots.
% xy          :Coordinates of the waypoints selected in the problem.
% last_depot  :Integer with the last depot of the UGV in its path. Usually, it is equal to home.
% uav_data    :Data of the UAV
% ugv_data    :Data of the UGV
% saveResults :Boolean to save figures in a daily log
% fullname    :Directory where save the figures

% OUTPUTS %
% data_sol    :list of structures with the cost of the paths of each solution.
%   struct {
%      f_ugv_d       :Cost of the UGV as distance travelled in meters.
%      f_ugv_t       :Cost of the UGV as time in seconds.
%      f_uav_d       :Cost of the UAV as distance travelled in meters. (SEARCH_UAV_PATH)
%      f_uav_t       :Cost of the UAV as time in seconds. (SEARCH_UAV_PATH)
%      ftotal_d      :Total cost in meters.
%      ftotal_t      :Total cost in seconds.
%   }
% path_sol    :paths of both UGV-UAV to reach the solutions computed in the algorithm.

%Initialize outputs
path_sol = [];
data_sol = [];

%UAV Data
R = uav_data.R;
theta = linspace(0,2*pi);

%HOME Location / Init Depot
home_x = home(1,1);
home_y = home(2,1);

%Waypoints Generated
x = xy(1,:);
y = xy(2,:);

% START Drawing 
if (showResults)
    fig1 = figure('Name','TERRA | FCURP-MRS | Results','Numbertitle','off','Position', [1600, 40, 900, 750]);
    subplot(2,2,1);
    plot(x,y,'blue.',home_x,home_y,'black*');
    title('New FCURP-MRS generated');
    axis equal;
else
    fig1 = [];
end
% END Drawing 

% [STEP 1] - Calculate the Voronoi's Diagram over the WPm not covered
[ vx_sol, vy_sol, wp_c ] = Voronoi_Covering_Optimization_nearVertex( x,y,home_x,home_y, R, theta, fig1, showResults);

% [STEP 2] - Calculate with the Greedy Algorithm Set Cover Problem, the minimun set of
%vértices, whose circles contains our universe of waypoints
[SolC,SolL,scp_table,sol_x,sol_y] = greedy_scp(x,y,home_x,home_y, vx_sol, vy_sol, wp_c, R, theta, fig1, showResults, saveResults, fullname);

%Find out if home is part of the vertices solution, in that case, put into
%the first position
enc = false;
for i=1:length(sol_x)
    if sol_x(i) == home_x
        enc = true;
        pos = i;
    end
end
if enc
    path_x = sol_x;
    path_x(pos) = path_x(1);
    path_x(1) = home_x;
    path_y = sol_y;
    path_y(pos) = path_y(1);
    path_y(1) = home_y;

else
    path_x = [home_x sol_x];
    path_y = [home_y sol_y];
end


%%%% [STEP 3,4 and 5] 
% - Compute shortest path among sol_x and sol_y with a genetic algorithm for the TSP
% - Optimizing the distances between the vertices (path_x, path_y) of the
%   solution and the average center computed. To optimize the cost of
%   travel of the UGV
% - Compute the shortest path for the UAV, to achieve the
% waypoints. For each subpath

[route, ugv_distance, xy] = tsp_ga_ugv(x, y, path_x, path_y, last_depot, showResults, 'Without GO', saveResults, fullname);

%UGV Path
ugv_path = [xy(route,1) xy(route,2)];
ugv_time = ugv_distance/ugv_data.Vugv;

%UAV Path
[uav_path1, uav_path2, uav_distance, uav_time, stops] = uav_compute_path(xy, route, home_x, home_y, x, y, wp_c, scp_table, SolL, [vx_sol;vy_sol], showResults, ugv_path, uav_data, 'Without GO', saveResults, fullname);

%1st Solution without GOA
total_time = ugv_time + uav_time;
total_distance = ugv_distance + uav_distance;

s = struct('f_ugv_d',ugv_distance,'f_ugv_t',ugv_time,'f_uav_d',uav_distance,'f_uav_t',uav_time,'ftotal_d',total_distance,'ftotal_t',total_time,'stops',stops);
data_sol = [data_sol ; s];
[ p_sol ] = build_matrix_solution( ugv_path, uav_path1, uav_path2 );
path_sol = [path_sol ; struct('path_solution',p_sol)];

%GOA Application
avg_x = [mean(path_x) home_x median(path_x)];
avg_y = [mean(path_y) home_y median(path_y)];
G = string({'GravityCenter','HomeCenter','MedianCenter'});

for i=1:length(avg_x)

    %GOA
    [vx_opt, vy_opt, figG, vres] = Gravitational_Optimization( home_x, home_y, x, y, vx_sol, vy_sol, avg_x(i), avg_y(i), wp_c, SolL, scp_table, R, showResults);
    
    %Genetic Algorithm to UGV Path
    [routeG, ugv_distance, xyG] = tsp_ga_ugv(x, y, [home_x vx_opt], [home_y vy_opt], last_depot, showResults, G(i), saveResults, fullname);
    ugv_path = [xyG(routeG,1) xyG(routeG,2)];
    ugv_time = ugv_distance/ugv_data.Vugv;
    
    %Search Algorithm to UAV Path
    [uav_path1, uav_path2, uav_distance, uav_time, stops] = uav_compute_path(xyG, routeG, home_x, home_y, x, y, wp_c, scp_table, SolL, vres, showResults, [xyG(routeG,1) xyG(routeG,2)], uav_data, G(i), saveResults, fullname);
    
    %i'st Solution
    total_time = ugv_time + uav_time;
    total_distance = ugv_distance + uav_distance;
    
    s = struct('f_ugv_d',ugv_distance,'f_ugv_t',ugv_time,'f_uav_d',uav_distance,'f_uav_t',uav_time,'ftotal_d',total_distance,'ftotal_t',total_time,'stops',stops);
    data_sol = [data_sol ; s];
    [ p_sol ] = build_matrix_solution( ugv_path, uav_path1, uav_path2 );
    path_sol = [path_sol ; struct('path_solution',p_sol)];

end

end
