% WELCOME TO THE TERRA2D ALGORITHM
%  This function executes the whole solution designed and developed to solve
%  the Energy Constrained UAV and Charging Station UGV Routing Problem (ECU-CSURP).
%
%  The goal is to visit all the waypoints with the UAV, supporting it with
%  UGV as a moving recharging station.

%  - INPUTS -
%  cfg_params     :configuration parameters
%  problem_params :ECU_CSURP Parameters
%  map_data       :xy coordinates
%  ugv_data       :ugv parameters
%  uav_data       :uav parameters

%  - OUTPUTS -
%  data_sol :list of structures with the cost of the paths of each solution.
%   struct {
%      f_ugv_d       :Cost of the UGV as distance travelled in meters.
%      f_ugv_t       :Cost of the UGV as time in seconds.
%      f_uav_d       :Cost of the UAV as distance travelled in meters. (SEARCH_UAV_PATH)
%      f_uav_t       :Cost of the UAV as time in seconds. (SEARCH_UAV_OPERATIONS)
%      ftotal_d      :Total cost in meters.
%   }
%  path_sol :paths of both UGV-UAV to reach the solutions computed in the algorithm.
%  figures  :list of figures with the solutions of each algorithm stage

function [data_sol, path_sol, figures] = TERRA(cfgParams, problem_params, uav_data, ugv_data)

%Init Output Variables
data_sol = [];
path_sol = [];
figures  = [];
figGo = [];

% (1) Compute Voronoi Diagram
%[ V1, wp_c, figV] = Voronoi_Covering_Optimization_nearVertex(problem_params, cfgParams);
[ V1, wp_c, figV] = Voronoi_Covering_TimeOptimization(problem_params, uav_data, cfgParams);

% (2) Compute Set Covering Problem
[SolL,scp_table,V2, figGr] = greedy_scp(problem_params, V1, wp_c, uav_data, cfgParams);

%Check If Home is a vertex of the solution
[Vsol] = checkHome(V2,problem_params.Home);

if (isempty(problem_params.Gp))
    % (3.A) UGV's Path without Gp
    [~,ugv_distance,ugv_path] = tsp_ga_ugv(Vsol, ugv_data.ugv_tsp);
    ugv_time = ugv_distance/ugv_data.Vugv;

    % (3.B) UAV's Path without Gp
    [uav_path1, uav_path2, uav_distance, uav_time, stops] = uav_compute_path(wp_c, scp_table, SolL, V1, cfgParams, ugv_path, uav_data);

    %1st Solution without GOA
    total_time = ugv_time + uav_time;
    total_distance = ugv_distance + uav_distance;
    s = struct('f_ugv_d',ugv_distance,'f_ugv_t',ugv_time,'f_uav_d',uav_distance,'f_uav_t',uav_time,'ftotal_d',total_distance,'ftotal_t',total_time,'stops',stops);
    data_sol = [data_sol ; s];
    [ p_sol ] = build_matrix_solution( ugv_path, uav_path1, uav_path2 );
    path_sol = [path_sol ; struct('path_solution',p_sol)];

else
    % (4) GOA Execution
    if (problem_params.Gp=='GravityCenter')
        avg_x = mean(Vsol(1,:));
        avg_y = mean(Vsol(2,:));
    elseif (problem_params.Gp=='HomeCenter')
        avg_x = problem_params.Home(1,1);
        avg_y = problem_params.Home(2,1);
    elseif (problem_params.Gp=='MedianCenter')
        avg_x = median(Vsol(1,:));
        avg_y = median(Vsol(2,:));
    end
    % (4.1) GOA 
    [vx_opt, vy_opt, figGo, vres] = Gravitational_Optimization( problem_params, V1, avg_x, avg_y, wp_c, SolL, scp_table, cfgParams);
    
    % (4.2) Genetic Algorithm to UGV Path
    [~,ugv_distance,ugv_path] = tsp_ga_ugv([[problem_params.Home(1,1) vx_opt];[problem_params.Home(2,1) vy_opt]], ugv_data.ugv_tsp);
    ugv_time = ugv_distance/ugv_data.Vugv;
    
    % (4.3) Search Algorithm to UAV Path
    [uav_path1, uav_path2, uav_distance, uav_time, stops] = uav_compute_path(wp_c, scp_table, SolL, vres, cfgParams, ugv_path, uav_data);
    
    %i'st Solution
    total_time = ugv_time + uav_time;
    total_distance = ugv_distance + uav_distance;
    s = struct('f_ugv_d',ugv_distance,'f_ugv_t',ugv_time,'f_uav_d',uav_distance,'f_uav_t',uav_time,'ftotal_d',total_distance,'ftotal_t',total_time,'stops',stops);
    data_sol = [data_sol ; s];
    [ p_sol ] = build_matrix_solution( ugv_path, uav_path1, uav_path2 );
    path_sol = [path_sol ; struct('path_solution',p_sol)];

end

figR = figure;
hold on
plot(problem_params.T(1,:),problem_params.T(2,:),'.blue');
plot(problem_params.Home(1,:),problem_params.Home(2,:),'*black')
plot(V2(1,:),V2(2,:),'+black');
plot(ugv_path(:,1),ugv_path(:,2),'-green');
[~,c]= size(V2);
theta = linspace(0,2*pi);
for i=1:c
    c_x(i,:) = problem_params.R*cos(theta) + V2(1,i);
    c_y(i,:) = problem_params.R*sin(theta) + V2(2,i);
    plot(c_x(i,:),c_y(i,:),'r:');
end
[~,c] = size(uav_path2);

for i=1:c
    plot(uav_path2(i).Coordinates(:,1),uav_path2(i).Coordinates(:,2),'-blue');
end

%plot(uav_path1(1,:),uav_path1(2,:),'-blue');

axis equal;
hold off;
if (isempty(problem_params.Gp))
    title('Final Solution without Gp');
else
    title(['Final Solution with Gp=' problem_params.Gp]);
end

%Save All figures
figures = [figures figV figGr figGo figR];

end

function [ugv_path] = checkHome(V2,Home)
    %Find out if home is part of the vertices solution, in that case, put into
    %the first position
    sol_x = V2(1,:);
    sol_y = V2(2,:);
    home_x = Home(1,:);
    home_y = Home(2,:);
    ugv_path = [];
    enc = false;
    for i=1:length(sol_x)
        if sol_x(i) == home_x && sol_y(i) == home_y
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
    ugv_path = [path_x;path_y];
end


