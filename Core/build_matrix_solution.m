function [ path_sol ] = build_matrix_solution( ugv_path, uav_path1, uav_path2 )
%%% Build Final Solution Path for the UGV and UAV
path_sol = [];
for i=1:length(ugv_path)
    path_sol = [path_sol struct('c_ugv',ugv_path(i,:),'c_uav1',[],'c_uav2',[])];
end

%Insert c_uav1 to the final path solution
[~,r] = size(path_sol);
[~,u] = size(uav_path1);
for i=1:r
    for j=1:u
        if (path_sol(i).c_ugv(1)==uav_path1(j).Coordinates(1,1) && path_sol(i).c_ugv(2)==uav_path1(j).Coordinates(1,2))
            path_sol(i).c_uav1 = uav_path1(j).Coordinates;
        end
    end
end

%Insert c_uav2 to the final path solution
[~,u] = size(uav_path2);

for i=1:r
    for j=1:u
        if (path_sol(i).c_ugv(1)==uav_path2(j).Coordinates(1,1) && path_sol(i).c_ugv(2)==uav_path2(j).Coordinates(1,2))
            path_sol(i).c_uav2 = uav_path2(j).Coordinates;
        end
    end
end


end

