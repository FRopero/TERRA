function [ vx_opt, vy_opt, figG, vres] = Gravitational_Optimization( problem_params, V1, avg_x, avg_y, wp_c, SolL, scp_table, cfgParams )
%This algorithm optimizes the UGV path for the FCURP-MRS
figG = [];
vres = zeros(2,length(V1(1,:)));

wp_c_aux = wp_c;
[r,c] = size(scp_table);
vx_opt = [];
vy_opt = [];
%theta = linspace(0,2*pi);

for i=1:c
    for k=1:length(SolL)
        a = [];
        b = [];
        if (i==SolL(k) && V1(1,i)~=problem_params.Home(1,:) && V1(2,i)~=problem_params.Home(2,:))
            Xm = [V1(1,i) avg_x]; %Vertice of the initial solution
            Ym = [V1(2,i) avg_y];
            if (V1(1,i)~=avg_x && V1(2,i)~=avg_y)
                for j=1:r
                    if (scp_table(j,i)==1 && wp_c_aux(1,j) ~= -1) %Waypoint j
                        a = [a wp_c(1,j)];
                        b = [b wp_c(2,j)];
                        %c_x(i,:) = problem_params.R*cos(theta) + wp_c(1,j);
                        %c_y(i,:) = problem_params.R*sin(theta) + wp_c(2,j);
                        wp_c_aux(1,j) = -1;
                    end
                end
                [sx,sy] = foptimus(a, b, Xm, Ym, problem_params.R);
            else
                sx = V1(1,i);
                sy = V1(2,i);
            end
            vx_opt = [vx_opt sx];
            vy_opt = [vy_opt sy];
            
            vres(:,i) = [sx;sy];
            
           elseif (i==SolL(k) && V1(1,i)==problem_params.Home(1,1) && V1(2,i)==problem_params.Home(2,1))
             vres(:,i) = [V1(1,i);V1(2,i)];
        end
    end
end

end

