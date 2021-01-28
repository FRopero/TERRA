function [ vx_sol, vy_sol, wp_c ] = Voronoi_Covering_Optimization_allVertex( x,y,home_x,home_y, R, theta, fig1, showResults )
%VORONOI_COVERING_OPTIMIZATION Summary of this function goes here
%   This function uses the voronoi diagram to optimize the problem of
%   covering all the waypoints received in the problem.

%Algorithm variables
wp_c = []; 
vx_p = [];
vy_p = [];
duplicates = 0;
near_wd = [];
wp_c_cnt = 0; 
wp_nc = [x;y];
wp_nc_cnt = length(wp_nc(1,:));
art_vert = [home_x;home_y];
v_used = [];

while ( 0 < wp_nc_cnt  )
    
    % [STEP 1] - Calculate the Voronoi's Diagram over the WPm not covered
    
    if (~isempty(vx_p)) %For n iterations
        near_v = [vx_p;vy_p];
        %Delete duplicates in near_v and add them to near_wd
        near_wd = [0;0];
        [~,c_nearv] = size(near_v);
        t = 1;
        for i=1:c_nearv
            val = near_v(:,i);
            in = 0;
            j = 1;
            while (j <= length(near_wd(1,:)))
                if (val==near_wd(:,j))
                    in = in + 1;
                end
                j = j + 1;
            end
            if (0==in)
                near_wd(:,t) = val;
                t = t + 1;
            end
        end
        
    end

    %Add these vertices to the next voronoi iteration
    wp_plus_v = [wp_nc near_wd];
    
    %Create the voronoi diagram for the n iteration
    [vx,vy] = voronoi(wp_plus_v(1,:),wp_plus_v(2,:));
        
    %Delete duplicates vertices
    p = 1;
    vertices = [vx(1,:); vy(1,:)];
    vertices_s = [0;0];
    for i=1:length(vertices)
        v = vertices(:,i);
        in = 0;
        [~,c] = size(vertices_s);
        for j=1:c
            if (v==vertices_s(:,j))
                in = in + 1;
            end
        end
        if in==0
            vertices_s(:,p) = v;
            p = p + 1;
        end
    end

    vx_p = vertices_s(1,:);
    vy_p = vertices_s(2,:);
    
    %Add artificial vertices created in the last iteration
    vx_p  = [art_vert(1,:) vx_p];
    vy_p = [art_vert(2,:) vy_p];
    
    %Compute the new set of uncovered waypoints
    [wx,wy] = size(wp_nc);
    for i=1:wy
        dmin = R;
        best_v = -1;
        wp_v = -1;
        for j=1:length(vx_p)
            if (wp_nc(1,i)~=-1)
                d = sqrt( ((wp_nc(1,i)-vx_p(j))^2) + ((wp_nc(2,i)-vy_p(j))^2) );
                if (d < dmin && best_v ~= 1) %%Si el home ya lo ha cogido, es la mejor opción
                    dmin = d;
                    best_v = j;
                    wp_v = i;
                end
            end
        end
        %Actualizar tablas con el vértice de la mínima distancia al wp_nc
        % The Waypoint wp_nc(x,i) was covered by unless one vertice!!
        if (best_v ~= -1)
            wp_c_cnt = wp_c_cnt + 1;
            wp_c(1,wp_c_cnt) = wp_nc(1,wp_v);
            wp_c(2,wp_c_cnt) = wp_nc(2,wp_v);
            
            vx_sol(1,wp_c_cnt) = vx_p(best_v);
            vy_sol(1,wp_c_cnt) = vy_p(best_v);
            
            fprintf('Waypoint[%.4f,%.4f] - Vértice[%.4f,%.4f]',wp_nc(1,wp_v),wp_nc(2,wp_v),vx_p(best_v),vy_p(best_v));
            fprintf('\n');
            wp_nc(1,wp_v) = -1;
            wp_nc(2,wp_v) = -1;
            wp_nc_cnt = wp_nc_cnt - 1;
        end
    end
    
    del = 0; %To count the number of wp covered this iteration
    for i=1:wy
        if (wp_nc(1,i)==-1)
            del = del + 1;
        end
    end
    
    %Update wp_nc Table for the next iteration
    if (del>0)
        tmp_wp_nc = wp_nc;
        wp_nc = zeros(wx,wy-del);
        for i=1:wx
            t_y = 0;
            for j=1:wy
                if (tmp_wp_nc(i,j)~=-1)
                    t_y = t_y + 1;
                    wp_nc(i,t_y) = tmp_wp_nc(i,j);
                end
            end
        end
    end
    
end

% START Drawing
if (showResults)
    figure(fig1);
    subplot(2,2,2);
    [vx_t,vy_t] = voronoi(x,y);
    plot(x,y,'blue.',vx_sol,vy_sol,'black+',vx_t,vy_t);
    hold on;
    plot(home_x,home_y,'black*');

    for i=1:length(vx_sol)
        c_x(i,:) = R*cos(theta) + vx_sol(i);
        c_y(i,:) = R*sin(theta) + vy_sol(i);
        plot(c_x(i,:),c_y(i,:),'r:');
    end
    title('[STEP 1] - Voronoi Diagram Stage');
    hold off;
    axis equal;
end
% END Drawing

end

