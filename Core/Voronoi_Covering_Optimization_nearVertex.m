function [ V1, wp_c, figV] = Voronoi_Covering_Optimization_nearVertex( problem_params, cfgParams)
%VORONOI_COVERING_OPTIMIZATION Summary of this function goes here
%   This function uses the voronoi diagram to optimize the problem of
%   covering all the waypoints received in the problem.
figV = [];

%Algorithm variables
V1 = [];
wp_c = []; 
vx_p = [];
vy_p = [];
duplicates = 0;
near_wd = [];
wp_c_cnt = 0; 
wp_nc = problem_params.T;
wp_nc_cnt = length(wp_nc(1,:));
art_vert = problem_params.Home;
R = problem_params.R;
v_used = [];

% if (cfgParams.printResults)
%     vis = 'on';
% else
%     vis = 'off';
% end
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


while ( 0 < wp_nc_cnt  )
    
    % [STEP 1] - Calculate the Voronoi's Diagram over the WPm not covered
    
    if (~isempty(vx_p)) %For n iterations
        
        %Search and catch the nearest vertice to each wp_nc and add to
        %near_v
        [~,c_wp] = size(wp_nc);
        [~,c_vx] = size(vx_p);
        near_v = [];
        t = 1;
        dminT = Inf;
        for i=1:c_wp
            dmin = sqrt( ((wp_nc(1,i)-vx_p(1))^2) + ((wp_nc(2,i)-vy_p(1))^2) );
            near_v(1,i) = vx_p(1);
            near_v(2,i) = vy_p(1);
            for j=1:c_vx
                d = sqrt( ((wp_nc(1,i)-vx_p(j))^2) + ((wp_nc(2,i)-vy_p(j))^2) );
                if (d < dmin)
                    dmin = d;
                    near_v(1,t) = vx_p(j);
                    near_v(2,t) = vy_p(j);
                    t = t + 1;
                end
            end
            if (dmin < dminT)
                dminT = dmin;
            end
        end
        
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
        
        [~,lv] = size(near_wd);
        [~,lu] = size(v_used);
        
        %Detect duplicates and so, infinite loop in voronoi diagram
        duplicates = 0;
        for i=1:lv
            for j=1:lu
                if (near_wd(:,i)==v_used(:,j))
                    duplicates = duplicates + 1;
                end
            end
        end
        v_used = [v_used near_wd];
    end

    %Add these vertices to the next voronoi iteration
    wp_plus_v = [wp_nc near_wd];
    
    %We can't do voronoi if we have 2 points (one wp and its nearest vertice)
    %or Only 2 wp and 1 vértice, and the distance between the WPs and their common nearest
    %vertice, is less than the R. This situation cause an infinite loop.
    if ((length(wp_plus_v(1,:))==2)||(duplicates > 0))
        %m is the intermediate point between the WPs. Each iteration,
        %this vertice is nearest to the first wp
        m_x = (wp_plus_v(1,1) + wp_plus_v(1,2)) / 2;
        m_y = (wp_plus_v(2,1) + wp_plus_v(2,2)) / 2;
        
        d = sqrt( ((wp_plus_v(1,1)-m_x)^2) + ((wp_plus_v(2,1)-m_y)^2) );
        
        while (d > R)
            m_x = (m_x + wp_plus_v(1,1)) / 2;
            m_y = (m_y + wp_plus_v(2,1)) / 2;
            d = sqrt( ((wp_plus_v(1,1)-m_x)^2) + ((wp_plus_v(2,1)-m_y)^2) );
        end
        art_vert = [art_vert [m_x;m_y]];
        
    else
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
    end
    
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
            
            %fprintf('Waypoint[%.4f,%.4f] - Vértice[%.4f,%.4f]',wp_nc(1,wp_v),wp_nc(2,wp_v),vx_p(best_v),vy_p(best_v));
            %fprintf('\n');
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

V1 = [vx_sol; vy_sol];

%Last iteration plot
% figLast = figure('visible',vis);
% figLast.Name = 'Voronoi_Solution';
% if (cfgParams.saveResults)
%     figV = [figV figLast];
% end
% hold on
% plot(problem_params.T(1,:),problem_params.T(2,:),'blue.',vx_sol,vy_sol,'green+');
% [~,c]= size(vx_sol);
% theta = linspace(0,2*pi);
% for i=1:c
%     c_x(i,:) = problem_params.R*cos(theta) + vx_sol(i);
%     c_y(i,:) = problem_params.R*sin(theta) + vy_sol(i);
%     plot(c_x(i,:),c_y(i,:),'r:');
% end
% hold off
% axis equal
% title('Voronoi Solution');

end

