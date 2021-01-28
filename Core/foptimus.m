function [sx,sy] = foptimus(a, b, Xm, Ym, R)
 % This function solve equations system to among the WPs Cirfunference's
 % Equations and the line between the points (avg_x,avg_y) and (vx,vy)
 
 % Cirfunference's Equation
 % (x-a)^2 - (y-b)^2 = R^2
 
 % Line's Equation between two points
 % (x - Xm1)/(Xm2 - Xm1) = (y - Ym1)/(Ym2 - Ym1)
 
 % INPUTS
 % a - vector with the a's center of the WPs circunferences C(a,b)
 % b - vector with the b's center of the WPs circuenfreces C(a,b)
 % Xm - vector with the X coordinates of the points of the vertices
 % Ym - vector with the Y coordinates of the points of the vertices
 
 % OUTPUT
 % sx - Coordinate X of the optimal point of this set(vertice,WP1,WP2, ..., WPN)
 % sy - Coordinate Y of the optimal point of this set(vertice,WP1,WP2, ..., WPN)
 
 
 %Variables of the Equation's Systems
 syms x y
 vars = [x y];
 
 %Build the Line's equation between the vertice and the average
 %center
 eqL = ((x-Xm(1))/(Xm(2)-Xm(1)))==((y-Ym(1))/(Ym(2)-Ym(1)));
 px = [];
 py = [];
 
 for i=1:length(a)
     eqC = ((((x-a(i))^2) + ((y-b(i))^2)) == R^2);
     eqns = [eqC, eqL];
     [solx, soly] = solve(eqns,vars);
     
     for j=1:length(solx)
         px = [px double(vpa(solx(j)))];
     end
     for j=1:length(soly)
         py = [py double(vpa(soly(j)))];
     end
 end
 

 %For each solution, catch the ones who cover all the waypoints
 p_x = [];
 p_y = [];
 for j=1:length(px)
     wp_c_cnt = 0;
     
     for k=1:length(a)
         d = sqrt( ((px(j)-a(k))^2) + ((py(j)-b(k))^2) );
         if (vpa(d) <= R)
             wp_c_cnt = wp_c_cnt + 1;
         end
     end
     if (wp_c_cnt == length(a))
         p_x = [p_x px(j)];
         p_y = [p_y py(j)];
     end
 end
 
 %For those which cover all the waypoints, select the nearest to the
 %gravity center
 dmin = Inf;
 minX = Inf;
 minY = Inf;
 for j=1:length(p_x)
     
     d = sqrt( ((p_x(j)-Xm(2))^2) + ((p_y(j)-Ym(2))^2) );
     if (vpa(d) < dmin)
         dmin = d;
         minX = p_x(j);
         minY = p_y(j);
     end
 end
 
 sx = minX;
 sy = minY;
 
end
