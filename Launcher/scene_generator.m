function [ new_problem ] = scene_generator(problem_params)
% This function creates a map ensuring that 'n' target points are distributed 
% in 'delta' groups inside an specific radius 2x'r', around a map 
% with an specific area. 
% Mandatory: 'n' must be multiple of 'delta'

% INPUTS
% problem_params

% OUTPUTS
% xy = coordinates of the target points

% Example:
% n = 12
% r = 3
% area = 50 km2
% delta = 2
% Solution: This algorithm will create 6 groups of 2 target points inside a
% r = 3 km of radius in a 50km2 area

new_problem = problem_params;
delta = problem_params.D;
area = problem_params.Area;
n = problem_params.N;
home = problem_params.Home;
r = problem_params.R;

groups = [];
xy = [];
theta = linspace(0,2*pi);
finish = false;

if (n>0)
    if (area>0)
        if (delta>0)
            n_group = n / delta;
            rng('shuffle','twister');
            while (delta > 0)
                %Random center of the group
                g = [(area*rand(1,1));(area*rand(1,1))];    
                a = (1.5)*rand(1,1)+0.5;
                sigma = a*r;
                
                % Target points following a normal distribution
                t = g + normrnd(0,sigma,[2,n_group]);
                xy = [xy  t];
                
                groups = [groups g];
                delta = delta - 1;

            end
        else
            disp('delta must be > 0')
        end
    else
       disp('area must be > 0') 
    end
else
   disp('n must be > 0') 
end

%Round to 3 decimals
f = 10.^3;
xy = round(f*xy)/f;

new_problem.T = xy;

if (~finish)
    %Draw
    figure;
    plot(xy(1,:),xy(2,:),'blue.',home(1,1),home(2,1),'*');
    hold on;
    [~,c] = size(groups);
    for i=1:c
        c_x(i,:) = r*sin(theta) + groups(1,i);
        c_y(i,:) = r*cos(theta) + groups(2,i);
        plot(groups(1,i),groups(2,i),'black+',c_x(i,:),c_y(i,:),'r:');
    end
    hold off;
    axis equal;
    title('Random map generated');
end
end



