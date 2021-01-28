function [ c_time, data_results, path_solution] = Test_2D(test_name, iterations, radius, N, A, D, vrep, xy)
% struct {
%   Vuav = 18                       :Speed of the UAV in km/h
%   Ct = 0                          :Charging time in h
%   Tl = 20/3600                    :Landing time in h
%   To = 10/3600                     :Taking off time in h
%   Tg = 5/3600                     :Time to execute the goal in h
%   R = [(Vuav * (Dr / Em)[s]) / 2] :Radius maximus in km
% }
uav_data = struct('Vuav',18,'Ct',0,'Tl',20/3600,'To',10/3600,'Tg',5/3600,'R',radius);
ugv_data = struct('Vugv',0.13); %Mars Curiosity Max Speed in km/h

%General Variables
it_completed = 0;
showResults = false;
saveResults = false; %Save Results only if showResults is True
saveData = false;
fullname = '';
c_time = [];
data_results = [];
path_solution = [];
%TERRA Parameters
home = [0.5;0.5];
last_depot = 1;% 1 = HOME; IT NEEDS TO BE CHOOSEN BEFORE THE SCP SOLUTION!!!!

for i=1:iterations
    tic
    %[xy] = map_generator(home,N,uav_data.R,A,D,1000);
    [xy] = map_generatorv2(home,N,uav_data.R,A,D);
    [~,c] = size(xy);
    if (c==N)
        [data_sol, path_sol] = TERRA(showResults, home, xy, last_depot, uav_data, ugv_data, saveResults, fullname);
        c_time = [c_time toc];
        data_results = [data_results; struct('It',data_sol)];
        path_solution = [path_solution; struct('It',path_sol)];
        sprintf('COMPLETED - Iteration %d',i)
        it_completed = it_completed + 1;
    end
end

if (saveData)
   content = 'D:\OneDrive\Universidad\Doctorado\Matlab\TERRA\Results'; 
   %PC-HOME = C:\Users\Fernando\OneDrive\Universidad\Doctorado\Matlab\TERRA\Results
   %LAPTOP = 'D:\OneDrive\Universidad\Doctorado\Matlab\TERRA\Results'
   fullname = [content,'\',test_name];
   mkdir(fullname);
   filename = [fullname,'\','data.mat'];
   save(filename);
end

fprintf('SUCCESSFULLY FINISHED Test :%s\n',test_name)

if (vrep)
    %MainLoop(path_sol(1).path_solution,toc,false);

end


