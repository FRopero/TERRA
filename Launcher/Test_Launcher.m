%                        -> 2D TERRA Launcher <-                          %

% This script launchs an execution of TERRA-2D over a random generated
% scenario to solve the Energy Constrained UAV and Charging Station UGV 
% Routing Problem (ECU-CSURP).
% For more information, please refer to the paper:
%   Ropero, F., Muñoz, P., & R-Moreno, M. D. TERRA: A path planning 
%   algorithm for cooperative UGV–UAV exploration. Engineering Applications 
%   of Artificial Intelligence, 78, 260-272, 2019.

% System Settings
if isunix
    slash='/';
elseif ispc
    slash='\';
end

% Check execution directory
if ~contains(pwd,'TERRA')
    disp('Execution error: You need to be in the "TERRA\" directory')
    return
end

% Configuration of the testing parameters
%  iterations   = nº of executions with different scenarios
%  printResults = LoL
%  saveResults  = xd
%  Vrep         = to launch the V-REP simulation
%  saveDir      = o.O
%  LKHdir       = directory where it is installed the TSP heuristic (LKH). Mandatory bars format 'C:\\Users\\'
format shortg
[c, ~] = clock;
saveDir = ['Results' slash 'Test-' num2str(c(3)) '.' num2str(c(2)) '.' num2str(c(1)) '.' num2str(c(4)) '.' num2str(c(5)) slash];
LKHdir = 'C:\\Users\\ferna\\OneDrive - Universidad de Alcala\\Universidad\\Doctorado\\Matlab\\TERRA\\EstimatingHeuristic\\';
cfgParams = struct('iterations',1,'printResults',true,'saveResults',true,'Vrep',false,'saveDir',saveDir,'slash',slash,'fullname','','LKHdir',LKHdir);

%Create directory if saveresults = true
if (cfgParams.saveResults)
    [status, msg, msgID] = mkdir(saveDir);
end

% ECU_CSURP Parameters
%  T    = Target Points of the scenario
%  N    = Number of target points
%  R    = Farthest distance the UAV can travel [m]
%  D    = Number of groups for the random map generator
%  Home = Home location
%  Area = map area
%  Gp   = gravitational point (if null, then no GOA is applied) - Gps = 'GravityCenter','HomeCenter','MedianCenter'
problem_params = struct('T',[],'R',0,'N',45,'D',1,'Home',[0.5;0.5],'Area',100,'Gp', []);

% UGV's TSP Genetic Algoritm Parameters (Current Used = Id. 2 from TERRA paper)
%  popSize     = size of the population of chromosomes.
%  tournaments = number of chromosomes to select for the tournament selection
%  mutOper     = mutation operator (1 = Flip; 2 = Swap; 3 = Slide)
%  mutRate     = mutation rate from 0.0 to 1.0 (0-100%)
%  crossOper   = crossover operator (1 = OX; 2 = CX; 3 = CBX)
%  eliteP      = elitism selection from 0 to 100%
ugv_tsp = struct('popSize',430,'tournaments',9,'mutOper',2,'mutRate',0.06,'crossOper',1,'eliteP',2.7);

% UGV Path Planning Parameters
ugv_data = struct('Vugv',0.4, 'ugv_tsp',ugv_tsp); %Mars Curiosity Max Speed in m/s

% UAV Path Planning Parameters
%   Tt   :Total time budget [budget] budget = flight seconds
%   Vuav :UAV's speed [m/s]
%   Tl   :Landing time [s]
%   To   :Taking off time [s]
%   R    :Radius [m]
uav_data = struct('Tt',308,'Vuav',0.5,'Tl',4,'To',4,'R',0);
problem_params.R = uav_data.Vuav * ((uav_data.Tt - uav_data.To - uav_data.Tl)/2); %m
uav_data.R = problem_params.R;

%Launch TERRA-2D
[solution] = Test_2D(problem_params, uav_data, ugv_data, cfgParams);


