% function [solution] = Test_2D(problem_params, dana_params, cfgParams)
%  This function tests TERRA2D for a particular ECU-CSURP scenario
%  Inputs:
%    problem_params - ECU_CSURP Parameters
%    uav_data       - data info of the uav
%    ugv_data       - data info of the ugv
%    cfgParams      - Configuration Parameters of the Test
%  Outputs:
%    solution       - struct with the ouput parameters      
function [solution] = Test_2D(problem_params, uav_data, ugv_data, cfgParams)
    cnt = 1;
    solution = [];
    
        for i=1:cfgParams.iterations
            %Create save dir for the results of this iteration
            itdir = [cfgParams.saveDir 'It-' num2str(cnt) cfgParams.slash];
            mkdir (itdir);

            % Generate Random 2D Scenario
            [problem_params] = scene_generator(problem_params);
            % Yu. et al (2018) Random Generator
            %[problem_params] = graphMakingNew([300 0 110 0], 1, problem_params);
            disp(['Computing solution for Scenario ' num2str(cnt) ' ...']);
            tic
            [data_sol, path_sol, figures] = TERRA(cfgParams, problem_params, uav_data, ugv_data);
            tc = toc;
            
            %Saving
            if (cfgParams.saveResults)
               [solution] = saveResults(cfgParams,tc,data_sol,path_sol,problem_params,uav_data,ugv_data,figures,itdir);
            end

            disp(['Execution ' num2str(cnt) ' completed.']);
            cnt = cnt + 1;
        end

    %Simulation Stage
    if (cfgParams.Vrep)
        MainLoop(path_sol(1).path_solution);
    end

end


%This function stacks the iteration output to the general output, and also, save
%the iteration output
function [solution] = saveResults(cfgParams,tc,data_sol, path_sol,problem_params,uav_data,ugv_data,figures, itdir)
    disp('Saving results ...');
    
    % Init solution struct
    %  cT      - array of the computational time taken by TERRA to compute a solution
    %  dataRes - array of the computational results of the TERRA execution (see definition in TERRA3D.m)
    %  pathRes - array of the UGV-UAV path planning solution of TERRA
    solution = struct('cT',[],'dataRes',[],'pathRes',[]);
    
    %Stacking
    solution.cT = tc;
    solution.dataRes =struct('It',data_sol);
    solution.pathRes = struct('It',path_sol);
    
    %Saving Data
    fpath = [itdir 'data.mat'];
    save(fpath,'tc','data_sol','path_sol','problem_params','uav_data','ugv_data','cfgParams');
    
    %Saving Figures
    if (~isempty(figures))
        [~,c] = size(figures);
        for i=1:c
            name = [itdir figures(i).Name '.fig'];
            figure(figures(i));
            savefig(name);
        end
    end
    
    if (~cfgParams.printResults) 
        close all;
    end
   
end


