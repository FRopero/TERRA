%   Syntax:
%   TSPsolution = LKH_TSP(CostMatrix,pars_struct,fname_tsp,LKHdir,TSPLIBdir)
%
%   This functions solves TSP problems using the Lin-Kernighan-Helsgaun
%   solver. It assumes that a compiled executable of the LKH solver as
%   found at: http://www.akira.ruc.dk/~keld/research/LKH/ is available at
%   the LKHdir directory. Furthermore a TSPLIB directory is assumed.
%   For the definition of the TSPLIB and the compilation of the LKH code
%   check the aforementioned url. 
%
%   Inputs:
%   CostMatrix      : the Cost Matrix of the (asymmetric) TSP. [e.g. it can be an NxN matrix of distances]
%   pars_struct     : parameters structure with
%                   : -> CostMatrixMulFactor (value that makes Cost Matrix
%                        almost integer. [eg. pars_struct.CostMatrixMulFactor = 1000; ]
%                     -> user_comment (a user comment for the problem) [optional]
%   fname_tsp       : the filename to save the tsp problem
%   LKHdir          : the directory of the LKH executable
%   TSPLIBdir       : the directory of the TSPLIB files
%   
%   Outputs:
%   TSPsolution     : the TSP solution
%   
%   Authors:
%   Kostas Alexis (kalexis@unr.edu)
%

function [TSPsolution,TSPcost] = LKH_TSP(CostMatrix,pars_struct,fname_tsp,LKHdir,TSPLIBdir)

CostMatrix_tsp = pars_struct.CostMatrixMulFactor*CostMatrix;
CostMatrix_tsp = floor(CostMatrix_tsp);
user_comment = pars_struct.user_comment; 

fileID = writeTSPLIBfile_FE(fname_tsp,CostMatrix_tsp,TSPLIBdir,user_comment);

%disp('### LKH problem set-up...');
%%%  Solve the TSP Problem via the LKH Heuristic
%disp('### Now solving the TSP problem using the LKH Heuristic...');
%start_lkh_time = cputime;

%lkh_cmd = [LKHdir 'lkh-3.0.5' ' ' TSPLIBdir fname_tsp '.par'];
lkh_cmd = ['"' LKHdir 'lkh-3.0.5' '"' ' ' '"' TSPLIBdir fname_tsp '.par' '"'];
[~,~] = system(lkh_cmd,'');
%end_lkh_time = cputime;
%disp('### ... done!');

%movefile([fname_tsp '.txt'],TSPLIBdir)
%solution_file = [TSPLIBdir fname_tsp '.txt']; 
[TSPsolution,TSPcost] = readSolution([fname_tsp '.txt']);

end

function [TSPsol,TSPcost] = readSolution(file)
TSPsol = [];
TSPcost = 0;
fid = fopen(file, 'r');
tline = fgetl(fid);
while ischar(tline)
    if (contains(tline,'NAME'))
        C = strsplit(tline,'.');
        TSPcost = str2double(C{2});
    elseif(contains(tline,'TOUR_SECTION'))
        tline = fgetl(fid);
        while ischar(tline)
            TSPsol = [TSPsol str2double(tline)];
            tline = fgetl(fid);
        end
        break;
    end
    tline = fgetl(fid);
end
fclose(fid);

end

function [fileID] = writeTSPLIBfile_FE(fname_tsp,CostMatrix,tsplib_dir,user_comment)
    dims_tsp = length(CostMatrix);
	name_line = ['NAME : '  fname_tsp  '\n'];
	type_line = ['TYPE: TSP' '\n'];
	comment_line = ['COMMENT : '  user_comment  '\n'];
	tsp_line = ['TYPE : TSP \n'];
	dimension_line = ['DIMENSION : '  int2str(dims_tsp)  '\n'];
	edge_weight_type_line = 'EDGE_WEIGHT_TYPE : EXPLICIT\n'; %# explicit only
	edge_weight_format_line = 'EDGE_WEIGHT_FORMAT: FULL_MATRIX\n';
%	display_data_type_line ='DISPLAY_DATA_TYPE: ' + 'NO_DISPLAY' + '\n' %# 'NO_DISPLAY'
	edge_weight_section_line = 'EDGE_WEIGHT_SECTION\n';
	eof_line = 'EOF\n';
	Cost_Matrix_STRline = {};
	for i=1:dims_tsp
		cost_matrix_strline = '';
		for j=1:dims_tsp-1
			cost_matrix_strline = {cost_matrix_strline ; {int2str(int32(CostMatrix(i,j)))  ' '}};
        end
		j = dims_tsp-1;
		cost_matrix_strline = {cost_matrix_strline ; int2str(int32(CostMatrix(i,j)))};
		cost_matrix_strline = {cost_matrix_strline  '\n'};
        if (isempty(Cost_Matrix_STRline))
            Cost_Matrix_STRline = cost_matrix_strline;
        else
            Cost_Matrix_STRline = {Cost_Matrix_STRline ; cost_matrix_strline};
        end
    end   
    file = [tsplib_dir fname_tsp '.tsp'];
    %edit(file);
	fileID = fopen(file, 'wt');
	%print name_line;
    fprintf(fileID,name_line);
	fprintf(fileID,comment_line);
	fprintf(fileID,tsp_line);
	fprintf(fileID,dimension_line);
	fprintf(fileID,edge_weight_type_line);
	fprintf(fileID,edge_weight_format_line);
	fprintf(fileID,edge_weight_section_line);
    
	for i=1:length(CostMatrix) 
        for j=1:length(CostMatrix)
            fprintf(fileID,'%g ',CostMatrix(i,j));
        end
        fprintf(fileID,'\n');
    end
	fprintf(fileID,eof_line);
	fclose(fileID);

    file2 = [tsplib_dir fname_tsp '.par'];
	fileID2 = fopen(file2, 'wt');
	problem_file_line = ['PROBLEM_FILE = ' tsplib_dir fname_tsp '.tsp\n']; %# remove pwd + tsplib_dir
	optimum_line = 'OPTIMUM = 378032\n';
	move_type_line = 'MOVE_TYPE = 5\n';
	patching_c_line = 'PATCHING_C = 3\n';
	patching_a_line = 'PATCHING_A = 2\n';
	runs_line = 'RUNS = 1\n';
	tour_file_line = ['TOUR_FILE = ' fname_tsp '.txt\n'];

	fprintf(fileID2,problem_file_line);
	fprintf(fileID2,optimum_line);
	fprintf(fileID2,move_type_line);
	fprintf(fileID2,patching_c_line);
	fprintf(fileID2,patching_a_line);
	fprintf(fileID2,runs_line);
	fprintf(fileID2,tour_file_line);
	fclose(fileID2);
end