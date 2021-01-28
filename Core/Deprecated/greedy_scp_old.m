function [solution_setsCell, solution_setsLabelsV, scp_table, sol_x, sol_y] = greedy_scp(x,y,home_x,home_y, vx_sol, vy_sol, wp_c, R, fig1, showResults, saveResults, fullname)
% MODIFIED VERSION [SolC,SolL] = GREEDYSCP(vx_sol, vy_sol, wp_c): where
% vx_sol is the 'x' coordinate of the vertexes of the solution and vy is the
% 'y' coordinate of the vertexes of the solution, and wp_c the list of all
% the covered waypoints. 
% The OUTPUT is the same of the previous version (see below).

% GREEDYSCP Greedy SCP algorithm:: function [solution_setsCell, solution_setsLabelsV] = greedyscp(setsCell, setsLabelsV) .
%	[SolC,SolL] = GREEDYSCP(C, L) if C is an array, creates a cell array SolC that is a solution of Set Cover Problem defined by C, where C{i} = S_i, an input set made by some of the elements we want to cover; 
%    SolC is made by the cells of C selected by the algorithm. The elements that we want to cover are indicates by numbers from 1 to n, where n is the number of elements we want to cover; 
%    therefore, C{i} is a vector of integers between 1 and n.
%
%	If C is a logical or numerical array of n rows, where C(j,i) > 0 iff element j is contained in set S_i, the output SolC will be a logical array made by the column of log(C) corresponding to the solution
%	
%	If a vector L of integer labels of the elements of C is provided, SolL contains the labels corresponding to SolC. Otherwise SolL contains the positions of elements of SolC in C. SolC and SolL elements are sorted in ascending order of SolL.
%
%	This is an implementation of the well-known greedy algorithm (Chvátal, 1979), with two small modifications:
%	* In case of more than one possible choice at one step, the biggest set is chosen.
%	* Once the solution is found, we check the selected sets to find a better cover solution, removing a set if is a subset of the union of the other set.
%	
%	If you use this code, please cite:
%	F. Gori, G. Folino, M.S.M. Jetten, E. Marchiori
%	"MTR: Taxonomic annotation of short metagenomic reads using clustering at multiple taxonomic ranks", Bioinformatics 2010.
%	doi = 10.1093/bioinformatics/btq649 
%

theta = linspace(0,2*pi);
%*=*=*=*=*
%* Input adjustment
%Get ready the SCP table without empty rows and empty columns and launch the greedy algorithm.
[wx,wy] = size(wp_c);
scp_table = zeros(wx,length(vx_sol));
%Creating NEW Vertices table with the wp_c
for i=1:wy
    for j=1:length(vx_sol)
        d = sqrt( ((wp_c(1,i)-vx_sol(j))^2) + ((wp_c(2,i)-vy_sol(j))^2) );
        if (d < R)
            scp_table(i,j) = 1;%j : to indicate the WP covered
        else
            scp_table(i,j) = 0;
        end
    end
end

%Building scp_table, without empty columns. Deleting the vertices which
%don't surround any waypoint
t = 1;
D = []; %Table containing the deleted vertices, because don't cover any wp
z = 1;
[~,c] = size(scp_table);
for i=1:c
    if (0<sum(scp_table(:,i)))
        scp_table_s(:,t) = scp_table(:,i);
        t = t + 1;
    else
        D(z) = i;
        z = z + 1;
    end
end

%Generate an adecuate labeling of the vertices, if it has been some vertice
%deleted
cnt = 0;
L = 1:c;
for i=1:c
    L_p(i) = L(i);
    for k=1:length(D)
        if D(k)==L(i)
            L_p(i) = 0;
            cnt = cnt + 1;
        end
    end
end

L_p = sort(L_p);
t = 1;
for i=(cnt+1):c
    L_s(t) = L_p(i);
    t = t + 1;
end

setsCell = scp_table_s;
setsLabelsV = L_s;

% Checking first input
cellL = iscell(setsCell) ;
if cellL
	if any(cellfun(@isempty, setsCell))
		disp('Error - Input array has some empty cells')
		return;
	end
	%Put cells of setsCell as rows vector
	setsCell = cellfun(@(x) x(:)', setsCell, 'UniformOutput', false) ;
	%setsCell must be a column
	setsN = numel(setsCell) ;
	if size(setsCell) ~= [setsN 1]
		setsCell = setsCell' ;
	end
	OriginalElementsV = unique([setsCell{:}]) ;
	elementsN = numel(OriginalElementsV) ;
	if (min(OriginalElementsV) == 1) && (max(OriginalElementsV) ==  elementsN)
		%to check that OriginalElementsV = 1 : elementsN, i just need to check that it has minimum 1 and numel(OriginalElementsV) = max(OriginalElementsV)
		setsCardinalitiesV = cellfun(@numel, setsCell) ;
	else
		disp('Error - some elements are missing')
		return;
	end
else
	A = setsCell ;
	if isnumeric(A)
		%disp('Warning : the input is a numeric array; trying to convert it to a logical array')
		A = logical(A) ;
	end	
	if islogical(A)
		if ~all(sum(A, 2))
			disp('Error - Input array has some empty rows!') 
			return;
		end
		[elementsN, setsN] = size(A) ;
		setsCardinalitiesV = sum(A, 1) ;
		if ~all(setsCardinalitiesV)
			disp('Error - Input array has some empty columns!') 
			return;
		end
	else
		disp('Error - Input array not of the type "cell", "logical" or "numeric"!') 
		return;
	end
end

% setsLabelsV must be a column
if nargin == 1
	%disp('Warning - Creating column labels'); 
	setsLabelsV = (1 : setsN)' ;
else %put in column
	if size(setsLabelsV) ~= [setsN 1]
		setsLabelsV = setsLabelsV' ;
	end
end

%*=*=*=*=*
%*Run SCP

if cellL
	[solution_setsCell, solution_setsLabelsV] = cellSCP(setsCell, setsLabelsV, setsCardinalitiesV, OriginalElementsV, elementsN, setsN) ;
else
	[solution_setsCell, solution_setsLabelsV] = matrixSCP(A, setsLabelsV, setsCardinalitiesV) ;
end

% START Drawing
if (showResults)
    figure(fig1);
    subplot(2,2,3)
    plot(x,y,'blue.');
    title('[STEP 2] - Set Cover Problem Stage');
    hold on;
    plot(home_x,home_y,'black*');
end
t = 1;
[~,c] = size(vx_sol);
for i=1:c
    for k=1:length(solution_setsLabelsV)
        if (i==solution_setsLabelsV(k))
            c_x(i,:) = R*cos(theta) + vx_sol(i);
            c_y(i,:) = R*sin(theta) + vy_sol(i);
            sol_x(t) = vx_sol(i);
            sol_y(t) = vy_sol(i);
            t = t + 1;
            if (showResults)
                plot(c_x(i,:),c_y(i,:),'r:',vx_sol(i),vy_sol(i),'black+');
            end
        end
    end
end
if (showResults)
    axis equal;
    hold off;
    if (saveResults)
        dir = strcat(fullname,'/','main.fig');
        saveas(gcf,dir,'fig');
    end
end
% END Drawing 

%end%end%end
%end%end%end

function [solution_setsCell, solution_setsLabelsV] = cellSCP(setsCell, setsLabelsV, setsCardinalitiesV, OriginalElementsV, elementsN, setsN)
%CELL Version

%*=*=*=*=*
%*Identify elements covered by just one set and select them

% numCoverVector(i) =  number of sets covering element i
numCoverVector = zeros(1, elementsN) ;
% posCoverVector(i) =  "last" set covering element i
posCovSetVector = zeros(1, elementsN) ;
for iSet = 1 : setsN 
    elementsCovLogical = ismember(OriginalElementsV, setsCell{iSet}) ;
    numCoverVector = numCoverVector + elementsCovLogical ;
	%Keep tack of the last set that cover each element, so that I'll know the set that uniquely covers a certain element
    posCovSetVector(elementsCovLogical) = iSet ;
end

%Select the sets that cover the elements just uniquely
UnCovLogical = numCoverVector == 1 ;
if any(UnCovLogical)
    toBeSelectVector = unique(posCovSetVector(UnCovLogical)) ;
    %display([num2str(sum(UnCovLogical)) ' elements are uniquely covered'])
    %display([num2str(numel(toBeSelectVector)) ' sets would be selected'])
    toBeSelectVector = unique(posCovSetVector(UnCovLogical)) ;
    solution_setsCell = setsCell(toBeSelectVector) ;
    solution_setsLabelsV = setsLabelsV(toBeSelectVector) ;  
    % update
    setsCell(toBeSelectVector) = [ ] ;
    setsLabelsV(toBeSelectVector) = [ ] ;
    setsCardinalitiesV(toBeSelectVector) = [ ] ;
    %
	elementsCoveredV = unique([solution_setsCell{:}]) ;
	%display(['Elements Covered: ' num2str(numel(elementsCoveredV))])
    % Equivalent: remainingElementsVector = setdiff(OriginalElementsV, elementsCoveredV) ;
    remainingElementsVector = OriginalElementsV(~ismember(OriginalElementsV, elementsCoveredV)) ;
end

%*=*=*=*=* 
%* Run greedy algorithm

if isempty(remainingElementsVector)
	%disp('No algorithm needed')
else
	%disp(' ')
	%Sort Sets in descending order according to their cardinality
	[setsCardinalitiesV, sortedIndexVector] = sort(setsCardinalitiesV, 'descend') ;
	setsLabelsV = setsLabelsV(sortedIndexVector) ;
	setsCell = setsCell(sortedIndexVector) ;
    %disp('Greedy algorithm: start')
	countN = 0 ;
    while ~isempty(remainingElementsVector) 
		%* Trick to speed up: if a set (the first, for instance) contains n unconvered elements, 
		% all the sets having cardinality < n will contain less uncovered elements than this set, so that I do not need to compute intVector for them because they can't be better that the set selected at the beginning. 
		% Since the sets are sorted according to their cardinality, I can just focus on the ones from 1 to focusedN.
        thresholdN = sum(ismember(setsCell{1}, remainingElementsVector)) ;
        
        indexFocusedSetsLogical = setsCardinalitiesV >= thresholdN ;
        focusedSetsCell = setsCell(indexFocusedSetsLogical) ;
        focusedLabelsVector = setsLabelsV(indexFocusedSetsLogical) ;
        focusedN = sum(indexFocusedSetsLogical) ;
		%{
        intVector = zeros(focusedNo, 1) ;
        for iSet = 1 : focusedN 
            intVector(iSet) = sum(ismember(focusedSetsCell{iSet}, remainingElementsVector)) ; 
        end
		%}
		%* Vectorized code
		intVector = cellfun(@(a,b) sum(ismember(a,b)), focusedSetsCell , num2cell(repmat(remainingElementsVector,focusedN,1),2)) ;
        
		%* selectedSetPosN is the position of the BIGGEST set among the ones with maximal number of uncovered elements. 
        [dummy, selectedSetPosN] = max(intVector) ;
        % Add to solution
        solution_setsCell = [solution_setsCell; focusedSetsCell{selectedSetPosN}] ;
        solution_setsLabelsV = [solution_setsLabelsV; focusedLabelsVector(selectedSetPosN)] ;
        % Update remainingElementsVector  
        remainingElementsVector = remainingElementsVector(~ismember(remainingElementsVector, focusedSetsCell{selectedSetPosN})) ;
        %remainingElementsVector = setdiff(remainingElementsVector, focusedSetsCell{selectedSetPosN}) ;
        % Remove selected set from list from which we pick the next ones
        relativePosV = find(indexFocusedSetsLogical) ;
        setsCell(relativePosV(selectedSetPosN)) = [ ] ;
        setsLabelsV(relativePosV(selectedSetPosN)) = [ ] ;
        setsCardinalitiesV(relativePosV(selectedSetPosN)) = [ ] ;
		%
		countN = countN + 1 ;
        %display(['Step ' num2str(countN) ', Elements Covered: ' num2str(sum(ismember(OriginalElementsV, unique([solution_setsCell{:}]))))])
    end
	%disp(' ')
	totalSetsN = numel(solution_setsLabelsV) ;
	%display(['Total sets: ' num2str(totalSetsN)])
	%disp('Algorithm terminated; check if a set is not necessary')
	setN = totalSetsN ;
	%disp(' ')
	for iSet  = totalSetsN : -1 : 1
		indexSetsL = true(setN,1) ;
		indexSetsL(iSet) = false ;
		elementsCoveredWithoutCurrentV = unique([solution_setsCell{indexSetsL}]) ;
		if numel(elementsCoveredWithoutCurrentV) == elementsN
			%disp(['Set ' num2str(iSet) ' is unnecessary'])
			solution_setsCell(iSet) = [ ] ;
			solution_setsLabelsV(iSet) = [ ] ;
			setN = setN - 1 ;
		end
	end
	if totalSetsN > setN
		%disp([num2str(totalSetsN - setN) ' sets were superfluous'])
		%Update number of sets
		totalSetsN = setN ;
		%display(['Total sets: ' num2str(totalSetsN)])
	else
		%disp('No superfluous set')
	end
end

%*=*=*=*=* 
%* Sort solution, according to the labels
[solution_setsLabelsV, sortedSolutionIndexVector] = sort(solution_setsLabelsV) ;
solution_setsCell = solution_setsCell(sortedSolutionIndexVector) ; 
%end-cellSCP-function

function [solution_A, solution_setsLabelsV] = matrixSCP(A, setsLabelsV, setsCardinalitiesV)
%MATRIX Version

%*=*=*=*=*
%*Identify elements covered by just one set and select them

%* UniquelyCoveredL(j) = true iff element j is covered by just one set of A
UniquelyCoveredL = sum(A, 2) == 1 ;
if any(UniquelyCoveredL)
	%* toBeSelectL(i) = true iff set i (i.e. A(:,i) is the only one that covers a certain
	toBeSelectL = logical(sum(A(UniquelyCoveredL,:), 1)) ;
    %display([num2str(sum(UniquelyCoveredL)) ' elements are uniquely covered'])
    %display([num2str(sum(toBeSelectL)) ' sets would be selected'])
	% add to solution
	solution_A = A(:, toBeSelectL) ;
	solution_setsLabelsV = setsLabelsV(toBeSelectL) ;  
	% update
	A(:, toBeSelectL) = [ ] ;
	setsLabelsV(toBeSelectL) = [ ] ;
	setsCardinalitiesV(toBeSelectL) = [ ] ;
	remainingElementsL = sum(solution_A, 2) == 0 ;
	%display(['Elements Covered: ' num2str(sum(~remainingElementsL))])
else
	% initialize solution
	solution_A = [ ] ;
	solution_setsLabelsV = [ ] ;  
	%
	remainingElementsL = true(size(UniquelyCoveredL)) ;
end

%*=*=*=*=* 
%* Run greedy algorithm

if ~any(remainingElementsL)
	%disp('No algorithm needed')
else
	%disp(' ')
	%remainingElementsL = remainingElementsL(:) ;
	%Sort Sets in descending order according to their cardinality
	[setsCardinalitiesV, sortedIndexVector] = sort(setsCardinalitiesV, 'descend') ;
	setsLabelsV = setsLabelsV(sortedIndexVector) ;
	A = A(:, sortedIndexVector) ;
    %disp('Greedy algorithm: start')
	countN = 0 ;
    while any(remainingElementsL)
		%* Trick to speed up: if a set (the first, for instance) contains n unconvered elements, 
		% all the sets having cardinality < n will contain less uncovered elements than this set, so that I do not need to compute intVector for them because they can't be better that the set selected at the beginning. 
        thresholdN = sum(A(:,1)&remainingElementsL) ;

        indexFocusedSetsLogical = setsCardinalitiesV >= thresholdN ;
        focusedSetsA = A(:,indexFocusedSetsLogical) ;
        focusedLabelsV = setsLabelsV(indexFocusedSetsLogical) ;

		% Count remaining elements covered by each set  of focusedSetA
		intV = (+remainingElementsL)' * (+focusedSetsA) ;

		%* selectedSetPosN is the position of the BIGGEST set among the ones with maximal number of uncovered elements. 
        [dummy, selectedSetPosN] = max(intV) ;
        % Add to solution
        solution_A = [solution_A focusedSetsA(:, selectedSetPosN)] ;
        solution_setsLabelsV = [solution_setsLabelsV ; focusedLabelsV(selectedSetPosN)] ;
        % Update remainingElementsL  
        remainingElementsL = remainingElementsL & ~focusedSetsA(:, selectedSetPosN) ;
        % Remove selected set from list from which we pick the next ones
        relativePosV = find(indexFocusedSetsLogical) ;
        A(:, relativePosV(selectedSetPosN)) = [ ] ;
        setsLabelsV(relativePosV(selectedSetPosN)) = [ ] ;
        setsCardinalitiesV(relativePosV(selectedSetPosN)) = [ ] ;
		%
		countN = countN + 1 ;
        %display(['Step ' num2str(countN) ', Elements Covered: ' num2str(sum(~remainingElementsL))])
	end
	%disp(' ')
	totalSetsN = numel(solution_setsLabelsV) ;
	%display(['Total sets: ' num2str(totalSetsN)])
	%disp('Algorithm terminated; check if a set is not necessary')
	setN = totalSetsN ;
	%disp(' ')
	for iSet  = totalSetsN : -1 : 1
		tempA = solution_A ;
		tempA(:, iSet) = [ ] ;
		if all(sum(tempA, 2)) % removing set iSet, all the elements are still covered
			%disp(['Set ' num2str(iSet) ' is unnecessary'])
			solution_A = tempA ;
			solution_setsLabelsV(iSet) = [ ] ;
			setN = setN - 1 ;
		end
	end
	if totalSetsN > setN
		%disp([num2str(totalSetsN - setN) ' sets were superfluous'])
		%Update number of sets
		totalSetsN = setN ;
		%display(['Total sets: ' num2str(totalSetsN)])
	else
		%disp('No superfluous set')
	end
end

%*=*=*=*=* 
%* Sort solution, according to the labels
[solution_setsLabelsV, sortedSolutionIndexVector] = sort(solution_setsLabelsV) ;
solution_A = solution_A(:, sortedSolutionIndexVector) ; 
%end-matrixSCP-function


