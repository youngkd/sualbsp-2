% Author: Kenneth Young 2017
% Setup Assembly Line Balancing and Scheduling Problem
% Constraint Programming model
% Attempt: 4
%	-adding circuit constraints
% NOTE: UNFINISHED
 
include "circuit.mzn";
include "cumulative.mzn";
include "disjunctive.mzn";
include "redefinitions.mzn";

0..1: full_output;
%full_output = 1;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% INSTANCE INITIALISATION
int: nTasks;
int: nPrecs;
int: maxLoad;  % maximum makespan

set of int: TASK;
set of int: PREC = 1..nPrecs;
set of int: TIME = 0..maxLoad;

array[TASK] of int: dur; % duration
array[TASK] of set of TASK: suc; % set of successors
array[TASK,TASK] of int: forwSU; % forward setup times
array[TASK,TASK] of int: backSU; % backward setup times
array[TASK] of set of TASK: followForw; % allowed followers in forward load
array[TASK] of set of TASK: followBack; % allowed followers in backward load
array[TASK] of set of TASK: precedeForw; % allowed preceders in forward load
array[TASK] of set of TASK: precedeBack; % allowed preceders in backward load

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% DECISION VARIABLES
array[TASK] of var TIME: s; % start time
array[TASK,TASK] of var TIME: spair; % start time pairings

array[TASK,TASK] of var bool: y; % forward direction following
array[TASK,TASK] of var bool: z; % backward direction following

var TIME: load;  % load

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% CONSTRAINTS

% Only one follower in either station load direction
constraint
	forall (
		i in TASK
	)(
		  sum( j in followForw[i] )( y[i,j] )
		+ sum( j in followBack[i] )( z[i,j] )
		== 1 
	);

% Only one preceder in either station load direction
constraint
	forall (
		j in TASK
	)(
		  sum( i in precedeForw[j] )( y[i,j] )
		+ sum( i in precedeBack[j] )( z[i,j] )
		== 1 
	);

% Exactly one backward setup
constraint
	sum( 
		i in TASK, j in followBack[i]
	)(
		z[i,j]
	) == 1
	;

% Precedence constraints
constraint
	forall ( 
		i in TASK, j in suc[i] 
	)(
		s[i] + dur[i] + forwSU[i,j]*y[i,j] <= s[j]
	);

% Forward station load respects setup times
constraint
	forall (
		i in TASK, j in followForw[i] 
	)(
		y[i,j] <-> ( s[i] + dur[i] + forwSU[i,j] == s[j] )
	);

% Backward station load respects station load
constraint
	forall (
		i in TASK
	)(
		  s[i] + dur[i]
		+ sum( 
			j in followBack[i]
		  )(
		  	backSU[i,j]*z[i,j]
		  )
		<= load
	);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% REDUNDANT CONSTRAINTS

% Cumulative Global
constraint
	cumulative(
		[ spair[i,j] 			| i in TASK, j in TASK ],
		[ dur[i] + forwSU[i,j] 	| i in TASK, j in TASK ],
		[ y[i,j]				| i in TASK, j in TASK ],
		1
	);

constraint
	circuit(
		s
	);

% linking constraint
%constraint
%	forall(i in TASK, j in TASK)
%	(
%		y[i,j] -> (s[i] == spair[i,j])
%	);

%constraint
%	forall(i in TASK, j in TASK where i != j)
%	(
%		spair[i,j] == spair[j,i]
%	);

constraint
	forall( 
		i in TASK, j in TASK
	)(
		s[i] == spair[i,j]
	);

%constraint
%	forall(i in TASK, j in TASK)
%	(
%		z[i,j] -> ()
%	);
	

% add sysmmetery breaking

% Fix some ordering variables to zero
constraint
	forall (
		i in TASK, j in TASK
	where
		not( j in followForw[i] )
	)(
		y[i,j] == 0
	);
constraint
	forall (
		i in TASK, j in TASK
	where
		not( j in followBack[i] )
	)(
		z[i,j] == 0
	);
	

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% OBJECTIVE
% solve minimize makespan;

ann: default_s;

% Basic Searches
ann: start_s = int_search(s, smallest, indomain_min, complete);   
ann: startpair_s = int_search( [ spair[i,j] | i in TASK, j in TASK ],
							smallest, indomain_min, complete);

% Sequential Searches
ann: start_Then_startpair = seq_search([ start_s, startpair_s ]);
ann: startpair_Then_start = seq_search([ startpair_s, start_s ]);

ann: start_io = int_search(s, input_order, indomain_min, complete);  
ann: start_sl = int_search(s, smallest_largest, indomain_min, complete);  
ann: start_ff = int_search(s, first_fail, indomain_min, complete);  

ann: io_Then_startpair = seq_search([ start_io, startpair_s ]);
ann: s_Then_startpair = seq_search([ start_s, startpair_s ]);
ann: sl_Then_startpair = seq_search([ start_sl, startpair_s ]);
ann: ff_Then_startpair = seq_search([ start_ff, startpair_s ]);

% Priority Searches
ann: priority_input_order = 
	priority_search(s,
					[seq_search([
								int_search([s[i]], input_order, indomain_min, complete),
								int_search([ spair[i,j] | j in TASK ],
										input_order, indomain_min, complete)])
					| i in TASK ],
					input_order, complete);

ann: priority_smallest = 
	priority_search(s,
					[seq_search([
								int_search([s[i]], input_order, indomain_min, complete),
								int_search([ spair[i,j] | j in TASK ],
										smallest, indomain_min, complete)])
					| i in TASK ],
					smallest, complete);

ann: priority_smallest_largest = 
	priority_search(s,
					[seq_search([
								int_search([s[i]], input_order, indomain_min, complete),
								int_search([ spair[i,j] | j in TASK ],
										smallest_largest, indomain_min, complete)])
					| i in TASK ],
					smallest_largest, complete);

ann: priority_first_fail = 
	priority_search(s,
					[seq_search([
								int_search([s[i]], input_order, indomain_min, complete),
								int_search([ spair[i,j] | j in TASK ],
										first_fail, indomain_min, complete)])
					| i in TASK ],
					first_fail, complete);

ann: my_search;

% Solve
solve :: my_search
minimize load;
% satisfy;


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% OUTPUT
output
if full_output == 0 then    
  ["load = " ++ show(load) ++ "\n"]
elseif full_output == 1 then
  ["load = " ++ show(load) ++ "\n"] ++
  ["start = " ++ show(s) ++ "\n"]
else
  [""]
endif
;
