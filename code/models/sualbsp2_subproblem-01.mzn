% Author: Kenneth Young 2017
% Setup Assembly Line Balancing and Scheduling Problem
% Constraint Programming model
% Attempt: 1
  
include "cumulative.mzn";
include "disjunctive.mzn";
include "redefinitions.mzn";

0..1: full_output;
% full_output = 1;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% INSTANCE INITIALISATION
int: nTasks;
int: nAssTasks;
int: nPrecs;
int: maxLoad;  % maximum makespan

set of int: TASK = 1..nTasks;
set of int: ASSTASK; % assigned tasks
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

array[TASK,TASK] of var bool: y; % forward direction following
array[TASK,TASK] of var bool: z; % backward direction following

var TIME: load;  % load

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% CONSTRAINTS

% Only one follower in either station load direction
constraint
	forall (
		i in ASSTASK
	)(
		  sum( j in followForw[i] where j in ASSTASK )( y[i,j] )
		+ sum( j in followBack[i] where j in ASSTASK )( z[i,j] )
		== 1 
	);

% Only one preceder in either station load direction
constraint
	forall (
		j in ASSTASK
	)(
		  sum( i in precedeForw[j] where i in ASSTASK )( y[i,j] )
		+ sum( i in precedeBack[j] where i in ASSTASK )( z[i,j] )
		== 1 
	);

% Exactly one backward setup
constraint
	sum( 
		i in ASSTASK, j in followBack[i]
	where
		j in ASSTASK
	)(
		z[i,j]
	) == 1
	;

% Precedence constraints
constraint
	forall ( 
		i in ASSTASK, j in suc[i] 
	where
		j in ASSTASK
	)(
		s[i] + dur[i] + forwSU[i,j]*y[i,j] <= s[j]
	);

% Forward station load respects setup times
constraint
	forall (
		i in ASSTASK, j in followForw[i] 
	where
		j in ASSTASK
	)(
		y[i,j] -> ( s[i] + dur[i] + forwSU[i,j] <= s[j] )
	);

% Backward station load respects station load
constraint
	forall (
		i in ASSTASK
	)(
		  s[i] + dur[i]
		+ sum( 
			j in followBack[i]
		  where
			j in ASSTASK 
		  )(
		  	backSU[i,j]*z[i,j]
		  )
		<= load
	);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% REDUNDANT CONSTRAINTS

% Cumulative Global (very loose)
%constraint
%	cumulative(
%		[ s[i] | i in ASSTASK ],
%		[ dur[i] | i in ASSTASK ],
%		[ 1 | i in ASSTASK ],
%		1
%	);

% Disjunctive Global (also loose)
%constraint 
%	disjunctive(
%		[ s[i] | i in ASSTASK ],
%		[ dur[i] | i in ASSTASK ]
%	);




% Fix the start times of unassigned tasks to zero
constraint
	forall (
		i in TASK
	where
		not( i in ASSTASK )
	)(
		s[i] == 0
	);

% Fix the ordering variables of unassigned tasks to zero
constraint
	forall (
		i in TASK, j in TASK
	where
		not( j in ASSTASK )
	)(
		y[i,j] == 0
	/\
		z[i,j] == 0
	/\
		y[j,i] == 0
	/\
		z[j,i] == 0
	);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% OBJECTIVE
% solve minimize makespan;

ann: default_s;

% Basic Searches
ann: start_s = int_search(s, first_fail, indomain_min, complete);   

ann: my_search;

% Solve
solve :: start_s
minimize load;
% satisfy;


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% OUTPUT
output
if full_output == 0 then    
  ["load = " ++ show(load) ++ "\n"]
elseif full_output == 1 then
  ["load = " ++ show(load) ++ "\n"] ++
  ["start = " ++ show(s) ++ "\n"] ++
  ["y = " ++ show(y) ++ "\n"] ++
  ["z = " ++ show(z) ++ "\n"]
else
  [""]
endif
;
