% Author: Kenneth Young 2017
% Setup Assembly Line Balancing and Scheduling Problem
% Constraint Programming model
% Attempt: 2
  
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

% disjunctive (loose)
constraint
	disjunctive(
		s,
		dur,
	);

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
ann: start_s = int_search(s, first_fail, indomain_min, complete);

% Required Initialisations
ann: startpair_s;
ann: start_Then_startpair;
ann: startpair_Then_start;
ann: start_io;
ann: start_sl; 
ann: start_ff;
ann: io_Then_startpair;
ann: s_Then_startpair;
ann: sl_Then_startpair;
ann: ff_Then_startpair;
ann: priority_input_order;
ann: priority_smallest;
ann: priority_smallest_largest;
ann: priority_first_fail;

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
