# Masters Research Project
# Kenneth Young
# FSBF MIP Model of the SUALBSP-2

# This file contains:
# 	-A MIP model of the SUALBSP-2
# 	-This model was adapted from Esmaeilbeigi et al. (2016),
# 	specifically their FSBF-2 model.

# Packages
import sys
import pdb 
import time
# import itertools
import csv
import argparse
import numpy as np
import networkx as nx
from gurobipy import *

# User-defined Packages
from ALB_instance_storage import AssemblyLineInstance

# initilise settings for argument parser
parser = argparse.ArgumentParser()
parser.add_argument('file', help='instance file')
parser.add_argument('-q', '--quiet', help='Some output', action='store_true')
parser.add_argument('-vq', '--very-quiet', help='Minimal output', action='store_true')
parser.add_argument('-s', '--statistics', help='Print statistics', action='store_true')
parser.add_argument('-c', '--check-solution', help='Check solution', action='store_true')
parser.add_argument('-H', '--human-readable', help='Human readable output', action='store_true')
parser.add_argument('-b', '--backwardSU-type', type=str, default='copy-forward-setups',
					help='Type of backward setup times to use. Options include:' 
						 '{just-forward-setups}')
parser.add_argument('-t', '--time-limit', type=float, default=1800,
					help='Optimisation time limit.')
parser.add_argument('-et', '--experiment-token', type=int, default=0,
					help='Indicator for which experiment is being run')
parser.add_argument('-v1', '--valid-ineq-1', help='Use this valid inequality', action='store_true')
parser.add_argument('-v2', '--valid-ineq-2', help='Use this valid inequality', action='store_true')
args = parser.parse_args()

# Define globals constants
if sys.platform == "win32":
	INST_DIR = 'instances\\'
elif sys.platform =="cygwin":
	INST_DIR = 'instances\\'
elif sys.platform == "darwin":
	INST_DIR = 'instances/'
elif sys.platform == "linux" or sys.platform == "linux2":
	INST_DIR = 'instances/'

# argument definitions
PRINT_STATISTICS = args.statistics
BACKWARD_SETUP_TYPE = args.backwardSU_type
VALID_INEQ_1 = args.valid_ineq_1
VALID_INEQ_2 = args.valid_ineq_2

CHECK_SOLUTION = 'tools/bin/check-solution' # create this file

TIMELIMIT = args.time_limit
EXPERIMENT_TOKEN = args.experiment_token

if args.very_quiet:
	args.quiet = True

class SolverFSBF:
	def __init__(self, inst):
		self.inst = inst
		# initialise the full MIP model
		self.model = Model('assemblyline')
		self.init_model_parameters()
		if args.quiet:
			self.model.setParam('OutputFlag', 0)
		if not args.very_quiet:
			print('Initialising the MIP...', end='')
		self.optimisation_times = []
		self.sequencing_times = []
		start = time.time()
		self.init_vars()
		self.create_objective()
		self.create_constraints()
		self.init_time = time.time() - start
		self.optimisation_times.append(self.init_time)
		if not args.very_quiet:
			print(' initialisation complete ({:.3f}s)'.format(self.init_time))

	def init_model_parameters(self):
		self.model.setParam('TimeLimit', TIMELIMIT)

	def init_vars(self):
		self.cycleTime = self.model.addVar(lb=self.inst.minCycleTime, 
										   ub=self.inst.maxCycleTime,
										   obj=1.0,
										   vtype=GRB.CONTINUOUS,
										   name='c')
		# initialise x variables
		self.xs = self.model.addVars([ (i,k) 
									   for i in self.inst.tasks
									   for k in self.inst.feasibleStations[i] ],
									 vtype=GRB.BINARY, name='x')
		
		# initialise y variables
		self.ys = self.model.addVars([ (i,j) 
									   for i in self.inst.tasks
									   for j in self.inst.followForw[i] ],
									 vtype=GRB.BINARY, name='y')

		# initialise w variables
		self.ws = self.model.addVars([ (i,j) 
									   for i in self.inst.tasks
									   for j in self.inst.followBack[i] ],
									 vtype=GRB.BINARY, name='w')

		# initialise z variables
		self.zs = self.model.addVars(self.inst.tasks, vtype=GRB.CONTINUOUS, name='z')

		# initialise start time variables, s
		self.ss = self.model.addVars(self.inst.tasks, vtype=GRB.CONTINUOUS, name='s')

		# initialise o variables
		self.os = self.model.addVars([ (i,k) 
									   for i in self.inst.tasks
									   for k in self.inst.feasibleStations[i] ],
									 vtype=GRB.BINARY, name='o')

		self.model.update()

	def create_objective(self):
		self.objective = self.cycleTime
		self.model.setObjective(self.objective, GRB.MINIMIZE)

	def create_constraints(self):
		# Each task i is assigned exactly one station
		self.model.addConstrs((self.xs.sum(i,'*') == 1
							   for i in self.inst.tasks), 'oneStationPerTask')

		# Encode the index of the stations which task i is assigned
		self.model.addConstrs((sum([k*self.xs[i,k] for k in self.inst.stations]) == self.zs[i]
							   for i in self.inst.tasks), 'encodeStationNums')

		# Each task i has exactly one successor (in forward and backward station loads)
		self.model.addConstrs((   self.ys.sum(i,'*') 
								+ self.ws.sum(i,'*') == 1
							   for i in self.inst.tasks), 'oneSuccessor')

		# Each task j has exactly one predecessor (in forward and backward station loads)
		self.model.addConstrs((   sum([self.ys.sum(i,j) for i in self.inst.followBack[j]]) 
								+ sum([self.ws.sum(i,j) for i in self.inst.precedeBack[j]]) == 1
							   for j in self.inst.tasks), 'onePredecessor')

		# Forward load: tasks contained in the same cycle are assigned the same station
		self.model.addConstrs((self.zs[j] - self.zs[i] <= self.inst.bigM * (1 - self.ys[i,j])
							   for i in self.inst.tasks 
							   for j in self.inst.followForw[i]), 'sameCycleForwA')
		self.model.addConstrs((self.zs[i] - self.zs[j] <= self.inst.bigM * (1 - self.ys[i,j])
							   for i in self.inst.tasks 
							   for j in self.inst.followForw[i]), 'sameCycleForwB')

		# Backward load: tasks contained in the same cycle are assigned the same station
		self.model.addConstrs((self.zs[j] - self.zs[i] <= self.inst.bigM * (1 - self.ws[i,j])
							   for i in self.inst.tasks 
							   for j in self.inst.followBack[i]), 'sameCycleBackA')
		self.model.addConstrs((self.zs[i] - self.zs[j] <= self.inst.bigM * (1 - self.ws[i,j])
							   for i in self.inst.tasks 
							   for j in self.inst.followBack[i]), 'sameCycleBackB')

		# Task i can only be the last task of station k if it is assigned to k
		self.model.addConstrs((self.os[i,k] <= self.xs[i,k]
							   for i in self.inst.tasks
							   for k in self.inst.feasibleStations[i]), 'onlyLastIfAssigned')

		# o[i,k] gets the value 1 only if task i is the last task of k
		self.model.addConstrs((self.ws.sum(i,'*') <= self.os.sum(i,'*')
							   for i in self.inst.tasks), 'ifBackSUthenLast')

		# In combination with 'oneSucc' and 'onePred' constraints, each station has only one o[i,k]==1
		self.model.addConstrs((sum([ self.os[i,k] for i in self.inst.feasibleTasks[k] ]) <= 1
							   for k in self.inst.stations), 'onlyOneLast')

		# Strengthening Knapsack constraint: Total load of each station is less than cycle time
		self.model.addConstrs((sum([ self.inst.procList[i]*self.xs[i,k] for i in self.inst.feasibleTasks[k] ]) <= self.cycleTime
							   for k in self.inst.stations), 'loadLessThanCycleTime')

		# The last task of each station finishes by the cycle time
		self.model.addConstrs((  sum([ self.inst.backSU[i][j]*self.ws[i,j] for j in self.inst.followBack[i] ])
							   + self.ss[i] + self.inst.procList[i] <= self.cycleTime
							   for i in self.inst.tasks), 'lastTaskFinishByCycleTime')

		# The number of backward setups is at least the number of stations
		self.model.addConstr(sum([ self.ws.sum(i,'*') for i in self.inst.tasks ]) <= self.inst.numStations,
							  'numBackSUsAtLeastNumStations')
		
		# Precedence Relations are respected in the forward direction
		self.model.addConstrs((  self.ss[i] + self.inst.maxCycleTime*(self.zs[i] - self.zs[j]) + self.inst.procList[i]
							   + self.inst.forwSU[i][j]*self.ys[i,j] <= self.ss[j]
							   for (i,j) in self.inst.precList), 'precedenceRelations')

		# Constrain start times of tasks following task i in the forward direction
		self.model.addConstrs((  self.ss[i] + (self.inst.procList[i] + self.inst.forwSU[i][j])
							   + (self.inst.maxCycleTime + self.inst.forwSU[i][j])*(self.ys[i,j] - 1) <= self.ss[j]
							   for i in self.inst.tasks
							   for j in set(self.inst.followForw[i]) - set(self.inst.precGraph.successors(i))),
							   '')

		# Bounds for the cycle time
		self.model.addConstr(self.cycleTime <= self.inst.maxCycleTime, 'cycleTimeUB')
		self.model.addConstr(self.cycleTime >= self.inst.minCycleTime, 'cycleTimeLB')
		
		# Valid Inequality:
		if VALID_INEQ_1:
			self.model.addConstrs(( self.ys[i,j] + self.ys[j,i] <= 1
									for i in self.inst.tasks
									for j in self.inst.followForw[i].intersection(self.inst.precedeForw[i])),
									'validiIneq1')

		# Valid Inequality: lower bound on the total line capacity
		if VALID_INEQ_2:
			self.model.addConstr(  sum([  sum([ self.inst.forwSU[i][j]*self.ys[i,j]
										for j in self.inst.followForw[i] ])
										for i in self.inst.tasks ])
								 + sum([ sum([ self.inst.backSU[i][j]*self.ws[i,j]
								 		for j in self.inst.followBack[i] ]) 
										for i in self.inst.tasks ])
								 + sum( self.inst.procList ) <= self.inst.numStations*self.cycleTime,
								 'lineCapacityLowerBound')

	def optimise(self):
		start = time.time()
		self.model.optimize()
		self.optimisation_times.append(time.time() - start)
		# # record the results for outputting
		

	def store_results_summary(self):
		# store assignment of tasks to stations
		self.taskAssignment = [None for k in self.inst.stations]
		for k in self.inst.stations:
			self.taskAssignment[k] = { i for i in self.inst.tasks if self.xs[i,k].x > 0.5 }
		# store start times of all tasks
		self.startTimes = [None for k in self.inst.stations]
		for k in self.inst.stations:
			self.startTimes[k] = [ round(self.ss[i].x) for i in sorted(self.taskAssignment[k])]
		# store load of each station
		self.stationLoad = [None for k in self.inst.stations]
		for k in self.inst.stations:
			try:
				self.stationLoad[k] = round(max([  self.ss[i].x + self.inst.procList[i] 
												  + sum([ self.inst.backSU[i][j]*self.ws[i,j].x 
														  for j in self.inst.followBack[i] ])
												  for i in self.taskAssignment[k] ]))
			except ValueError:
				self.stationLoad[k] = 0

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
	# OUTPUT METHODS
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
	def process_solution_statistics(self):
		# pdb.set_trace()
		if self.model.solcount == 0:
			self.solFeasible = 0
			self.solOptimal = 0
		else:
			self.solFeasible = 1

		# find out how model terminated. with feasible sol? with optimal sol?
		if self.solFeasible == 1:
			if self.model.status == 2:
				self.solOptimal = 1
				self.store_results_summary()		
			else:
				self.solOptimal = 0
				self.store_results_summary()
			self.optimalCycleTime = round(self.model.objval)
		else:
			self.optimalCycleTime = 0

		# get gap value at time of termination (either succesfful or timeout)
		self.gap = round(100*self.model.mipgap,2)
		self.statsTotalRuntime = round(sum(self.optimisation_times),4)
		# get number of nodes
		self.statsTotalNodes = round(self.model.nodecount)


	def print_solution(self):
		if args.human_readable:
			if not args.very_quiet:
				print('\n! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ')
				print('! \tSOLUTION ')
			else:
				print('\n',end='')
			print('! Cycle Time:\t{}'.format(round(self.model.objval)))
			if not args.very_quiet:
				for k in self.inst.stations:
					print('! Station {}'.format(k))
					print('!    Load = \t{}'.format(self.stationLoad[k]))
					print('!    Tasks = \t{}'.format(sorted(self.taskAssignment[k])))
					print('!    Starts = \t{}'.format(self.startTimes[k]))
		else:
			print(self.model.objval)

	def save_solution(self, results_file):
		# add stuff to a text file
		# if args.check_solution:
		# 	print('checking solution:')
		# 	os.system('%s %s %s' % (CHECK_SOLUTION, self.inst.filename, results_file))
		# pdb.set_trace()	
		with open(results_file, 'w', newline='') as csvfile:
			results = csv.writer(csvfile)
			results.writerow([self.solFeasible, self.solOptimal,
								self.optimalCycleTime, self.gap,
								self.statsTotalRuntime, self.statsTotalNodes])

	def print_statistics(self):
		if args.human_readable:
			print('\n! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ')
			print('! \tSOLUTION STATISTICS ')
			print('! Feasible Solution:\t{}'.format(self.solFeasible))
			print('! Optimal Solution:\t{}'.format(self.solOptimal))
			print('! Nodes Explored:\t{}'.format(int(self.statsTotalNodes)))
			if self.solFeasible:
				print('! Gap:\t\t\t{:.2f}'.format(self.gap))
			else:
				print('! Gap:\t\t\tNA')
			print('\n! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ')
			print('! \tRUNTIME STATISTICS ')
			print('! Init time:\t{:.4f}'.format(self.init_time))
			print('! Total:\t{:.4f}'.format(self.statsTotalRuntime))
			print('! Maximum:\t{:.4f}'.format(max(self.optimisation_times)))
			print('! Average:\t{:.4f}'.format(sum(self.optimisation_times)/len(self.optimisation_times)))
		else:
			print(self.solOptimal)
			print(self.nodesExplored)
			print(self.gap)
			print(self.init_time)
			print(sum(self.optimisation_times))
			print(max(self.optimisation_times))
			print(sum(self.optimisation_times)/len(self.optimisation_times))

# Script to create instance class, run the solver and output the solution
if __name__ == '__main__':
	# start total runtime timer
	start = time.time()

	filename = args.file # retrieve filename of instance to solve
	if args.human_readable:
		print('instance:', filename)
	else:
		print(filename)

	# store assembly line instance data
	inst = AssemblyLineInstance(INST_DIR,filename)
	# pdb.set_trace()

	# create Solver for given instance and optimise it
	s = SolverFSBF(inst)
	if not args.very_quiet:
		print('Solving the MIP...')
	s.optimise()

	# output
	s.process_solution_statistics()
	if s.solFeasible == 1:
		s.print_solution()
		
	s.save_solution('summary_results_{}.txt'.format(EXPERIMENT_TOKEN))
	if PRINT_STATISTICS:
		s.print_statistics()

	# print total runtime
	end = time.time()

	if not args.quiet:
		if args.human_readable:
			print('\n! Language runtime:\t{:.4f}'.format(end-start-sum(s.optimisation_times)))
		else:
			print(end-start-sum(s.optimisation_times))
	if args.human_readable:
		if args.very_quiet:
			print('')
		print('! Total runtime:\t{:.4f}'.format(end-start))
	else:
		print(end-start)

# EOF #