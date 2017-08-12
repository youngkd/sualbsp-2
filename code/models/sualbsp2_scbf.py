# Masters Research Project
# Kenneth Young
# SCBF MIP Model of the SUALBSP-2

# This file contains:
# 	-A MIP model of the SUALBSP-2
# 	-This model was adapted from Esmaeilbeigi et al. (2016),
# 	specifically their SCBF-2 model.

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

class SolverSCBF:
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

		# initialise start time variables, s
		self.ss = self.model.addVars(self.inst.tasks, vtype=GRB.CONTINUOUS, name='s')

		# initialize q variables
		self.qs = self.model.addVars([ (i,j)
										for i in self.inst.tasks
										for j in set(self.inst.tasks).difference(self.inst.allPredecessors[i].union(self.inst.allSuccessors[i]))
										if i != j],
									vtype=GRB.BINARY, name='q')

		# initialize r variables
		self.rs = self.model.addVars(self.inst.tasks, vtype=GRB.CONTINUOUS, name='r')

		self.model.update()

	def create_objective(self):
		self.objective = self.cycleTime
		self.model.setObjective(self.objective, GRB.MINIMIZE)

	def create_constraints(self):
		# (4) Each task i has exactly one successor (in forward and backward station loads)
		self.model.addConstrs((   self.ys.sum(i,'*') 
								+ self.ws.sum(i,'*') == 1
							   for i in self.inst.tasks), 'oneSuccessor')

		# (5) Each task j has exactly one predecessor (in forward and backward station loads)
		self.model.addConstrs((   sum([self.ys.sum(i,j) for i in self.inst.followBack[j]]) 
								+ sum([self.ws.sum(i,j) for i in self.inst.precedeBack[j]]) == 1
							   for j in self.inst.tasks), 'onePredecessor')

		# (30)
		self.model.addConstrs((   self.ss[i] + self.inst.procList[i]
								+ sum([ self.inst.backSU[i][j]*self.ws[i,j]
										for j in self.inst.followBack[i] ])
								<= self.cycleTime 
								for i in self.inst.tasks),
								'(30)')

		# (44)
		self.model.addConstrs(( self.rs[i] + 1 <= self.rs[j]
								for (i,j) in self.inst.precList ),
								'(44)')

		# (54)
		self.model.addConstrs(( self.qs[i,j] + self.qs[j,i] == 1
								for i in self.inst.tasks
								for j in set(self.inst.tasks).difference(self.inst.allPredecessors[i].union(self.inst.allSuccessors[i]))
								if i<j ),
								'(54)')

		# (55)
		self.model.addConstrs(( self.rs[i] + 1 + ( self.inst.numTasks
													-len(self.inst.allSuccessors[i])
													-len(self.inst.allPredecessors[j]) )*(self.qs[i,j]-1)
								<= self.rs[j]
								for i in self.inst.tasks
								for j in set(self.inst.tasks).difference(self.inst.allPredecessors[i].union(self.inst.allSuccessors[i]))
								if i!=j ),
								'(55)')

		# pdb.set_trace()

		# (56)
		self.model.addConstrs(( self.rs[j] - 1 + ( self.inst.numTasks
													-len(self.inst.allSuccessors[j])
													-len(self.inst.allPredecessors[i]) 
													-2 )*(self.ys[i,j]-1)
								<= self.rs[i]
								for i in self.inst.tasks
								for j in self.inst.followForw[i] ),
								'(56)')

		# (57)
		self.model.addConstrs(( self.ys[i,j] <= self.qs[i,j]
								for i in self.inst.tasks
								for j in self.inst.followForw[i].difference(self.inst.allSuccessors[i]) ),'(57)')

		# (58)
		self.model.addConstrs(( self.ws[j,i] <= self.qs[i,j]
								for i in self.inst.tasks
								for j in self.inst.precedeBack[i].difference(self.inst.allSuccessors[i])
								if i!=j ),'(58)')

		# (60)
		self.model.addConstrs((   self.ss[i] + (self.inst.procList[i] + self.inst.forwSU[i][j])
								+ (self.cycleTime + self.inst.forwSU[i][j])*(self.ys[i,j] - 1)
								<= self.ss[j]
								for i in self.inst.tasks
								for j in self.inst.followForw[i] ),'(60)')

		# (61)
		self.model.addConstr( sum([ sum([ self.ws[i,j]
											for j in self.inst.followBack[i] ])
											for i in self.inst.tasks ])
								== self.inst.numStations ,'(61)')

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
		# find the last task of each station
		finalTasks = []
		self.stationLoad = []
		for i in self.inst.tasks:
			for j in self.inst.followBack[i]:
				if self.ws[i,j].x > 0.5:
					finalTasks.append(i)
					self.stationLoad.append(round(self.ss[i].x + self.inst.procList[i] + self.inst.backSU[i][j]*self.ws[i,j].x,0))


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
					# print('!    Tasks = \t{}'.format(sorted(self.taskAssignment[k])))
					# print('!    Starts = \t{}'.format(self.startTimes[k]))
		else:
			print(self.model.objval)

	def save_solution(self, results_file):
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
			print(self.status)
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
	s = SolverSCBF(inst)
	if not args.very_quiet:
		print('Solving the MIP...')
	s.optimise()

	# output
	s.process_solution_statistics()
	if s.solFeasible == 1:
		s.print_solution()

	# check solution consistency
	# if CHECK_SOLUTION:
	# 	print('\nVerifying solution:')
	# 	check_solution_SCBF(s.optimalCycleTime, s.taskAssignment, s.gs, s.hs, inst)

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