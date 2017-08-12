# Masters Research Project
# Kenneth Young
# SSBF MIP Model of the SUALBSP-2

# This file contains:
# 	-A MIP model of the SUALBSP-2
# 	-This model was adapted from Esmaeilbeigi et al. (2016),
# 	specifically their SSBF-2 model.

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
from solChecker import *

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

CHECK_SOLUTION = args.check_solution # create this file

TIMELIMIT = args.time_limit
EXPERIMENT_TOKEN = args.experiment_token

if args.very_quiet:
	args.quiet = True

class SolverSSBF:
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
		# self.model.setParam('Threads',1)

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

		# initialize g variables
		self.gs = self.model.addVars([ (i,j,k)
									   for k in self.inst.stations
									   for i in self.inst.feasibleTasks[k]
									   for j in self.inst.feasibleTasks[k].intersection(self.inst.followForw[i]) ],
									 vtype=GRB.BINARY, name='g')

		# initialize h variables
		self.hs = self.model.addVars([ (i,j,k)
									   for k in self.inst.stations
									   for i in self.inst.feasibleTasks[k]
									   for j in self.inst.feasibleTasks[k].intersection(self.inst.followBack[i]) ],
									 vtype=GRB.BINARY, name='h')

		# initialise z variables
		self.zs = self.model.addVars(self.inst.tasks, vtype=GRB.CONTINUOUS, name='z')

		# initialize r variables
		self.rs = self.model.addVars(self.inst.tasks, vtype=GRB.CONTINUOUS, name='r')

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

		# pdb.set_trace()

		# (39)
		self.model.addConstrs((   sum([ self.gs[i,j,k] 
										for j in self.inst.feasibleTasks[k].intersection(self.inst.followForw[i]) ])
								+ sum([ self.hs[i,j,k] 
										for j in self.inst.feasibleTasks[k].intersection(self.inst.followBack[i]) ])
								== self.xs[i,k]
								for i in self.inst.tasks
								for k in self.inst.feasibleStations[i] ),
								'(39)')

		# (40)
		self.model.addConstrs((   sum([ self.gs[i,j,k] 
										for i in self.inst.feasibleTasks[k].intersection(self.inst.precedeForw[j]) ])
								+ sum([ self.hs[i,j,k] 
										for i in self.inst.feasibleTasks[k].intersection(self.inst.precedeBack[j]) ])
								== self.xs[j,k]
								for j in self.inst.tasks
								for k in self.inst.feasibleStations[j] ),
								'(40)')

		# (41)
		self.model.addConstrs(( sum([ sum([ self.hs[i,j,k]
											for j in self.inst.feasibleTasks[k].intersection(self.inst.followBack[i]) ])
											for i in self.inst.feasibleTasks[k] ])
								== 1
								for k in self.inst.stations ),
								'(41)')

		# (43)
		self.model.addConstrs(( self.rs[i] + 1 + 
								(sum([ self.gs[i,j,k]-1
								 for k in self.inst.feasibleStations[i].intersection(self.inst.feasibleStations[j]) ]) )*( self.inst.numTasks
																										-len(self.inst.allSuccessors[i])
																										-len(self.inst.allPredecessors[j]) )
								<= self.rs[j]
								for i in self.inst.tasks
								for j in self.inst.followForw[i] ),
								'(43)')

		# (44)
		self.model.addConstrs(( self.rs[i] + 1 <= self.rs[j]
								for (i,j) in self.inst.precList ),
								'(44)')

		# (45)
		self.model.addConstrs(( self.zs[i] <= self.zs[j]
								for (i,j) in self.inst.precList ),
								'(45)')

		# (46)
		self.model.addConstrs((   sum([ self.inst.procList[i]*self.xs[i,k]
										for i in self.inst.feasibleTasks[k] ])
								+ sum([ sum([ self.inst.forwSU[i][j]*self.gs[i,j,k]
											for j in self.inst.feasibleTasks[k].intersection(self.inst.followForw[i]) ])
											for i in self.inst.feasibleTasks[k] ])

								+ sum([ sum([ self.inst.backSU[i][j]*self.hs[i,j,k]
											for j in self.inst.feasibleTasks[k].intersection(self.inst.followBack[i]) ])
											for i in self.inst.feasibleTasks[k] ])
								<= self.cycleTime 
								for k in self.inst.stations),
								'(46)')

		# (48)
		self.model.addConstrs(( sum([ self.xs[i,k]
									for i in self.inst.feasibleTasks[k].difference(set([j])) ])
								<= (self.inst.numTasks - self.inst.numStations + 1)*(1 - self.hs[j,j,k])
								for k in self.inst.stations
								for j in self.inst.feasibleTasks[k] ),
								'(48)')

		# (51) Bounds for r variable
		self.model.addConstrs(( self.rs[i] <= self.inst.numTasks - len(self.inst.allSuccessors[i])
								for i in self.inst.tasks ),
								'(51a)')
		self.model.addConstrs(( self.rs[i] >= 1 + len(self.inst.allPredecessors[i])
								for i in self.inst.tasks ),
								'(51b)')

		# Bounds for the cycle time
		self.model.addConstr(self.cycleTime <= self.inst.maxCycleTime, 'cycleTimeUB')
		self.model.addConstr(self.cycleTime >= self.inst.minCycleTime, 'cycleTimeLB')

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
		# store the ordering of the tasks
		self.taskSequence = [[] for k in self.inst.stations]
		for k in self.inst.stations:
			for i in self.taskAssignment[k]:
				self.taskSequence[k].append(round(self.rs[i].x))
		# pdb.set_trace()

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
					# print('!    Load = \t{}'.format(self.stationLoad[k]))
					print('!    Tasks = \t{}'.format(sorted(self.taskAssignment[k])))
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
	s = SolverSSBF(inst)
	if not args.very_quiet:
		print('Solving the MIP...')
	s.optimise()

	# output
	s.process_solution_statistics()
	if s.solFeasible == 1:
		s.print_solution()

	# check solution consistency
	if CHECK_SOLUTION:
		print('\nVerifying solution:')
		check_solution_SSBF(s.optimalCycleTime, s.taskAssignment, s.gs, s.hs, inst)

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