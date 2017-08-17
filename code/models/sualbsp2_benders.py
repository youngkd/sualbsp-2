# Project: 	Masters Research
# Author: 	Kenneth Young
# Title:	Benders decomposition of the SUALBSP-2

# This file contains:
# 	-Classes and methods for iterating through Benders
#	-MIP: Gurobi is the MIP solver used
#	-CP: chuffed compiled with MiniZinc is the CP solver used
#	-TSP: currently no dedicated TSP solver has been utilised

# Example call from the command line:
#	python sualbsp2_benders.py -H -s -q -gc -lc -ic "InstancePath/InstanceFile.alb"
# Explanation:
#	This will run the Benders algorithm with human readable output (-H),
#	output statistics (-s), quiet output (-q), global cuts (-gc), logic cuts (-lc)
#	and simple infer cuts (-ic) on the given instance file.

# Packages
import sys
import pdb 
import time
# import itertools
import csv
import re
import ast
import argparse
import numpy as np
import networkx as nx
from gurobipy import *

# User-defined Functionality
from ALB_instance_storage import AssemblyLineInstance
from callback_SubTourElim import *
from solChecker import *

# initilise settings for argument parser
parser = argparse.ArgumentParser()
parser.add_argument('file', help='Instance file')
parser.add_argument('-q', '--quiet', help='Some output', action='store_true')
parser.add_argument('-vq', '--very-quiet', help='Minimal output', action='store_true')
parser.add_argument('-s', '--statistics', help='Print statistics', action='store_true')
parser.add_argument('-c', '--check-solution', help='Check solution', action='store_true')
parser.add_argument('-H', '--human-readable', help='Human readable output', action='store_true')
parser.add_argument('-ste', '--sub-tour-elimination', type=int, default=0,
					help='Use sub-tour elimination lazy constraint generation when '
						 'solving the relaxed master problem. Specify the size of '
						 'sub-tours to eliminate.')
parser.add_argument('-nc', '--nogoods', help='Use nogood cuts', action='store_true')
parser.add_argument('-gb', '--global-bounds', help='Use global bounds', action='store_true')
parser.add_argument('-ic', '--infer-cuts', help='Use simple infer cuts', action='store_true')
parser.add_argument('-ic2', '--smart-infer-cuts', help='Use smart infer cuts', action='store_true')
parser.add_argument('-ic3', '--smartest-infer-cuts', help='Use smartest infer cuts', action='store_true')
parser.add_argument('-lc', '--logic-cuts', help='Use logic cuts for infeasible assignment', action='store_true')
parser.add_argument('-b', '--backwardSU-type', type=str, default='copy-forward-setups',
					help='Type of backward setup times. Options include: ' 
						 "'from-data'(default), 'just-forward-setups'")
parser.add_argument('-mi', '--max-iterations', type=int, default=sys.maxsize/2,
					help='Maximum number of Benders iterations')
parser.add_argument('-mpt', '--master-problem-type', type=str, default='ass',
					help='Type of master-problems to use. Options include:'
						  "'sched' and 'ass'(default)")
parser.add_argument('-spt', '--sub-problem-type', type=str, default='opt',
					help='Type of sub-problems to use. Options include:'
						  "'opt'(default) and 'feas'")
parser.add_argument('-sps', '--sub-problem-solver', type=str, default='mip',
					help='Type of sub-problem solver. Options include:'
						 "'mip'(default), 'cp2', 'cp3', 'tsp-solver'")
parser.add_argument('-t', '--time-limit', type=float, default=1800,
					help='Optimisation time limit')
parser.add_argument('-cps', '--cp-search', type=str, default='start_s',
					help='Search strategy to use when using CP to solve'
						 'the scheduing sub-problems. Options include:'
						  "'default', 'start_s'(default) and others")
parser.add_argument('-et', '--experiment-token', type=int, default=0,
					help='Indicator for which experiment is being run')
parser.add_argument('-ws', '--warm-start', help='Use warm-starts for RMP', action='store_true')
parser.add_argument('-tl', '--thread-limit', help='Restrict MIP-SP solver to 1 thread', action='store_true')
args = parser.parse_args()

# Define globals constants
if sys.platform == "win32":
	INST_DIR = 'instances\\'
	CHUFFED_DIR = 'chuffed\\windows\\'
elif sys.platform =="cygwin":
	INST_DIR = 'instances\\'
	CHUFFED_DIR = 'chuffed\\windows\\'
elif sys.platform == "darwin":
	INST_DIR = 'instances/'
	CHUFFED_DIR = 'chuffed/mac/'
elif sys.platform == "linux" or sys.platform == "linux2":
	INST_DIR = 'instances/'
	CHUFFED_DIR = 'chuffed/unix/'

# store the given arguments as globals
PRINT_STATISTICS = args.statistics
SIZE_OF_SUB_TOURS_TO_ELIMNATE = args.sub_tour_elimination
USE_NOGOODS = args.nogoods
USE_GLOBAL_BOUNDS = args.global_bounds
USE_INFER_CUTS = args.infer_cuts
USE_INFER_CUTS_SMART = args.smart_infer_cuts
USE_INFER_CUTS_SMARTEST = args.smartest_infer_cuts
USE_LOGIC_CUTS = args.logic_cuts
RMP_TYPE = args.master_problem_type
SUB_PROBLEM_TYPE = args.sub_problem_type
SUB_PROBLEM_SOLVER = args.sub_problem_solver
MAX_BENDERS_ITERATIONS = args.max_iterations
# BACKWARD_SETUP_TYPE = args.backwardSU_type
TIMELIMIT = args.time_limit
CHECK_SOLUTION = args.check_solution
SEARCH = args.cp_search
EXPERIMENT_TOKEN = args.experiment_token
WARM_START = args.warm_start
SP_THREAD_LIMIT = args.thread_limit

if args.very_quiet:
	args.quiet = True

# Class defining the instance of a particular station and its sub-problem
class Station:
	def __init__(self, inst, stationNum, tasks, curCycleTime, bestCycleTimeUB):
		self.inst = inst
		self.stationNum = stationNum
		self.tasks = tasks
		self.curCycleTime = curCycleTime
		self.bestCycleTimeUB = bestCycleTimeUB
		self.calculate_tour_maximum_naive()
		self.initialise()

	def __str__(self):
		return 'Station {} tasks: {}'.format(self.stationNum,self.tasks)

	def __repr__(self):
		return 'Station {} tasks: {}'.format(self.stationNum,self.tasks)

	def initialise(self):
		# create set of precedence relations
		self.precList = [ (i,j) for (i,j) in self.inst.precList
								if i in self.tasks
								and j in self.tasks ]

		self.construct_sets_of_allowed_followers_and_preceders_for_each_task()

	def construct_sets_of_allowed_followers_and_preceders_for_each_task(self):
		self.followForw = [ {j for j in self.inst.followForw[i]
							   if i in self.tasks
							   and j in self.tasks}
							for i in self.inst.tasks ]
		self.precedeForw = [ {j for j in self.inst.precedeForw[i]
								if i in self.tasks
								and j in self.tasks}
							 for i in self.inst.tasks ]
		self.followBack = [ {j for j in self.inst.followBack[i]
							   if i in self.tasks
							   and j in self.tasks}
							for i in self.inst.tasks ]
		self.precedeBack = [ {j for j in self.inst.precedeBack[i]
							   if i in self.tasks
							   and j in self.tasks}
							for i in self.inst.tasks ]

	def calculate_tour_minimum_naive(self):
		self.tourLB = 0

	def calculate_tour_maximum_naive(self):
		if len(self.tasks)==0:
			self.naiveLoadUB=0
		# elif len(self.tasks)==1:
		# 	self.naiveLoadUB = self.inst.procList[list(self.tasks)[0]]
		else:
			sortedTasks = sorted(list(self.tasks))
			pairedTasks = zip(sortedTasks[:-1], sortedTasks[1:])
			# pdb.set_trace()

			# assuming all tasks in down in serial
			sumOfProcList = sum([self.inst.procList[i] for i in self.tasks])
			# pdb.set_trace()
			# assuming tasks are completed in input order
			naiveFeasibleOrderingCost = sum( [ self.inst.forwSU[i][j] 
				   for (i,j) in pairedTasks ] ) + self.inst.backSU[sortedTasks[-1]][sortedTasks[0]]
										
			# together these give a naive feaible maximum
			self.naiveLoadUB = sumOfProcList + naiveFeasibleOrderingCost

	def reindex_precedence_relations(self):
		self.oldIndexedTasks = sorted(list(self.tasks))
		# make precedence graph
		self.originalPrecList = [ (i,j) for (i,j) in self.inst.precList
										if i in self.tasks
										if j in self.tasks ]

		# make list of reindexed precedence relations
		self.precList = [None for i in range(len(self.originalPrecList))]
		for i in range(len(self.originalPrecList)):
			prec = self.originalPrecList[i]
			self.precList[i] = (self.reindexedTasks[self.oldIndexedTasks.index(prec[0])],
								self.reindexedTasks[self.oldIndexedTasks.index(prec[1])] )

		# store precedence relations in an alternative manner
		self.altPrecList = [set([j for (idash,j) in self.precList if i==idash]) for i in self.reindexedTasks]

	def create_precedence_graph_for_station(self):
		# initialise the graph
		self.precGraph = nx.DiGraph()
		self.precGraph.add_nodes_from(self.reindexedTasks)
		# add non-trivial edges
		for prec in self.precList:
			self.precGraph.add_edge(prec[0],prec[1])

	def construct_sets_of_allowed_followers_and_preceders_for_each_task_assigned_this_station(self):
		self.stationFollowForw = [ (   set(self.reindexedTasks)
									 - (self.allSuccessors[i-1] - set(self.precGraph.successors(i)) )
									 - self.allPredecessors[i-1])
									 - set([i])
									for i in self.reindexedTasks ]
		self.stationPrecedeForw = [ set([ j for j in self.reindexedTasks 
									if i in self.stationFollowForw[j-1] ])
									for i in self.reindexedTasks ]
		self.stationFollowBack = [ set(self.reindexedTasks) - self.allSuccessors[i-1]
									for i in self.reindexedTasks ]
		self.stationPrecedeBack = [ set([ j for j in self.reindexedTasks 
									if i in self.stationFollowBack[j-1] ])
									for i in self.reindexedTasks ]

	def find_all_successors(self, node, allSuccessors):
		# (Recursive function)
		# finds all successors of a given node
		for succ in self.precGraph.successors(node):
			if succ not in allSuccessors:
				allSuccessors.add(succ)
				self.find_all_successors(succ, allSuccessors)
		return allSuccessors

	def find_all_predecessors(self, node, allPredecessors):
		# (Recursive function)
		# finds all predecessors of a given node
		for pred in self.precGraph.predecessors(node):
			if pred not in allPredecessors:
				allPredecessors.add(pred)
				self.find_all_predecessors(pred, allPredecessors)
		return allPredecessors

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
	# MIP SUB-PROBLEM
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

	def initialise_MIP(self, time_remaining):
		self.model = Model('station[%d]' %(self.stationNum))
		start = time.time()
		if args.quiet:
			self.model.setParam('LogToConsole', 0)
		self.model.setParam('TimeLimit', time_remaining)
		if SP_THREAD_LIMIT:
			self.model.setParam('Threads', 1)
		# create big-M value
		self.bigM = self.inst.maxCycleTime # self.curCycleTime
		self.init_MIP_vars()
		self.create_MIP_objective()
		self.create_MIP_constraints()
		self.init_time = time.time() - start

	def init_MIP_vars(self):
		self.load = self.model.addVar(ub=self.naiveLoadUB,
										obj=0.0,
										vtype=GRB.INTEGER,
										name='l')

		# Initialise y variables: Forward Sequencing
		self.ys = self.model.addVars([ (i,j)
									   for i in self.tasks
									   for j in self.followForw[i] ],
									 vtype=GRB.BINARY, name='y')

		# Initialise z variables: Backward Sequencing
		self.zs = self.model.addVars([ (i,j)
									   for i in self.tasks
									   for j in self.followBack[i] ],
									 vtype=GRB.BINARY, name='z')

		# Initialise s variables: Start times
		self.ss = self.model.addVars(self.tasks,
									 ub=self.naiveLoadUB, # ub=self.curCycleTime,
									 vtype=GRB.INTEGER,
									 name='s')

		self.model.update()

	def create_MIP_objective(self):
		# Satisfiability objective
		# self.model.setObjective()

		# Minimisation objective
		self.objective = self.load
		self.model.setObjective(self.objective, GRB.MINIMIZE)

	def create_MIP_constraints(self):
		# Each task has exactly one sucessor (in forward and backward direction)
		self.model.addConstrs(( self.ys.sum(i,'*') + self.zs.sum(i,'*') == 1
								for i in self.tasks), 'oneSuccessor')

		# Each task has exactly one predecessr ('')
		self.model.addConstrs((   sum([ self.ys[i,j] for i in self.precedeForw[j]
													 if (i,j) in self.ys])
								+ sum([ self.zs[i,j] for i in self.precedeBack[j]
													 if (i,j) in self.zs ]) == 1
								for j in self.tasks), 'onePredecessor')

		# There is exactly one backward setup
		self.model.addConstr( sum([ self.zs.sum(i,'*') for i in self.tasks ]) == 1,
							  'oneBackSU' )

		# Precedence Relations are respected with this station
		self.model.addConstrs((   self.ss[i] + self.inst.procList[i] 
								+ self.inst.forwSU[i][j]*self.ys[i,j] <= self.ss[j]
								for (i,j) in self.precList), 'precedenceRelations')

		# Forward Load: If y[i,j]==1 then task i finishes before task j starts
		self.model.addConstrs((   self.ss[i] + self.inst.procList[i]
								+ self.inst.forwSU[i][j] <= self.ss[j] + self.bigM*(1 - self.ys[i,j])
								for i in self.tasks
								for j in self.followForw[i]), 'forwardLoadStartTimes')

		# Backward Load: Last task of this station finishes by the station load
		self.model.addConstrs((   self.ss[i] + self.inst.procList[i]
								+ sum([ self.inst.backSU[i][j]*self.zs[i,j]
										for j in self.followBack[i] ])
								<= self.load
								for i in self.tasks), 'backwardLoadStartTime')

		if USE_LOGIC_CUTS:
			# Station Load: Best upper bound on the load from global bounds
			self.model.addConstr(self.load <= self.bestCycleTimeUB,
								'stationLoadUB')
		else:
			self.model.addConstr(self.load <= self.naiveLoadUB,
								'stationLoadUB')

	def solve_MIP(self):
		self.model.optimize()

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
	# CP SUB-PROBLEM
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

	def initialise_CP(self):
		start = time.time()
		if SUB_PROBLEM_SOLVER == 'cp2':
			tmp= 'sualbsp2_subproblem-02'
		elif SUB_PROBLEM_SOLVER == 'cp3':
			tmp= 'sualbsp2_subproblem-03'
		else:
			tmp= 'sualbsp2_subproblem-03'
		self.dznFile = 'subprob{}'.format(EXPERIMENT_TOKEN)
		self.statsFile = 'CPstats{}'.format(EXPERIMENT_TOKEN)
		self.solFile = 'CPsol{}'.format(EXPERIMENT_TOKEN)
		self.modelFile = tmp
		self.fullOutput = 1
		self.searchStrat = SEARCH
		# self.CPtimelimit = 600

		# create datazinc file
		self.store_station_data()
		self.write_dzn_sub_problem_file()

		# flatten datazinc file to a flatzinc file
		os.system("mzn2fzn -Gchuffed -D \"my_search = {}; full_output = {};\" "
				"{}.mzn {}.dzn".format(self.searchStrat,
									self.fullOutput,
									self.modelFile,
									self.dznFile))
		self.init_time = time.time() - start

	def store_station_data(self):
		# change indexing of assigned tasks
		self.reindexedTasks = range(1,len(self.tasks)+1)
		self.reindex_precedence_relations()
		self.create_precedence_graph_for_station()

		# find all predecessors and successors of each task
		self.allPredecessors = []
		for i in self.reindexedTasks:
			self.allPredecessors.append(self.find_all_predecessors(i,set()))
		self.allSuccessors = []
		for i in self.reindexedTasks:
			self.allSuccessors.append(self.find_all_successors(i,set()))

		self.construct_sets_of_allowed_followers_and_preceders_for_each_task_assigned_this_station()

	def OLD_OLD_write_dzn_sub_problem_file(self):
		# create datazinc file for cp solver
		with open(self.dznFile+'.dzn', 'w') as f:
			# store maximum load value
			f.write('% max load\nmaxLoad={};\n\n'.format(self.bestCycleTimeUB))

			# store total number of tasks and number of assigned tasks
			f.write('% total number of tasks\nnTasks={};\n\n'.format(self.inst.numTasks))
			f.write('% number of assigned tasks\nnAssTasks={};\n\n'.format(len(self.tasks)))

			# store list of assigned tasks
			f.write('% assigned tasks\nASSTASK={};\n\n'.format(set([ j+1 
																	for j in self.tasks ])))

			# store processing times
			f.write('% processing times\ndur={};\n\n'.format(self.inst.procList))

			# store all precedence relations
			f.write('% precedence relations\nnPrecs={};\nsuc = ['.format(len(self.inst.precList)))
			for i in self.inst.tasks:
				if self.inst.altPrecList[i]==set():
					if i == self.inst.tasks[-1]:
						f.write('{}')
					else:
						f.write('{}, ')
				else:
					f.write('{}, '.format(set([ j+1 for j in self.inst.altPrecList[i] ])))
			f.write('];\n\n')

			# store forward setup times
			f.write('% forward setup times\nforwSU=[')
			for i in self.inst.tasks:
				if i == self.inst.tasks[-1]:
					f.write('| {} |];\n\n'.format(str(list(self.inst.forwSU[i]))[1:-1]))
				else:
					f.write('| {}\n\t'.format(str(list(self.inst.forwSU[i]))[1:-1]))

			# store backward setup times
			f.write('% backward setup times\nbackSU=[')
			for i in self.inst.tasks:
				if i == self.inst.tasks[-1]:
					f.write('| {} |];\n\n'.format(str(list(self.inst.backSU[i]))[1:-1]))
				else:
					f.write('| {}\n\t'.format(str(list(self.inst.backSU[i]))[1:-1]))

			# store the set of forward followers
			f.write('% other sets to define\nfollowForw=[')
			for i in self.inst.tasks:
				if self.inst.followForw[i]==set():
					if i == self.inst.tasks[-1]:
						f.write('{}')
					else:
						f.write('{}, ')
				else:
					f.write('{}, '.format(set([ j+1 
												for j in self.inst.followForw[i] ])))
			f.write('];\n\n')

			f.write('followBack=[')
			for i in self.inst.tasks:
				if self.inst.followBack[i]==set():
					if i == self.inst.tasks[-1]:
						f.write('{}')
					else:
						f.write('{}, ')
				else:
					if i == self.inst.tasks[-1]:
						f.write('{}'.format(set([ j+1 
													for j in self.inst.followBack[i] ])))
					else:
						f.write('{}, '.format(set([ j+1 
													for j in self.inst.followBack[i] ])))
			f.write('];\n\n')

			f.write('precedeForw=[')
			for i in self.inst.tasks:
				if self.inst.precedeForw[i]==set():
					if i == self.inst.tasks[-1]:
						f.write('{}')
					else:
						f.write('{}, ')
				else:
					if i == self.inst.tasks[-1]:
						f.write('{}'.format(set([ j+1
													for j in self.inst.precedeForw[i] ])))
					else:
						f.write('{}, '.format(set([ j+1 
													for j in self.inst.precedeForw[i] ])))
			f.write('];\n\n')

			f.write('precedeBack=[')
			for i in self.inst.tasks:
				if self.inst.precedeBack[i]==set():
					if i == self.inst.tasks[-1]:
						f.write('{}')
					else:
						f.write('{}, ')
				else:
					if i == self.inst.tasks[-1]:
						f.write('{}'.format(set([ j+1
													for j in self.inst.precedeBack[i] ])))
					else:
						f.write('{}, '.format(set([ j+1
													for j in self.inst.precedeBack[i] ])))

			f.write('];\n\n')

			# pdb.set_trace()

	def OLD_write_dzn_sub_problem_file(self):
		# create datazinc file for cp solver
		with open(self.dznFile+'.dzn', 'w') as f:
			sortedTasks = sorted(list(self.tasks))
			# store maximum load value
			f.write('% max load\nmaxLoad={};\n\n'.format(self.bestCycleTimeUB))

			# store total number of tasks and number of assigned tasks
			f.write('% number of assigned tasks\nnTasks={};\n\n'.format(len(self.tasks)))

			# store list of assigned tasks
			f.write('% assigned tasks\nTASK={};\n\n'.format(set(range(1, len(self.tasks)+1))))

			# store processing times
			f.write('% processing times\ndur={};\n\n'.format([ self.inst.procList[i] for i in self.tasks ]))

			# store all precedence relations
			f.write('% precedence relations\nnPrecs={};\nsuc = ['.format(len([ (i,j) for (i,j) in self.inst.precList
																					 if i in self.tasks
																					 if j in self.tasks ])))
			# pdb.set_trace()
			for i in sortedTasks:
				if set([ index+1 for index,j in enumerate(sortedTasks)
								if j in self.inst.altPrecList[i]
								]) == set():
					if i == sortedTasks[-1]:
						f.write('{}')
					else:
						f.write('{}, ')
				else:
					if self.inst.altPrecList[i]==set():
						if i == sortedTasks[-1]:
							f.write('{}')
						else:
							f.write('{}, ')
					else:
						if i == sortedTasks[-1]:
							f.write('{}')
						else:
							f.write('{}, '.format(set([ index+1 for index,j in enumerate(sortedTasks)
																if j in self.inst.altPrecList[i]
															])))

			f.write('];\n\n')

			# store forward setup times
			assForwSU = self.inst.forwSU[list(self.tasks),:][:,list(self.tasks)]
			f.write('% forward setup times\nforwSU=[')
			for index,i in enumerate(self.tasks):
				if i == list(self.tasks)[-1]:
					f.write('| {} |];\n\n'.format(str(list(assForwSU[index]))[1:-1]))
				else:
					f.write('| {}\n\t\t'.format(str(list(assForwSU[index]))[1:-1]))

			# store backward setup times
			assBackSU = self.inst.backSU[list(self.tasks),:][:,list(self.tasks)]
			f.write('% backward setup times\nbackSU=[')
			for index,i in enumerate(self.tasks):
				if i == list(self.tasks)[-1]:
					f.write('| {} |];\n\n'.format(str(list(assBackSU[index]))[1:-1]))
				else:
					f.write('| {}\n\t\t'.format(str(list(assBackSU[index]))[1:-1]))

			# pdb.set_trace()
			# store the set of forward followers
			f.write('% other sets to define\nfollowForw=[')
			for i in sortedTasks:
				if self.inst.followForw[i]==set():
					if i == sortedTasks[-1]:
						f.write('{}')
					else:
						f.write('{}, ')
				else:
					if i == sortedTasks[-1]:
						f.write('{}')
					else:
						f.write('{}, '.format(set([ index+1 
													for index,j in enumerate(sortedTasks)
													if j in self.inst.followForw[i]
													])))
			f.write('];\n\n')

			f.write('followBack=[')
			for i in sortedTasks:
				if self.inst.followBack[i]==set():
					if i == list(self.inst.tasks)[-1]:
						f.write('{}')
					else:
						f.write('{}, ')
				else:
					if i == sortedTasks[-1]:
						f.write('{}'.format(set([ index+1 
													for index,j in enumerate(sortedTasks)
													if j in self.inst.followBack[i]
													])))
					else:
						f.write('{}, '.format(set([ index+1 
													for index,j in enumerate(sortedTasks)
													if j in self.inst.followBack[i]
													])))
			f.write('];\n\n')
			# pdb.set_trace()

			f.write('precedeForw=[')
			for i in sortedTasks:
				if set([ index+1 
							for index,j in enumerate(self.tasks)
							if j in self.inst.precedeForw[i]
							]) == set():
					if i == sortedTasks[-1]:
						f.write('{}')
					else:
						f.write('{}, ')
				else:
					if self.inst.precedeForw[i]==set():
						if i == sortedTasks[-1]:
							f.write('{}')
						else:
							f.write('{}, ')
					else:
						if i == sortedTasks[-1]:
							f.write('{}'.format(set([ index+1 
														for index,j in enumerate(sortedTasks)
														if j in self.inst.precedeForw[i]
														])))
						else:
							f.write('{}, '.format(set([ index+1 
														for index,j in enumerate(sortedTasks)
														if j in self.inst.precedeForw[i]
														])))
			f.write('];\n\n')

			f.write('precedeBack=[')
			for i in sortedTasks:
				if self.inst.precedeBack[i]==set():
					if i == sortedTasks[-1]:
						f.write('{}')
					else:
						f.write('{}, ')
				else:
					if i == sortedTasks[-1]:
						f.write('{}'.format(set([ index+1 
													for index,j in enumerate(sortedTasks)
													if j in self.inst.precedeBack[i]
													])))
					else:
						f.write('{}, '.format(set([ index+1 
													for index,j in enumerate(sortedTasks)
													if j in self.inst.precedeBack[i]
													])))

			f.write('];\n\n')

		# pdb.set_trace()

	def write_dzn_sub_problem_file(self):
		# create datazinc file for cp solver
		with open(self.dznFile+'.dzn', 'w') as f:
			# store the appropriate upper bound on the station load if using logic cuts
			if USE_LOGIC_CUTS:
				self.maxLoad = min(self.bestCycleTimeUB, self.naiveLoadUB)
			else:
				# self.maxLoad = self.inst.maxCycleTime
				# ^^^ using max cycle time is ridiculous and leads to terrible space
				# complexity issues for each sub-problem.
				# Calculate a new upper bound:
				self.maxLoad = self.naiveLoadUB

			# store maximum load value
			f.write('% max load\nmaxLoad={};\n\n'.format(self.maxLoad))

			# store total number of tasks and number of assigned tasks
			f.write('% number of assigned tasks\nnTasks={};\n\n'.format(len(self.tasks)))

			# store list of assigned tasks
			f.write('% assigned tasks\nTASK={};\n\n'.format(set(range(1, len(self.tasks)+1))))

			# store processing times
			f.write('% processing times\ndur={};\n\n'.format([ self.inst.procList[i] for i in self.oldIndexedTasks ]))

			# store all precedence relations
			f.write('% precedence relations\nnPrecs={};\nsuc = ['.format(len(self.precList)))
			# pdb.set_trace()
			for i in self.reindexedTasks:
				if self.altPrecList[i-1]==set():
					if i == self.reindexedTasks[-1]:
						f.write('{}')
					else:
						f.write('{}, ')
				else:
					f.write('{}, '.format(set([ j for j in self.altPrecList[i-1] ]))) 
			f.write('];\n\n')

			# store forward setup times
			assForwSU = self.inst.forwSU[self.oldIndexedTasks,:][:,self.oldIndexedTasks]
			f.write('% forward setup times\nforwSU=[')
			for index,i in enumerate(self.tasks):
				if i == list(self.tasks)[-1]:
					f.write('| {} |];\n\n'.format(str(list(assForwSU[index]))[1:-1]))
				else:
					f.write('| {}\n\t\t'.format(str(list(assForwSU[index]))[1:-1]))

			# store backward setup times
			assBackSU = self.inst.backSU[self.oldIndexedTasks,:][:,self.oldIndexedTasks]
			f.write('% backward setup times\nbackSU=[')
			for index,i in enumerate(self.tasks):
				if i == list(self.tasks)[-1]:
					f.write('| {} |];\n\n'.format(str(list(assBackSU[index]))[1:-1]))
				else:
					f.write('| {}\n\t\t'.format(str(list(assBackSU[index]))[1:-1]))

			# pdb.set_trace()
			# store the set of forward followers
			f.write('% other sets to define\nfollowForw=[')
			for i in self.reindexedTasks:
				if self.stationFollowForw[i-1]==set():
					if i == self.reindexedTasks[-1]:
						f.write('{}')
					else:
						f.write('{}, ')
				else:
					if i == self.reindexedTasks[-1]:
						f.write('{}'.format(set([ j
												for j in self.stationFollowForw[i-1] ])))
					else:
						f.write('{}, '.format(set([ j
												for j in self.stationFollowForw[i-1] ])))

			f.write('];\n\n')

			f.write('followBack=[')
			for i in self.reindexedTasks:
				if self.stationFollowBack[i-1]==set():
					if i == self.reindexedTasks[-1]:
						f.write('{}')
					else:
						f.write('{}, ')
				else:
					if i == self.reindexedTasks[-1]:
						f.write('{}'.format(set([ j
												for j in self.stationFollowBack[i-1] ])))
					else:
						f.write('{}, '.format(set([ j
												for j in self.stationFollowBack[i-1] ])))
			f.write('];\n\n')

			f.write('precedeForw=[')
			for i in self.reindexedTasks:
				if self.stationPrecedeForw[i-1]==set():
					if i == self.reindexedTasks[-1]:
						f.write('{}')
					else:
						f.write('{}, ')
				else:
					if i == self.reindexedTasks[-1]:
						f.write('{}'.format(set([ j
												for j in self.stationPrecedeForw[i-1] ])))
					else:
						f.write('{}, '.format(set([ j
												for j in self.stationPrecedeForw[i-1] ])))
			f.write('];\n\n')

			f.write('precedeBack=[')
			for i in self.reindexedTasks:
				if self.stationPrecedeBack[i-1]==set():
					if i == self.reindexedTasks[-1]:
						f.write('{}')
					else:
						f.write('{}, ')
				else:
					if i == self.reindexedTasks[-1]:
						f.write('{}'.format(set([ j
												for j in self.stationPrecedeBack[i-1] ])))
					else:
						f.write('{}, '.format(set([ j
												for j in self.stationPrecedeBack[i-1] ])))
			f.write('];\n\n')

		# pdb.set_trace()

	def solve_CP(self, time_remaining):
		# call the CP model from  the command line
		os.system("{0}fzn-chuffed {1}.fzn --time-out {4} -f --verbosity 2 2> {2}.txt "
				"| solns2out --output-time -o {3}.txt {1}.ozn".format(CHUFFED_DIR,
																		self.modelFile,
																		self.statsFile,
																		self.solFile,
																		round(time_remaining)));

	def store_station_solution_CP(self):
		with open(self.solFile+'.txt','r') as f:
			solOutput = f.readlines()

		solOutput = [ x.strip()
						for x in solOutput ]

		# find status of station sub-problem
		if '==========' in solOutput:
			self.status = 1
				# find output value of station load
			loadLine = [ elem
							for elem in solOutput
							if 'load' in elem ][0]
			self.stationLoad = int( re.findall(r'\d+', loadLine)[0] )
			# find start times of the assigned tasks
			startLine = [ elem 
							for elem in solOutput
							if 'start' in elem ][0]
			self.startTimesList = ast.literal_eval(startLine[8:])
		elif '% Time limit exceeded!' in solOutput:
			self.status = 'timeout'
		else:
			self.status = 0

	def store_station_statistics_CP(self):
		with open(self.statsFile+'.txt','r') as f:
			statsOutput = f.readlines()

		statsOutput = [ x.strip() 
						for x in statsOutput ]
		# pdb.set_trace()
		# find the number of nodes explored
		nodeLine = [ elem 
					for elem in statsOutput 
					if 'node' in elem ]
		if nodeLine != []:
			nodeString = nodeLine[0]
			self.nodesExplored = int( re.findall(r'\d+', nodeString)[0] )
		else:
			self.nodesExplored = 0

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
	# TSP SUB-PROBLEM
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

	def initialise_TSP(self):
		return True

	def solve_TSP(self):
		return True

	def get_station_load_TSP(self):
		return True

	def get_start_times_TSP(self):
		return True

	def get_status_TSP(self):
		return True

# Solver used to iterate through Benders for a given isntance
class Solver:
	def __init__(self, inst):
		self.inst = inst
		self.optimisation_times = []
		self.master_times = []
		self.sequencing_solve_times = []
		self.sequencing_overhead_times = []
		self.all_solutions_ever = [ [] for k in self.inst.stations]
		self.initialise_statistics()
		self.initialise_cut_sets()
		self.initialise()
		self.bigM = self.inst.maxCycleTime
		self.bestCycleTimeUB = self.inst.maxCycleTime
		self.bestCycleTimeLB = self.inst.minCycleTime

	def initialise(self):
		# define Gurobi model for the master
		self.model = Model('assemblyline')
		if args.quiet:
			self.model.setParam('OutputFlag', 0)
		if not args.very_quiet:
			print('Initialising the master problem... ', end='', flush=True)
		# time the initialisation of the master
		start = time.time()
		# initialise the master problem
		if RMP_TYPE == 'sched':
			self.init_vars()
			self.create_objective()
			self.create_constraints()
		elif RMP_TYPE == 'ass':
			self.init_vars_ass()
			self.create_objective_ass()
			self.create_constraints_ass()
		self.init_time = round(time.time() - start,4)
		self.optimisation_times.append(self.init_time)
		if not args.very_quiet:
			print('complete ({:.3f}s).'.format(self.init_time))
		# # set PreCrush to 1 to we can add lazy constaints in callbacks
		# self.model.setParam('PreCrush', 1)
		# set LazyConstraints to 1 to we can add lazy constaints in callbacks
		self.model.setParam('LazyConstraints', 1)

	def initialise_cut_sets(self):
		if USE_NOGOODS:
			self.numNoGoods = 0
			self.noGoods = []
		else:
			self.numNoGoods = '-'
		if USE_LOGIC_CUTS:
			self.numLogicCuts = 0
			self.logicCuts = []
		else:
			self.numLogicCuts = '-'
		if USE_INFER_CUTS:
			self.numInfAssCutsSimple = 0
			self.infAssCutsSimple = []
		else:
			self.numInfAssCutsSimple = '-'
		if USE_INFER_CUTS_SMART:
			self.numInfAssCutsSmart = 0
			self.infAssCutsSmart = []
		else:
			self.numInfAssCutsSmart = '-'
		if USE_INFER_CUTS_SMARTEST:
			self.numInfAssCutsSmartest = 0
			self.infAssCutsSmartest = []
		else:
			self.numInfAssCutsSmartest = '-'
		if USE_GLOBAL_BOUNDS:
			self.numGlobalUB = 0
			self.globalUB = []
			self.numGlobalLB = 0
			self.globalLB = []
		else:
			self.numGlobalUB = '-'
			self.numGlobalLB = '-'
		self.numTotalCuts = 0

	def initialise_statistics(self):
		# node statistics
		self.statsMasterNodes = np.empty([0],dtype=int)
		self.statsTotalMasterNodes = 0
		self.statsFinalMasterNodes = 0

		self.statsSubProblemNodes = [np.empty([0],dtype=int)]
		self.statsTotalSubProblemNodes = 0
		self.statsFinalSubProblemNodes = 0
		self.statsTotalNodes = 0

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
	# MASTER DEFINITION #1 (relaxed scheduling component)
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
	def init_vars(self):
		self.cycleTime = self.model.addVar(lb=self.inst.minCycleTime, 
										   ub=self.inst.maxCycleTime,
										   obj=0.0,
										   vtype=GRB.INTEGER,
										   name='c')

		# Initialise x variables: Station assignment
		self.xs = self.model.addVars([ (i,k)
									   for i in self.inst.tasks
									   for k in self.inst.feasibleStations[i] ],
									 vtype=GRB.BINARY, name='x')

		# Initialise y variables: Forward Sequencing
		self.ys = self.model.addVars([ (i,j,k)
									   for i in self.inst.tasks
									   for j in self.inst.followForw[i]
									   for k in self.inst.feasibleStations[i] ],
									 vtype=GRB.BINARY, name='y')

		# Initialise z variables: Backward Sequencing
		self.zs = self.model.addVars([ (i,j,k)
									   for i in self.inst.tasks
									   for j in self.inst.followBack[i]
									   for k in self.inst.feasibleStations[i] ],
									 vtype=GRB.BINARY, name='z')

		# Initialise xi variables: Sub-sequence Setup Lower Bound
		self.xis = self.model.addVars(self.inst.stations,
									  ub=self.inst.maxCycleTime,
									  vtype=GRB.INTEGER,
									  name='xi')

		self.model.update()

	def reinitialise_master(self, time_remaining):
		if not WARM_START:
			# remove previous RMP solution as the warm-starting solution
			self.model.reset()
		# limit the current relaxed master's runtime
		self.model.setParam('TimeLimit', time_remaining)

	def create_objective(self):
		# create the objective of the master problem for this iteration
		self.objective = self.cycleTime
		self.model.setObjective(self.objective, GRB.MINIMIZE)

	def create_constraints(self):
		# Each task i is assigned exactly one station
		self.model.addConstrs(( self.xs.sum(i,'*') == 1
								for i in self.inst.tasks), 'consOneStationPerTask')

		# Each task has exactly one successor (in the forward and backward loads together)
		self.model.addConstrs(( self.ys.sum(i,'*',k) + self.zs.sum(i,'*',k) == self.xs[i,k]
								for i in self.inst.tasks
								for k in self.inst.feasibleStations[i]), 'consOneSuccessor')

		# Each task has exactly one predecessor (in the forw and back loads together)
		self.model.addConstrs((   sum([ self.ys[i,j,k] for i in self.inst.precedeForw[j] ]) 
								+ sum([ self.zs[i,j,k] for i in self.inst.precedeBack[j] ]) == self.xs[j,k]
								for j in self.inst.tasks
								for k in self.inst.feasibleStations[j] ), 'consOnePredecessor')

		# Each station has at least one backward setup (relaxation: can have more than k backward setups)
		self.model.addConstrs(( sum([ self.zs.sum(i,'*',k)
									  for i in self.inst.feasibleTasks[k] ]) >= 1
								for k in self.inst.stations ), 'consAtLeastOneBackwardSU')

		# Precedence relations are respected between stations (relaxation: not necessarily within stations)
		self.model.addConstrs((    sum([ k*self.xs[i,k] for k in self.inst.feasibleStations[i] ])
								<= sum([ k*self.xs[j,k] for k in self.inst.feasibleStations[j] ])
								for (i,j) in self.inst.precList), 'precedenceRelations')

		# Each station load respects the cycle time (relaxation: relaxed setup cost)
		self.model.addConstrs((   sum([ self.inst.procList[i]*self.xs[i,k]
										for i in self.inst.tasks ])
								+ self.xis[k] <= self.cycleTime # 
								for k in self.inst.stations), 'consCycleTimeGreaterThanStationLoads')

		# Fix relaxed setup time for each station
		self.model.addConstrs(( self.xis[k] ==
								  sum([  sum([ self.inst.forwSU[i][j]*self.ys[i,j,k]
											 for j in self.inst.followForw[i] ])
									   + sum([ self.inst.backSU[i][j]*self.zs[i,j,k]
											 for j in self.inst.followBack[i] ])
									   for i in self.inst.feasibleTasks[k] ])
								for k in self.inst.stations), 'consFixRelaxedSetupTime')

		# Bounds for the cycle time
		self.consUB = self.model.addConstr(self.cycleTime <= self.inst.maxCycleTime, 'consCycleTimeUB')
		self.consLB = self.model.addConstr(self.cycleTime >= self.inst.minCycleTime, 'consCycleTimeLB')

	def optimise(self):
		# find optimal solution to current relaxed master (maybe use callbacks)
		if SIZE_OF_SUB_TOURS_TO_ELIMNATE > 0:
			# first add needed information to the model object
			self.model._tasks = self.inst.tasks
			self.model._numStations = self.inst.numStations
			self.model._xs = self.xs
			self.model._ys = self.ys
			self.model._zs = self.zs
			# optimise the master using callbacks
			self.model.optimize(callback_sub_tour_elimination)
		else:
			self.model.optimize()

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
	# MASTER DEFINITION #2 (strictly assignment problem)
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
	def init_vars_ass(self):
		self.cycleTime = self.model.addVar(lb=self.inst.minCycleTime, 
										   ub=self.inst.maxCycleTime,
										   obj=0.0,
										   vtype=GRB.INTEGER,
										   name='c')

		# Initialise x variables: Station assignment
		self.xs = self.model.addVars([ (i,k)
									   for i in self.inst.tasks
									   for k in self.inst.feasibleStations[i] ],
									 vtype=GRB.BINARY, name='x')

		self.model.update()

	def reinitialise_master_ass(self, time_remaining):
		# if not WARM_START:
		# 	# remove previous RMP solution as the warm-starting solution
		# 	self.model.reset()
		# limit the current relaxed master's runtime
		self.model.setParam('TimeLimit', time_remaining)

	def create_objective_ass(self):
		# create the objective of the master problem for this iteration
		self.objective = self.cycleTime
		self.model.setObjective(self.objective, GRB.MINIMIZE)

	def create_constraints_ass(self):
		# Each task i is assigned exactly one station
		self.model.addConstrs(( self.xs.sum(i,'*') == 1
								for i in self.inst.tasks), 'consOneStationPerTask')

		# Precedence relations are respected between stations (relaxation: not necessarily within stations)
		self.model.addConstrs((    sum([ k*self.xs[i,k] for k in self.inst.feasibleStations[i] ])
								<= sum([ k*self.xs[j,k] for k in self.inst.feasibleStations[j] ])
								for (i,j) in self.inst.precList), 'precedenceRelations')

		# Each station load respects the cycle time (relaxation: relaxed setup cost)
		self.model.addConstrs((   sum([ self.inst.procList[i]*self.xs[i,k]
										for i in self.inst.tasks ])
								<= self.cycleTime
								for k in self.inst.stations), 'consCycleTimeGreaterThanStationLoads')

		# NEW! (2017-06-17) Each station must be assigned at least 1 task
		self.model.addConstrs((  self.xs.sum('*',k) >= 1
								for k in self.inst.stations), 'minOneTaskAtEachStation')

		# Bounds for the cycle time
		self.consUB = self.model.addConstr(self.cycleTime <= self.inst.maxCycleTime, 'consCycleTimeUB')
		self.consLB = self.model.addConstr(self.cycleTime >= self.inst.minCycleTime, 'consCycleTimeLB')

	def optimise_ass(self):
		# find optimal solution to current relaxed master (maybe use callbacks)
		if SIZE_OF_SUB_TOURS_TO_ELIMNATE > 0:
			# first add needed information to the model object
			self.model._tasks = self.inst.tasks
			self.model._numStations = self.inst.numStations
			self.model._xs = self.xs
			self.model._ys = self.ys
			self.model._zs = self.zs
			# optimise the master using callbacks
			self.model.optimize(callback_sub_tour_elimination)
		else:
			self.model.optimize(master_callback)

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
	# BENDERS CUTS DEFINITION
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
	def add_nogood_cut(self, assignment):
		# nogood cut
		# for the current FULL assignment and cycle time, this
		#	assignment solution is now removed
		self.curCycleTime = round(self.curCycleTime)
		if not args.very_quiet:
			print('\n CUT: Nogood cut #{} added:'.format(self.numNoGoods))
			print('   [c = {}, assignment = ...]'.format(self.curCycleTime))
		# pdb.set_trace()
		self.model.addConstr(self.curCycleTime + 1 - sum([sum([ (1 - self.xs[i,k])
																for i in self.taskAssignment[k]])
																for k in self.inst.stations])
								<= self.cycleTime,
							 'NoGoodCut[{}]'.format(self.numNoGoods))
		self.numNoGoods += 1
		self.noGoods.append({'cycleTime':self.curCycleTime})

	def add_logic_cut_infeasible_assignment(self, k, tasks):
		# logic cut
		count = len(tasks)
		if not args.very_quiet:
			print('\n   CUT: Logic cut #{} added (infeas. assignment):'.format(self.numLogicCuts))
			print('  \t[Station {}, tasks = {}]'.format(k,tasks))
		self.model.cbLazy(sum([ (1 - self.xs[i,k]) for i in tasks ]) >= 1,
							 'LogicCut[{}]'.format(self.numLogicCuts))
		self.logicCuts.append({'stationNum':k, 
							   'tasks':tasks, 
							   'cycleTime':self.bestCycleTimeUB})
		self.numLogicCuts += 1

	def add_infer_cut_infeasible_assignment_simple(self, k, tasks):
		# infer cut
		count = len(tasks)
		if not args.very_quiet:
			print('   CUT: Infer cut #{} added (simple):'.format(self.numInfAssCutsSimple))
			print('  \t[Station {}: tasks = {} implies c >= {}]'.format(k,tasks,round(self.curStationLoad[k])))
		self.model.cbLazy(self.cycleTime >=   self.curStationLoad[k] 
											   - self.bigM*(count - sum([ self.xs[i,k] for i in tasks ])),
							 'InferCut[{}]'.format(self.numInfAssCutsSimple))
		self.infAssCutsSimple.append({'stationNum':k, 
									  'tasks':tasks, 
									  'cycleTime':self.curStationLoad[k]})
		self.numInfAssCutsSimple += 1

	def add_infer_cut_infeasible_assignment_smart(self, k, tasks):
		# infer cut
		count = len(tasks)
		maxSetup = [None for i in self.inst.tasks]
		minSetup = [None for i in self.inst.tasks]
		burdenUB = [None for i in self.inst.tasks]
		minLinkingSetup = [None for i in self.inst.tasks]
		for i in tasks:
			# calculate upper and lower bounds on setup time
			maxSetup[i] =	max( 
								[self.inst.forwSU[j][i] 
									for j in self.inst.precedeForw[i] ] + 
								[self.inst.backSU[j][i]
									for j in self.inst.precedeBack[i] ] 
						) + max( 
								[self.inst.forwSU[i][j]
									for j in self.inst.followForw[i] ] +
								[self.inst.backSU[i][j]
									for j in self.inst.followBack[i] ] 
						)

			# broken version that provides a *too* tight upper bound
			# maxSetup[i] = 	max( 
			# 					[self.inst.forwSU[j][i] 
			# 						for j in list(set(self.inst.precedeForw[i]).intersection(tasks)) ] + 
			# 					[self.inst.backSU[j][i]
			# 						for j in list(set(self.inst.precedeBack[i]).intersection(tasks)) ] 
			# 			) + max( 
			# 					[self.inst.forwSU[i][j]
			# 						for j in list(set(self.inst.followForw[i]).intersection(tasks)) ] +
			# 					[self.inst.backSU[i][j]
			# 						for j in list(set(self.inst.followBack[i]).intersection(tasks)) ] 
			# 			)

			# calculate the minimum setup that can occur between the tasks on station k if
			# task i is removed from the sequence. Call this the 'linking' setup
			remainingTasks = set(self.inst.tasks).difference(set([i]))

			try:
				# to prevent the min setup being 0 we remove the backward setup from i to i from consideration
				# this breaks/excludes the trivial case when only one task is assigned a station
				minLinkingSetup[i] = min(
									[self.inst.forwSU[idash][j]
										for idash in list(remainingTasks)
										for j in list(set(self.inst.followForw[idash]).intersection(remainingTasks)) ] +
									[self.inst.backSU[idash][j]
										for idash in list(remainingTasks)
										for j in list(set(self.inst.followBack[idash]).intersection(remainingTasks))
										if j!=idash ]
								)

				# broken version which provides a *too* tight lower bound
				# minLinkingSetup[i] = min(
				# 					[self.inst.forwSU[idash][j]
				# 						for idash in list(remainingTasks)
				# 						for j in list(set(self.inst.followForw[idash]).intersection(remainingTasks)) ] +
				# 					[self.inst.backSU[idash][j]
				# 						for idash in list(remainingTasks)
				# 						for j in list(set(self.inst.followBack[idash]).intersection(remainingTasks))
				# 						if j!=idash ]
				# 				)
			except ValueError:
				# we handle the trival case poorly so when this arises,
				# set the min cost to 0 (this is theoretically correct)
				minLinkingSetup[i] = 0

			# for this cut upper bound on burden is "duration + maxSU - minSU"
			burdenUB[i] = self.inst.procList[i] + maxSetup[i] - minLinkingSetup[i]

		# pdb.set_trace()
		if not args.very_quiet:
			print('   CUT: Infer cut #{} added (smart):'.format(self.numInfAssCutsSmart))
			print('  \t[Station {}: tasks = {} implies c >= {}]'.format(k,tasks,round(self.curStationLoad[k])))
		self.model.cbLazy(self.cycleTime >=   self.curStationLoad[k] 
											   - sum([ burdenUB[i]*(1 - self.xs[i,k])
														for i in tasks ]),
							 'InferCut2[{}]'.format(self.numInfAssCutsSmart))
		self.infAssCutsSmart.append({'stationNum':k, 
									  'tasks':tasks, 
									  'cycleTime':self.curStationLoad[k]})
		self.numInfAssCutsSmart += 1

	def add_infer_cut_infeasible_assignment_smartest(self, k, tasks):
		# infer cut
		count = len(tasks)
		# setupCosts = [[0] for i in self.inst.tasks]
		maxSetup = [None for i in self.inst.tasks]
		minSetup = [None for i in self.inst.tasks]
		burdenUB = [None for i in self.inst.tasks]
		burdenLB = [None for i in self.inst.tasks]
		otherTasks = set(self.inst.tasks).difference(set(tasks))
		for i in tasks:
			# calculate upper bounds on setup time
			maxSetup[i] =	max( 
								[self.inst.forwSU[j][i] 
									for j in self.inst.precedeForw[i] ] + 
								[self.inst.backSU[j][i]
									for j in self.inst.precedeBack[i] ] 
						) + max( 
								[self.inst.forwSU[i][j]
									for j in self.inst.followForw[i] ] +
								[self.inst.backSU[i][j]
									for j in self.inst.followBack[i] ] 
						)
			burdenUB[i] = self.inst.procList[i] + maxSetup[i]

		for i in list(otherTasks):
			try:
				# to prevent the min setup being 0 we remove the backward setup from i to i from consideration
				# this breaks/excludes the trivial case when only one task is assigned a station
				minSetup[i] =	min( 
									[self.inst.forwSU[j][i] 
										for j in self.inst.precedeForw[i] ] + 
									[self.inst.backSU[j][i]
										for j in self.inst.precedeBack[i]
										if j!=i ] 
							) + min( 
									[self.inst.forwSU[i][j]
										for j in self.inst.followForw[i] ] +
									[self.inst.backSU[i][j]
										for j in self.inst.followBack[i]
										if j!=i ] 
							)

				# broken version which provides a *too* tight lower bound
				# minLinkingSetup[i] = min(
				# 					[self.inst.forwSU[idash][j]
				# 						for idash in list(remainingTasks)
				# 						for j in list(set(self.inst.followForw[idash]).intersection(remainingTasks)) ] +
				# 					[self.inst.backSU[idash][j]
				# 						for idash in list(remainingTasks)
				# 						for j in list(set(self.inst.followBack[idash]).intersection(remainingTasks))
				# 						if j!=idash ]
				# 				)
			except ValueError:
				# we handle the trival case poorly so when this arises,
				# set the min cost to 0 (this is theoretically correct)
				minSetup[i] = 0

			burdenLB[i] = self.inst.procList[i] + minSetup[i]

		# pdb.set_trace()
		if not args.very_quiet:
			print('   CUT: Infer cut #{} added (smartest):'.format(self.numInfAssCutsSmartest))
			print('  \t[Station {}: tasks = {} implies c >= {}]'.format(k,tasks,round(self.curStationLoad[k])))
		self.model.cbLazy(self.cycleTime >=   self.curStationLoad[k] 
											   - sum([ burdenUB[i]*(1 - self.xs[i,k])
											   			for i in tasks ])
											   + sum([ burdenLB[i]*self.xs[i,k]
											   			for i in otherTasks ]),
							 'InferCut3[{}]'.format(self.numInfAssCutsSmartest))
		self.infAssCutsSmartest.append({'stationNum':k, 
										'tasks':tasks, 
										'cycleTime':self.curStationLoad[k]})
		self.numInfAssCutsSmartest += 1

	def add_global_bounds(self, allowGlobalUB):
		# method to add all global bounds after all sub-problems have completed

		if allowGlobalUB:
			# upper bound
			cycleTimeUB = round(max([ self.curStationLoad[k] for k in self.inst.stations ]))
			# pdb.set_trace()
			if cycleTimeUB < self.bestCycleTimeUB:
				self.bestCycleTimeUB = cycleTimeUB
				self.add_global_upper_bound()

		# !~~~~ CURRENTLY NOT INCLUDING LOWER BOUND CUTS ~~~~~~~!
		# lower bound
		# cycleTimeLB = round(self.model.objval, 2)
		# if cycleTimeLB > self.bestCycleTimeLB:
		# 	# update best lower bound and add a cut
		# 	self.bestCycleTimeLB = cycleTimeLB
		# 	self.add_global_lower_bound()
		# update the model with the changed constraints
		# self.model.update()

	def add_global_upper_bound(self):
		# global bound
		if not args.very_quiet:
			print('\n BOUND: Global UB #{} added: [c <= {}]'.format(self.numGlobalUB,
																	self.bestCycleTimeUB))
		# change the rhs
		# self.consUB.setAttr('rhs', self.bestCycleTimeUB)
		self.model.cbLazy(self.cycleTime <= self.bestCycleTimeUB, 'consCycleTimeUB')

		self.globalUB.append(self.bestCycleTimeUB)
		self.numGlobalUB += 1

		# store current Benders iteration for future referral
		self.mostRecentUpperBoundIter = self.bendersIter

	def add_global_lower_bound(self):
		# global bound
		if not args.very_quiet:
			print('\n BOUND: Global LB #{} added: [c >= {}]'.format(self.numGlobalLB,
																	self.bestCycleTimeLB))
		# change the rhs
		self.consLB.setAttr('rhs', self.bestCycleTimeLB)

		self.globalLB.append(self.bestCycleTimeLB)
		self.numGlobalLB += 1

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
	# BENDERS IMPLEMENTATION
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
	def benders_optimise_with_optimality_sub_problems(self, benders_gap=0.01):
		startBenders = time.time()
		self.time_limit_exceeded = False
		self.master_timed_out = False
		self.bendersIter = 0
		doneBenders = False
		self.gap = []
		self.stations = [None for k in self.inst.stations]
		self.taskAssignment = [None for k in self.inst.stations]
		self.stationSatisfiesCurCycleTime = [True for k in self.inst.stations]
		self.stationFeasible = [True for k in self.inst.stations]
		self.curStationLoad = [None for k in self.inst.stations]
		self.startTimes = [None for k in self.inst.stations]

		while not doneBenders:
			# early termination consitions
			if self.bendersIter >= MAX_BENDERS_ITERATIONS:
				sys.exit('Terminating. Maximum number of Benders iterations exceeded.')

			# define the time used up until this relaxed master
			self.RMP_time_used = round(time.time()-startBenders,4)
			if self.RMP_time_used > TIMELIMIT:
				doneBenders = True
				self.time_limit_exceeded = True
				break
			# default to allowing a global UB this iteration
			allowGlobalUB = True
			if not args.very_quiet:
				print('\n{:.1f}/{} seconds elapsed'.format(self.RMP_time_used,TIMELIMIT))
				print('===============================')
				print('Master %d: ' %(self.bendersIter), end='', flush=True)

			# Master optimisation for current Benders iteration
			self.solve_master_problem(TIMELIMIT - self.RMP_time_used)

			self.gap.append(round(float((self.bestCycleTimeUB - self.curCycleTime)/self.curCycleTime)*100,4))

			if time.time() - startBenders > TIMELIMIT:
				self.master_timed_out = True
				doneBenders = True
				self.time_limit_exceeded = True
				break

			if not args.very_quiet:
				print('\tCycle: \t{}'.format(round(self.curCycleTime)))
				print('\t\tUB: \t{}'.format(self.bestCycleTimeUB))
				print('\t\tGap: \t{:.2f} %\n'.format(self.gap[self.bendersIter]))

			# !~~~~~ this should probably check if gap <= 0 as we might skip a cycle time value right?
			# if gap found is 0 then we have already found a feasible solution to the sub-problems
			if self.gap[self.bendersIter] == 0:
				# ignore the current master solution and take the old one instead
				# pdb.set_trace()
				for k in self.inst.stations:
					self.taskAssignment[k] = self.all_solutions_ever[k][self.mostRecentUpperBoundIter]['tasks']
					self.all_solutions_ever[k].append({'tasks': self.taskAssignment[k]})
			else:
				# store current assignment
				if self.bendersIter > 0:
					self.statsSubProblemNodes.append(np.empty([0],dtype=int))
				for k in self.inst.stations:
					self.taskAssignment[k] = { i for i in self.inst.tasks if self.xs[i,k].x > 0.5 }
					self.all_solutions_ever[k].append({'tasks': self.taskAssignment[k]})

			# solve each sub-problem, adding cuts to master
			for k in self.inst.stations:
				if not args.very_quiet:
					print(' Station %d' %(k), end='', flush=True)

				# define the time used up until starting this sub-problem
				self.SP_time_used = round(time.time() - startBenders,4)
				# check if we are out of time before starting each sub-problem
				if self.SP_time_used > TIMELIMIT:
					doneBenders = True
					self.time_limit_exceeded = True
					break

				# solve the current station's sub-problem
				result = self.solve_sub_problem(k)
				# check if time-limit is exceeded
				if self.time_limit_exceeded:
					doneBenders = True
					break

				# if we have already processed ths assignment before move onto next sub problem
				if result == True:
					allowGlobalUB = False
			# if exceeded time limit after sovling a sub-problem exit benders and output
			if self.time_limit_exceeded:
				break

			if not False in self.stationFeasible:
				self.mostRecentFeasibleCycleTime = round(max(self.curStationLoad),4)

			# stopping condition: continuing until all sub-problem solutions <= master solution
			if False in self.stationSatisfiesCurCycleTime or False in self.stationFeasible:
				# add cuts and iterate if we arentt done

				# ** don't do a nogood cut for optimality sub-problems **
				if USE_NOGOODS:
					self.add_nogood_cut(self.taskAssignment)
					# do not do any global upper bounds if we used a nogood
					allowGlobalUB = False

				if USE_GLOBAL_BOUNDS:
					self.add_global_bounds(allowGlobalUB)
				self.bendersIter += 1
			else:
				# self.debug_final_result()
				# pdb.set_trace()
				doneBenders = True

		self.benders_time = time.time() - startBenders
		# self.optimisation_times.append(benders_time)
		# record the total number of cuts

	def benders_optimise_with_feasibility_sub_problems(self, benders_gap=0.01):
		# TO BE COMPLETED
		print('USE_NOGOODS', USE_NOGOODS)
		self.bendersIter = 0
		doneBenders = False
		while not doneBenders:
			# solve master
			# blah blah
			self.curCycleTime = 0
			for k in self.inst.stations:
				# solve sub-problem
				# blah blah
				if USE_NOGOODS:
					self.add_nogood_cut(self.taskAssignment)
					# do not do any global upper bounds if we used a nogood
					allowGlobalUB = False

			doneBenders = True

	def benders_optimise_with_master_callbacks(self, benders_gap=0.01):
		self.startBenders = time.time()
		self.time_limit_exceeded = False
		self.master_timed_out = False
		self.bendersIter = 0
		self.doneBenders = False
		self.gap = []
		self.stations = [None for k in self.inst.stations]
		self.taskAssignment = [None for k in self.inst.stations]
		self.stationSatisfiesCurCycleTime = [True for k in self.inst.stations]
		self.stationFeasible = [True for k in self.inst.stations]
		self.curStationLoad = [None for k in self.inst.stations]
		self.startTimes = [None for k in self.inst.stations]

		self.allowGlobalUB = True

		# define the time used up until the first relaxed master
		self.RMP_time_used = round(time.time()-self.startBenders,4)
		
		if not args.very_quiet:
			print('\n{:.1f}/{} seconds elapsed'.format(self.RMP_time_used,TIMELIMIT))
			print('===============================')
			print('Master %d: ' %(self.bendersIter), end='', flush=True)

		# begin Benders iteration
		self.solve_master_problem_with_callbacks(TIMELIMIT - self.RMP_time_used)
		pdb.set_trace()

		if self.model.status == GRB.TIME_LIMIT:
			self.curCycleTime = self.inst.maxCycleTime

		self.statsMasterNodes = np.append(self.statsMasterNodes, int(self.model.nodecount))

		self.benders_time = time.time() - self.startBenders

	def solve_master_problem_with_callbacks(self, timeRemaining):
		# add updated timelimit
		self.reinitialise_master_ass(timeRemaining)

		self.startMaster = time.time()
		self.optimise_ass()
		
	def solve_master_problem(self, timeRemaining):
		# add updated timelimit
		self.reinitialise_master_ass(timeRemaining)

		startMaster = time.time()
		if RMP_TYPE == 'sched':
			self.optimise()
		elif RMP_TYPE == 'ass':
			self.optimise_ass()
		self.master_times.append(time.time() - startMaster)
		self.optimisation_times.append(self.master_times[-1])
		# if self.model.getAttr('status') != 2:
		# 	pdb.set_trace()
		# 	sys.exit('\nError: Master was not solved optimally')

		# if the master didnt time out then update current cycle time
		if self.model.status != 9:
			self.curCycleTime = round(self.model.objVal,4)
		else:
			self.curCycleTime = self.inst.maxCycleTime
		self.statsMasterNodes = np.append(self.statsMasterNodes, int(self.model.nodecount))

	def solve_sub_problem(self, k):
		# initialise sub-problem and solve it
		self.stations[k] = Station(self.inst, k, self.taskAssignment[k], self.curCycleTime, self.bestCycleTimeUB)
		# check if assignment is new, and don't solve if we already have
		[newAssignment, isFeasible] = self.is_assignment_new(k)
		logicallyInfeasibleAssignment = False

		# solve or retrieve old solution
		if newAssignment:
			# solve new sub-problem
			logicallyInfeasibleAssignment = self.call_sub_problem_solver(k)
			# check if time-limit is exceeded
			if self.time_limit_exceeded:
				return
			if USE_LOGIC_CUTS and logicallyInfeasibleAssignment == True:
				isFeasible = False
				satisfiesCurCycleTime = False
				# cut all solutions with the current assignment for this station
				self.add_logic_cut_infeasible_assignment(k,self.taskAssignment[k])
				self.store_sub_problem_result(k, satisfiesCurCycleTime, isFeasible)
				return logicallyInfeasibleAssignment
		else:
			# if old assignment was feasible then just copy the results
			if isFeasible:
				# print station load
				if not args.very_quiet:
					print(' Load = \t{}'.format(round(self.curStationLoad[k])))
				# check if old load satisfies current cycle time
				if self.curStationLoad[k] <= self.curCycleTime:
					satisfiesCurCycleTime = True
				else:
					# pdb.set_trace()
					satisfiesCurCycleTime = False
				self.store_sub_problem_result(k, satisfiesCurCycleTime, isFeasible)
				# return immediately, not doing any inference cuts
				return logicallyInfeasibleAssignment

		# print station load
		if not args.very_quiet and isFeasible:
			try:
				print(' Load = \t{}'.format(round(self.curStationLoad[k])))
			except TypeError:
				pdb.set_trace()

		approxStationLoad = round(self.curStationLoad[k], 2)
		approxCycleTime = round(self.curCycleTime, 2)
		# if station load is greater than master cycle time then this solution cannot be optimal
		if approxStationLoad > approxCycleTime:
			satisfiesCurCycleTime = False
			# add infer cuts for this station if we haven't already made a logic cut
			if USE_INFER_CUTS and logicallyInfeasibleAssignment != True:
				self.add_infer_cut_infeasible_assignment_simple(k, self.taskAssignment[k])
			if USE_INFER_CUTS_SMART and logicallyInfeasibleAssignment != True:
				self.add_infer_cut_infeasible_assignment_smart(k, self.taskAssignment[k])
			if USE_INFER_CUTS_SMARTEST and logicallyInfeasibleAssignment != True:
				self.add_infer_cut_infeasible_assignment_smartest(k, self.taskAssignment[k])

		else:
			satisfiesCurCycleTime = True
		self.store_sub_problem_result(k, satisfiesCurCycleTime, isFeasible)

		return logicallyInfeasibleAssignment

	def call_sub_problem_solver(self, k):
		logicallyInfeasibleAssignment = False
		if SUB_PROBLEM_SOLVER == 'mip':
			self.stations[k].initialise_MIP(TIMELIMIT-self.SP_time_used)
			# if args.quiet:
			# 	self.stations[k].model.setParam('OutputFlag', 0)
			startSequencing = time.time()
			self.stations[k].solve_MIP()
			# record station sub-problem optimisation time
			self.sequencing_solve_times.append(time.time() - startSequencing)
			self.optimisation_times.append(self.sequencing_solve_times[-1])
			self.sequencing_overhead_times.append(self.stations[k].init_time)
			self.optimisation_times.append(self.sequencing_overhead_times[-1])
			# if MIP solved optimally then store results otherwise return infeasible
			if self.stations[k].model.getAttr('Status') == 2:
				# optimal
				self.curStationLoad[k] = round(self.stations[k].model.objval,4)
				self.startTimes[k] = [ round(self.stations[k].ss[i].x) for i in self.taskAssignment[k]]
				self.statsSubProblemNodes[self.bendersIter] = np.append(self.statsSubProblemNodes[self.bendersIter],
																		int(self.stations[k].model.nodecount))
			elif self.stations[k].model.getAttr('Status') == 9:
				# exceeded time limit given to sub-problem
				self.time_limit_exceeded = True
				return
			else:
				# do a logic cut of the infeasible assignment
				logicallyInfeasibleAssignment = True

		elif 'cp' in SUB_PROBLEM_SOLVER:
			self.stations[k].initialise_CP()
			startSequencing = time.time()
			self.stations[k].solve_CP(TIMELIMIT-self.SP_time_used)
			# record station sub-problem optimisation time
			self.sequencing_solve_times.append(time.time() - startSequencing)
			self.optimisation_times.append(self.sequencing_solve_times[-1])
			self.sequencing_overhead_times.append(self.stations[k].init_time)
			self.optimisation_times.append(self.sequencing_overhead_times[-1])
			self.stations[k].store_station_solution_CP()
			# pdb.set_trace()
			self.stations[k].store_station_statistics_CP()
			# if CP solved optimally then store results otherwise return infeasible
			if self.stations[k].status == 1:
				self.curStationLoad[k] = self.stations[k].stationLoad
				# pdb.set_trace()
				self.startTimes[k] = [ self.stations[k].startTimesList[index] 
										for index,i in enumerate(self.taskAssignment[k]) ]
				self.statsSubProblemNodes[self.bendersIter] = np.append(self.statsSubProblemNodes[self.bendersIter],
																		self.stations[k].nodesExplored)
			elif self.stations[k].status == 'timeout':
				self.time_limit_exceeded = True
				return
			else:
				# do a logic cut of the infeasible assignment
				logicallyInfeasibleAssignment = True

		elif SUB_PROBLEM_SOLVER == 'tsp':
			self.stations[k].initialise_TSP()
			self.stations[k].solve_TSP()
			# if TSP solver found optimality then store results otherwise return infeasible
			if self.stations[k].getStatusTSP() == 'optimal':
				self.curStationLoad[k] = self.stations[k].getStationLoadTSP() # create this function
				self.startTimes[k] = self.stations[k].getStartTimesTSP() # create this function
			else:
				# do a logic cut of the infeasible assignment
				logicallyInfeasibleAssignment = True
		else:
			sys.exit('\n\nError: Typo in command line argument or sub-problem solver requested is not-supported.\n')

		return logicallyInfeasibleAssignment

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
	# SUPPORT FUNCTIONALITY
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
	def store_sub_problem_result(self, k, satisfiesCurCycleTime, feasible):
		self.stationSatisfiesCurCycleTime[k] = satisfiesCurCycleTime
		self.stationFeasible[k] = feasible
		self.all_solutions_ever[k][self.bendersIter].update( {'satisfiesCurCycleTime':satisfiesCurCycleTime,
															  'feasible':feasible,
															  'cycleTime':self.curStationLoad[k],
															  'startTimes':self.startTimes[k]} )

	def is_assignment_new(self, k):
		# check the current assigment of station k against all previous
		# if a previous assignment is identical then return False, o/w True
		new = True
		feasible = True
		for oldSolution in self.all_solutions_ever[k][:-1][::-1]:
			if oldSolution['tasks'] == self.taskAssignment[k]:
				new = False
				# this assignment is not new so it won't be solved, so we must store 
				# the results of the previous sub-problem solution.
				self.curStationLoad[k] = oldSolution['cycleTime']
				self.startTimes[k] = oldSolution['startTimes']
				feasible = oldSolution['feasible']
				break
		return new, feasible

	def OLD_debug_final_result(self):
		print('\n~~Debugging~~')
		print('sum(y):',sum([self.ys[i,j,k].x for i in self.inst.tasks for j in self.inst.followForw[i] for k in self.inst.stations]))
		print('sum(z):',sum([self.zs[i,j,k].x for i in self.inst.tasks for j in self.inst.followBack[i] for k in self.inst.stations]))		
		print('station loads:',[round(self.all_solutions_ever[k][self.bendersIter]['cycleTime']) for k in self.inst.stations])

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
	# OUTPUT METHODS
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
	def check_for_feasibility_and_optimality(self):
		# check if we have a feasible solution. i.e. if Benders has completed at least one full iteration
		self.solFeasible = 1
		if self.bendersIter == 0 and self.time_limit_exceeded:
			if self.model.solcount == 0 or None in self.curStationLoad:
				self.solFeasible = 0

		#define the status of the Benders algorithm
		if self.time_limit_exceeded:
			self.solOptimal = 0
		else:
			self.solOptimal = 1 # one means optimal

	def process_solution_statistics(self):
		if self.solFeasible == 1:
			if self.master_timed_out:
				# if the current master timeed out then take the previous master's solution
				#calculate cycle time as it may be different to master's solution
				self.optimalCycleTime = self.mostRecentFeasibleCycleTime
				self.gap[self.bendersIter] = self.gap[self.bendersIter-1]

				# calculate how many unique solutions of each sup-problem was found
				# self.uniqueSolutions = [None for k in self.inst.stations]
				# #only check uniqueness on previous Benders iters
				# for k in self.inst.stations:
				# 	stationSolutions = []
				# 	for i in range(self.bendersIter):
				# 		stationSolutions.append(tuple(self.all_solutions_ever[k][i]['tasks']))
				# 	self.uniqueSolutions[k] = set(stationSolutions)

				# calculate nodes statistics for master
				self.statsFinalMasterNodes = self.statsMasterNodes[-1]
				self.statsTotalMasterNodes = self.statsMasterNodes.sum()
				self.statsAvgMasterNodes = float(self.statsTotalMasterNodes/(self.bendersIter+1))

				if SUB_PROBLEM_SOLVER == 'mip' or 'cp' in SUB_PROBLEM_SOLVER:
					# calculate nodes statistics for sub problems
					self.statsFinalSubProblemNodes = 0
					self.statsTotalSubProblemNodes = sum([i.sum() for i in self.statsSubProblemNodes])
					self.statsAvgSubProblemNodes = float(self.statsTotalSubProblemNodes/((self.bendersIter)*self.inst.numStations))

			else:
				#calculate cycle time as it may be different to master's solution
				self.optimalCycleTime = self.mostRecentFeasibleCycleTime

				# calculate how many unique solutions of each sup-problem was found
				# self.uniqueSolutions = [None for k in self.inst.stations]
				# for k in self.inst.stations:
				# 	stationSolutions = []
				# 	for i in range(self.bendersIter+1):
				# 		stationSolutions.append(tuple(self.all_solutions_ever[k][i]['tasks']))
				# 	self.uniqueSolutions[k] = set(stationSolutions)

				# calculate nodes statistics for master
				self.statsFinalMasterNodes = self.statsMasterNodes[-1]
				self.statsTotalMasterNodes = self.statsMasterNodes.sum()
				self.statsAvgMasterNodes = float(self.statsTotalMasterNodes/(self.bendersIter+1))

				if SUB_PROBLEM_SOLVER == 'mip' or 'cp' in SUB_PROBLEM_SOLVER:
					# calculate nodes statistics for sub problems
					self.statsFinalSubProblemNodes = self.statsSubProblemNodes[-1].sum()
					self.statsTotalSubProblemNodes = sum([i.sum() for i in self.statsSubProblemNodes])
					self.statsAvgSubProblemNodes = float(self.statsTotalSubProblemNodes/((self.bendersIter+1)*self.inst.numStations))
		else:
			self.optimalCycleTime = 0

			# calculate nodes statistics for master
			self.statsFinalMasterNodes = self.statsMasterNodes[-1]
			self.statsTotalMasterNodes = self.statsMasterNodes.sum()
			self.statsAvgMasterNodes = self.statsTotalMasterNodes

			self.statsTotalSubProblemNodes = 0
			if SUB_PROBLEM_SOLVER == 'mip' or 'cp' in SUB_PROBLEM_SOLVER:
				# calculate nodes statistics for sub problems
				self.statsFinalSubProblemNodes = 0
				self.statsTotalSubProblemNodes = 0
				self.statsAvgSubProblemNodes = 0

		# calculate total number of cuts used
		cutAmountsList = [self.numGlobalLB, self.numNoGoods,
					 	  self.numGlobalUB, self.numInfAssCutsSimple,
						  self.numLogicCuts, self.numInfAssCutsSmart,
						  self.numInfAssCutsSmartest]
		self.numTotalCuts = sum(filter( lambda i: isinstance(i, int) , cutAmountsList))

		self.statsTotalNodes = self.statsTotalMasterNodes + self.statsTotalSubProblemNodes

		# calculate runtime stats
		self.statsTotalRuntime = round(sum(self.optimisation_times),4)
		self.statsMasterRuntime = round(sum(self.master_times),4)
		self.statsSubProbSolvetime = round(sum(self.sequencing_solve_times),4)
		self.statsSubProblemOverhead = round(sum(self.sequencing_overhead_times),4)
		self.statsSubProbRuntime = round(self.statsSubProbSolvetime + self.statsSubProblemOverhead,4)

	def print_solution(self):
		# output
		if args.human_readable:
			if not args.very_quiet:
				print('\n! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ')
				print('! \tSOLUTION ')
			else:
				print('\n',end='')
			print('! Cycle Time:\t{}'.format(round(self.optimalCycleTime)))
			if not args.very_quiet: 
				for k in self.inst.stations:
					print('! Station {}'.format(k))
					print('!   Load = \t{}'.format(round(self.curStationLoad[k])))
					print('!   Tasks = \t{}'.format(sorted(self.taskAssignment[k])))
					print('!   Starts = \t{}'.format(self.startTimes[k]))
		else:
			print(self.model.objval)

	def save_solution(self, results_file):
		# pdb.set_trace()
		with open(results_file, 'w', newline='') as csvfile:
			results = csv.writer(csvfile)
			results.writerow([self.solFeasible, self.solOptimal, self.optimalCycleTime, self.gap[self.bendersIter], 
							  self.statsTotalRuntime, self.init_time, self.statsMasterRuntime, self.statsSubProbRuntime, 
							  self.statsSubProbSolvetime, self.statsSubProblemOverhead,
							  self.bendersIter, self.statsTotalNodes, self.statsTotalMasterNodes, self.statsTotalSubProblemNodes,
							  self.numTotalCuts, self.numNoGoods, self.numGlobalLB, self.numGlobalUB,
							  self.numInfAssCutsSimple, self.numInfAssCutsSmart, self.numInfAssCutsSmartest, self.numLogicCuts])

	def print_statistics(self):
		if args.human_readable:
			print('\n! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ')
			print('! \tRUNTIME STATISTICS ')
			print('! Init time:\t{:.4f}'.format(self.init_time))
			print('! Total:\t{:.4f}'.format(self.statsTotalRuntime))
			print('! Maximum:\t{:.4f}'.format(max(self.optimisation_times)))
			print('! Average:\t{:.4f}'.format(sum(self.optimisation_times)/len(self.optimisation_times)))
			print('\n! Master times:')
			print('!   Total:\t{:.4f} ({:5.2f} %)'.format(self.statsMasterRuntime,
													   100*self.statsMasterRuntime/self.statsTotalRuntime))
			print('!   Maximum:\t{:.4f}'.format(max(self.master_times)))
			print('!   Average:\t{:.4f}'.format(self.statsMasterRuntime/len(self.master_times)))
			print('\n! Sequencing times:')
			if self.solFeasible:
				print('!   Total:\t{:.4f} ({:5.2f} %)'.format(self.statsSubProbRuntime,
														100*self.statsSubProbRuntime/self.statsTotalRuntime))
				print('!   Maximum:\t{:.4f}'.format(max(self.sequencing_solve_times + self.sequencing_overhead_times)))
				print('!   Average:\t{:.4f}'.format(self.statsSubProbSolvetime/len(self.sequencing_solve_times)))
				print('!   Solve time:\t{:.4f} ({:5.2f} %)'.format(self.statsSubProbSolvetime,
														   100*self.statsSubProbSolvetime/self.statsTotalRuntime))
				print('!   Overhead:\t{:.4f} ({:5.2f} %)'.format(self.statsSubProblemOverhead,
																100*self.statsSubProblemOverhead/self.statsTotalRuntime))
			else:
				print('!   No feasible solution found.')

			print('\n! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ')
			print('! \tSOLUTION STATISTICS ')
			print('! Feasible Solution:\t{}'.format(self.solFeasible))
			print('! Optimal Solution:\t{}'.format(self.solOptimal))
			print('! Master Iterations:\t{}'.format(self.bendersIter+1))
			print('! Gap:\t\t\t{:.2f}%'.format(self.gap[self.bendersIter]))
			# print('\n! Number of times sub-problem solved:')
			if self.solFeasible:
				# for k in self.inst.stations:
				# 	print('!   Station {}:\t\t{}'.format(k, len(self.uniqueSolutions[k])))
				print('\n! Total nodes explored:\t{}'.format(self.statsTotalNodes))
				print('! Master Nodes')
				print('!   TotalRMP:\t\t{}'.format(self.statsTotalMasterNodes))
				print('!   Final Master:\t{}'.format(self.statsFinalMasterNodes))
				print('!   AverageRMP:\t\t{:.2f}'.format(self.statsAvgMasterNodes))
				if SUB_PROBLEM_SOLVER == 'mip' or 'cp' in SUB_PROBLEM_SOLVER:
					print('! Sub-Problem')
					print('!   TotalSP:\t\t{}'.format(self.statsTotalSubProblemNodes))
					print('!   Final SPs Total:\t{}'.format(self.statsFinalSubProblemNodes))
					print('!   AverageSP:\t\t{:.2f}'.format(self.statsAvgSubProblemNodes))
			else:
				print('!   No feasible solution found.')
			print('\n! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ')
			print('! \tCUT STATISTICS ')
			print('! Total # Cuts:\t\t{}'.format(self.numTotalCuts))
			print('! Avg. Cuts per iter.:\t{:.2f}'.format(self.numTotalCuts/(self.bendersIter+1))) # don't include final iteration
			print('! Nogood Cuts:\t\t{}'.format(self.numNoGoods))
			print('! Global Lower Bounds:\t{}'.format(self.numGlobalLB))
			print('! Global Upper Bounds:\t{}'.format(self.numGlobalUB))
			print('! Simple Infer Cuts:\t{}'.format(self.numInfAssCutsSimple))
			print('! Smart Infer Cuts:\t{}'.format(self.numInfAssCutsSmart))
			print('! Smartest Infer Cuts:\t{}'.format(self.numInfAssCutsSmartest))
			print('! Logic Infeas. Cuts:\t{}'.format(self.numLogicCuts))
		else:
			print(self.init_time)
			print(sum(self.optimisation_times))
			print(max(self.optimisation_times))
			print(sum(self.optimisation_times)/len(self.optimisation_times))
			print(sum(self.master_times))
			print(max(self.master_times))
			print(sum(self.master_times)/len(self.master_times))
			print(sum(self.sequencing_solve_times))
			print(max(self.sequencing_solve_times))
			print(sum(self.sequencing_solve_times)/len(self.sequencing_solve_times))
			print(sum(self.statsSubProblemOverhead))
			print(self.bendersIter+1)
			for k in self.inst.stations:
				print(len(self.uniqueSolutions[k]))
			print(self.numTotalCuts)
			print(self.numTotalCuts/(self.bendersIter+1))
			print(self.numNoGoods)
			print(self.numGlobalLB)
			print(self.numGlobalUB)
			print(self.numInfAssCutsSimple)
			print(self.numInfAssCutsSmart)
			print(self.numInfAssCutsSmartest)
			print(self.numLogicCuts)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# BENDERS CALLBACK FUNCTIONALITY
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
def master_callback(model, where):
	# callback function to return to previously explored master B&B tree
	# adds the Benders cuts as lazy constraints

	# begin callback when a new incumbent MIP sol is found
	if where == GRB.Callback.MIPSOL:

		curMIPUB = round(s.model.cbGet(GRB.Callback.MIPSOL_OBJBST), 6)
		curMIPLB = round(s.model.cbGet(GRB.Callback.MIPSOL_OBJBND), 6)
		curMIPGap = 100*(curMIPUB - curMIPLB)/curMIPLB

		# move to sub problems once RMP is optimal (ie gap==0)
		if curMIPGap == 0.0:

			# store current RMP stats
			store_current_RMP_stats()

			if time.time() - s.startBenders > TIMELIMIT:
				s.master_timed_out = True
				s.doneBenders = True
				s.time_limit_exceeded = True

			if not args.very_quiet:
				print('\tCycle: \t{}'.format(round(s.curCycleTime)))
				print('\t\tUB: \t{}'.format(s.bestCycleTimeUB))
				print('\t\tGap: \t{:.2f} %\n'.format(s.gap[s.bendersIter]))

			# if gap found is 0 then we have already found a feasible solution to the sub-problems
			if s.gap[s.bendersIter] == 0:
				# ignore the current master solution and take the old one instead
				for k in s.inst.stations:
					s.taskAssignment[k] = s.all_solutions_ever[k][s.mostRecentUpperBoundIter]['tasks']
					s.all_solutions_ever[k].append({'tasks': s.taskAssignment[k]})
			else:
				# store current assignment
				if s.bendersIter > 0:
					s.statsSubProblemNodes.append(np.empty([0],dtype=int))
				for k in s.inst.stations:
					s.taskAssignment[k] = { i for i in s.inst.tasks if s.model.cbGetSolution(s.xs[i,k]) > 0.2 }
					s.all_solutions_ever[k].append({'tasks': s.taskAssignment[k]})

			# solve each sub-problem, adding cuts to master
			iterate_over_stations()

			if not False in s.stationFeasible:
				s.mostRecentFeasibleCycleTime = round(max(s.curStationLoad),4)

			# stopping condition: continuing until all sub-problem solutions <= master solution
			if False in s.stationSatisfiesCurCycleTime or False in s.stationFeasible:
				# add cuts and iterate if we arent done
				if USE_GLOBAL_BOUNDS:
					s.add_global_bounds(s.allowGlobalUB)
				s.bendersIter += 1
			else:
				# s.debug_final_result()
				# pdb.set_trace()
				s.doneBenders = True
			
def store_current_RMP_stats():
	s.master_times.append(time.time() - s.startMaster)
	s.optimisation_times.append(s.master_times[-1])

	s.curCycleTime = round(s.model.cbGetSolution(s.cycleTime),4)
	# s.statsMasterNodes = np.append(s.statsMasterNodes, int(s.model.nodecount))

	s.gap.append(round(float((s.bestCycleTimeUB - s.curCycleTime)/s.curCycleTime)*100,4))

def iterate_over_stations():
	for k in s.inst.stations:
		if not args.very_quiet:
			print(' Station %d' %(k), end='', flush=True)

		# define the time used up until starting this sub-problem
		s.SP_time_used = round(time.time() - s.startBenders,4)
		# check if we are out of time before starting each sub-problem
		if s.SP_time_used > TIMELIMIT:
			s.doneBenders = True
			s.time_limit_exceeded = True
			break

		# solve the current station's sub-problem
		result = s.solve_sub_problem(k)
		# check if time-limit is exceeded
		if s.time_limit_exceeded:
			doneBenders = True
			break

		# if we have already processed ths assignment before move onto next sub problem
		if result == True:
			s.allowGlobalUB = False


# Script to create instance class, run the solver and output the solution
if __name__ == '__main__':
	# start total runtime timer
	startAll = time.time()

	filename = args.file # retrieve filename of instance to solve
	if args.human_readable:
		print('Instance:', filename)
	else:
		print(filename)

	# store assembly line instance data
	if not args.very_quiet:
		print('Importing data... ', end='', flush=True)
	inst = AssemblyLineInstance(INST_DIR,filename)
	if not args.very_quiet:
		print('completed.')

	# Check for top-level infeasibility of the instance
	# if inst.obvious_infeasibility():
	# 	if args.human_readable:
	# 		print('infeasible!!!!')
	# 	else:
	# 		print('\n'.join('0'*9))
	# else:

	# create Solver for given instance and optimise it
	s = Solver(inst)

	s.benders_optimise_with_master_callbacks()

	# if SUB_PROBLEM_TYPE == 'opt':
	# 	s.benders_optimise_with_optimality_sub_problems()
	# elif SUB_PROBLEM_TYPE == 'feas':
	# 	s.benders_optimise_with_feasibility_sub_problems()

	# store output
	s.check_for_feasibility_and_optimality()
	s.process_solution_statistics()
	if s.solFeasible == 1:
		s.print_solution()

	# check solution consistency
	if CHECK_SOLUTION:
		print('\nVerifying solution:')
		check_solution_benders(s.optimalCycleTime, s.taskAssignment, s.startTimes, s.curStationLoad, inst)

	s.save_solution('summary_results_{}.txt'.format(EXPERIMENT_TOKEN))
	if PRINT_STATISTICS:
		s.print_statistics()

	# print total runtime
	end = time.time()
	if not args.quiet:
		print('\n! Language runtime:\t{:.4f}'.format(end-startAll-sum(s.optimisation_times)))
	if args.human_readable:
		if args.very_quiet:
			print('')
		print('\n! Total runtime:\t{:.4f}\n'.format(end-startAll))
	else:
		print(end-startAll)

	# pdb.set_trace()

# EOF #