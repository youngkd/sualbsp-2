# Project: 	Masters Research
# Author: 	Kenneth Young
# Title:	Instance Storage for SUALBSP-2

# This file contains:
# 	-Classes and methods for storing intance data for the type-2 SUALBSP

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

# Define globals constants
if sys.platform == "win32":
	INST_DIR = 'instances\\'
elif sys.platform =="cygwin":
	INST_DIR = 'instances\\'
elif sys.platform == "darwin":
	INST_DIR = 'instances/'
elif sys.platform == "linux" or sys.platform == "linux2":
	INST_DIR = 'instances/'

# Class defining the instance of the overall assembly line
class AssemblyLineInstance:
	def __init__(self, instDir, instFilename):
		if instDir is None or instFilename is None:
			return None
		self.instDir = instDir
		self.instFilename = instFilename

		# import the instance
		self.import_instance_data()
		# self.create_backward_setups_array()

		# create the precedence graph for this instance
		self.create_precedence_graph()

		# create list of sets defining feasible assignments
		self.create_feasible_station_sets()
		self.create_feasible_task_sets()

		# find all predecessors and successors of each task
		self.allPredecessors = []
		for i in self.tasks:
			self.allPredecessors.append(self.find_all_predecessors(i,set()))
		self.allSuccessors = []
		for i in self.tasks:
			self.allSuccessors.append(self.find_all_successors(i,set()))

		# pdb.set_trace()
		# construct lists of sets defining allowed following and preceding of tasks
		self.construct_sets_of_allowed_followers_and_preceders_for_each_task()

		# calculate naive upper and lower bounds on the cycle time
		self.calculate_cycle_time_minimum_naive()
		self.calculate_cycle_time_maximum_naive()

		# define big-M value
		self.bigM = self.numStations # really? ...is that it?

	def import_instance_data_OLD(self):
		with open(self.instFilename) as f:
			# retrieve data
			instData = csv.reader(f)
			# instData = csv.reader(open(self.instDir+self.instFilename))

			self.numTasks = int(next(instData)[0])	# num of tasks
			self.numPrecs = int(next(instData)[0])	# num of precedence relations
			self.numStations = int(next(instData)[0])	# given num of stations

			# define index sets
			self.tasks = range(self.numTasks)
			self.precedences = range(self.numPrecs)
			self.stations = range(self.numStations)

			# read in processing times
			self.procList = []
			for task in self.tasks:
				self.procList.append(int(next(instData)[1]))

			# read in precedence relations
			self.precList = []
			for i in self.precedences:
				Prec = next(instData)
				self.precList.append((int(Prec[0]),int(Prec[1])))

			# read in forward setup times
			self.forwSU = []
			for taski in self.tasks:
				forwardSetupTimesList = next(instData)
				for taskj in range(len(forwardSetupTimesList)):
					forwardSetupTimesList[taskj] = int(forwardSetupTimesList[taskj])
				self.forwSU.append(forwardSetupTimesList)

	def import_instance_data(self):
		with open(self.instFilename) as f:
			# retrieve data
			instData = csv.reader(f)
			
			next(instData)
			self.numTasks = int(next(instData)[0])
			# define index set of tasks
			self.tasks = range(self.numTasks)

			# remove unneeded lines
			for row in range(2):
				next(instData)

			# read processing times
			self.procList = []
			for task in self.tasks:
				self.procList.append(int(next(instData)[0].split(' ')[1]))
			# remove unneeded lines
			for row in range(2):
				next(instData)

			# read precedence relations
			self.precList = []
			while True:
				prec = next(instData)
				# pdb.set_trace()
				if prec == []:
					break
				self.precList.append((int(prec[0]) - 1, int(prec[1]) - 1))
			# define index set of precedences
			self.numPrecs = len(self.precList)
			self.precedences = range(self.numPrecs)
			next(instData) # remove unneeded line

			# store precedence relations in an alternative manner
			# pdb.set_trace()
			self.altPrecList = [set([j for (idash,j) in self.precList if i==idash]) for i in self.tasks]

			# read in forward setup times
			forwardSetups = []
			while True:
				forwardSetupRel = next(instData)
				if forwardSetupRel == []:
					break
				forwardSetups.append(forwardSetupRel)

			# initialise all forward setup costs to infinity
			self.forwSU = np.full((self.numTasks, self.numTasks), 10000000, dtype='int64')
			np.fill_diagonal(self.forwSU, 0) # setup to itself is zero
			maxForwSU = 0
			for setupRelation in forwardSetups:
				preceder = int(setupRelation[0]) - 1
				tmp = setupRelation[1].split(':')
				follower = int(tmp[0]) - 1
				setupTime = int(tmp[1])
				if setupTime > maxForwSU:
					maxForwSU = setupTime
				self.forwSU[preceder,follower] = setupTime
			# change the infinity values to larger than the largest forward setup
			self.forwSU[self.forwSU > maxForwSU] = 1.2*maxForwSU

			next(instData) # remove unneeded line

			# read in backward setup times
			backwardSetups = []
			while True:
				backwardSetupRel = next(instData)
				if backwardSetupRel == []:
					break
				backwardSetups.append(backwardSetupRel)

			# initialise all backward setup costs to infinity
			self.backSU = np.full((self.numTasks, self.numTasks), 10000000, dtype='int64')
			np.fill_diagonal(self.backSU, 0) # setup to itself is zero
			maxBackSU = 0
			for setupRelation in backwardSetups:
				preceder = int(setupRelation[0]) - 1
				tmp = setupRelation[1].split(':')
				follower = int(tmp[0]) - 1
				setupTime = int(tmp[1])
				if setupTime > maxBackSU:
					maxBackSU = setupTime
				self.backSU[preceder,follower] = setupTime
			# change the infinity values to larger than the largest backward setup
			self.backSU[self.backSU > maxBackSU] = 1.2*maxBackSU

			# remove unneeded lines
			for i in range(3):
				next(instData)

			# read in number of stations (from optimal value in data file)
			self.numStations = int(next(instData)[0])
			# define index set of stations
			self.stations = range(self.numStations)

			# pdb.set_trace()

	def create_backward_setups_array(self):
		# create the backsetup times for the old type of input instances
		if BACKWARD_SETUP_TYPE == "copy-forward-setups":
			# Backward setup times equal forward setup times for now
			self.backSU = self.forwSU

	def create_precedence_graph(self):
		# initialise the graph
		self.precGraph = nx.DiGraph()
		self.precGraph.add_nodes_from(self.tasks)
		# add non-trivial edges
		for prec in self.precList:
			self.precGraph.add_edge(prec[0],prec[1])

	def calculate_cycle_time_minimum_naive(self):
		# assuming tasks can be perfectly divided across stations
		exactlyEvenCycleTime = np.ceil(sum(self.procList)/self.numStations)
		maxProcTime = max( self.procList[i] + self.backSU[i][i] for i in self.tasks )
		self.minCycleTime = max( maxProcTime, exactlyEvenCycleTime )

	def calculate_cycle_time_maximum_naive(self):
		# assuming all tasks in down in serial, ie. assigned the same station
		sumOfProcList = sum(self.procList)
		# assuming tasks are completed in input order
		naiveFeasibleOrderingCost = sum( [ self.forwSU[i][i+1] 
										   for i in self.tasks[:-1] ] ) + \
									self.backSU[self.numTasks-1][0]
		# together these give a naive feaible maximum
		self.maxCycleTime = sumOfProcList + naiveFeasibleOrderingCost

	def create_feasible_station_sets(self):
		# naively define each task to be possibly assigned any station
		self.feasibleStations = [ set(self.stations) for i in self.tasks ]

	def create_feasible_task_sets(self):
		# naively define each stations to be possibly assigned any task
		self.feasibleTasks = [ set(self.tasks) for i in self.stations ]

	def construct_sets_of_allowed_followers_and_preceders_for_each_task(self):
		self.followForw = [ (   set(self.tasks)
								 - (self.allSuccessors[i] - set(self.precGraph.successors(i)) )
								 - self.allPredecessors[i])
								 - set([i])
								for i in self.tasks ]
		self.precedeForw = [ set([ j for j in self.tasks 
									 if i in self.followForw[j] ])
								for i in self.tasks ]
		self.followBack = [ set(self.tasks) - self.allSuccessors[i]
								for i in self.tasks ]
		self.precedeBack = [ set([ j for j in self.tasks 
									 if i in self.followBack[j] ])
								for i in self.tasks ]

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

	def is_instance_obviousl_infeasible(self):
		return False


# EOF #