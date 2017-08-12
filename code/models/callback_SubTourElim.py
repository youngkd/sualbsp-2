# Project: 	Masters Research
# Author: 	Kenneth Young
# Title:	Sub-tour elimination functions for the SUALBSP-2

# This file contains:
# 	-Function to eliminate sub-tours found in the relaxed master 
#	 solution of my Benders decomposition of the SUALBSP-2
# NOTE: UNFINISHED

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

# User-defined Functionality
from ALB_instance_storage import AssemblyLineInstance

# Define globals constants
if sys.platform == "win32":
	INST_DIR = 'instances\\'
elif sys.platform =="cygwin":
	INST_DIR = 'instances\\'
elif sys.platform == "darwin":
	INST_DIR = 'instances/'
elif sys.platform == "linux" or sys.platform == "linux2":
	INST_DIR = 'instances/'


# Callback Funtions - use lazy constraints to eliminate sub-tours
# Applicable to the relaxed master problem
def callback_sub_tour_elimination(model, where):
	if where == GRB.callback.MIPSOL: # perform sub-tour elimination when a new MIP solution is found
		# add STE constraints for each station
		for k in range(model._numStations):
			pdb.set_trace()
			tasks = { i for i in model._tasks if model._xs[i,k].x > 0.5 }
			n = len(tasks)
			selected = []
			# make a list of edges selected in the solution
			pdb.set_trace()
			for i in range(n-1):
				forwSol = model.cbGetSolution([model._ys[i,j] for j in tasks
															  if [i,j] in model._ys])
				selected += [(i,j) for j in range(n) if forwSol[j] > 0.5]
			backSol = model.cbGetSolution([model.ys[i,j] for j in range(n)])
			# find the shortest cycle in the selected edge list
			tour = find_sub_tour(n, selected)
			if len(tour) < n:
				# add a subtour elimination constraint
				expr = 0
				for i in range(len(tour)):
					for j in range(i+1, len(tour)):
						expr += model._vars[tour[i], tour[j]]
				model.cbLazy(expr <= len(tour)-1)

def find_sub_tour(n, edges):
	visited = [False]*n
	cycles = []
	lengths = []
	selected = [[] for i in range(n)]
	for x,y in edges:
		selected[x].append(y)
	while True:
		current = visited.index(False)
		thiscycle = [current]
		while True:
			visited[current] = True
			neighbors = [x for x in selected[current] if not visited[x]]
			if len(neighbors) == 0:
				break
			current = neighbors[0]
			thiscycle.append(current)
		cycles.append(thiscycle)
		lengths.append(len(thiscycle))
		if sum(lengths) == n:
			break
	return cycles[lengths.index(min(lengths))]

# EOF #