# Project: 	Masters Research
# Author: 	Kenneth Young
# Title:	Solution Verifier for SUALBSP-2

# This file contains:
# 	-Functions to confirm the correctness of a solution
#	 to the type-2 SUALBSP

# Packages
import sys
import pdb 
import time
# import itertools
import csv
import argparse
import numpy as np
import networkx as nx


# Define globals constants
if sys.platform == "win32":
	INST_DIR = 'instances\\'
elif sys.platform =="cygwin":
	INST_DIR = 'instances\\'
elif sys.platform == "darwin":
	INST_DIR = 'instances/'
elif sys.platform == "linux" or sys.platform == "linux2":
	INST_DIR = 'instances/'


def check_solution_benders(cycleTime, assignment, startTimes, load, inst):
	correct=True
	print('\n!Cycle Time: \t{}\n'.format(cycleTime))
	for k in inst.stations:
		print('! Station {} Load: \t{}'.format(k,round(load[k])))
		for index,i in enumerate(list(assignment[k])[:-1]):
			follower=index+1
			# pdb.set_trace()
			endTime=startTimes[k][index]+inst.procList[i]+inst.forwSU[i][list(assignment[k])[follower]]
			print('! -Task {}: [{},{}]'.format(i,startTimes[k][index],endTime))

			if endTime > cycleTime:
				correct=False

		i=list(assignment[k])[-1]
		endTime=startTimes[k][-1]+inst.procList[i]+inst.backSU[i][list(assignment[k])[0]]
		print('! -Task {}: [{},{}]\n'.format(i,startTimes[k][-1],endTime))

		if endTime > cycleTime:
			correct=False

	# pdb.set_trace()

	if correct:
		print('! ~~~ Solution is correct ~~~\n')
	else:
		print('! ~~~ ERROR: Solution is incorrect! ~~~\n')

def check_solution_SSBF(cycleTime, assignment, gs, hs, inst):
	correct=True

	startTimes = [None for i in inst.tasks]

	print('\n!Cycle Time: \t{}\n'.format(cycleTime))
	for k in inst.stations:
		print('! Station {}:'.format(k))
		for i in assignment[k]:
			for j in assignment[k]:
				try:
					if hs[i,j,k].x > 0.5:
						first = j
						last = i
				except KeyError:
					continue
		# pdb.set_trace()
		
		startTimes[first] = 0
		curTask = first

		# find the start times
		if first != last:
			while True:
				for j in inst.followForw[curTask]:
					if gs[curTask,j,k].x > 0.5:
						nextTask = j
				endTime = startTimes[curTask] + inst.procList[curTask] + inst.forwSU[curTask][nextTask]
				print('! -Task {}: [{},{}]'.format(curTask,startTimes[curTask],endTime))
				curTask = nextTask
				startTimes[curTask] = endTime

				if endTime > cycleTime:
					correct=False

				# exit wen we find the last task
				if curTask == last:
					endTime = startTimes[curTask] + inst.procList[curTask] + inst.backSU[curTask][first]
					print('! -Task {}: [{},{}]\n'.format(curTask,startTimes[curTask],endTime))

					if endTime > cycleTime:
						correct=False

					break
		else:
			endTime = startTimes[curTask] + inst.procList[curTask]
			print('! -Task {}: [{},{}]\n'.format(curTask,startTimes[curTask],endTime))

			if endTime > cycleTime:
				correct=False

		# for index,i in enumerate(list(assignment[k])[:-1]):
		# 	follower=index+1
		# 	pdb.set_trace()
		# 	endTime=startTimes[k][index]+inst.procList[i]+inst.forwSU[i][list(assignment[k])[follower]]
		# 	print(' -Task {}: [{},{}]'.format(i,startTimes[k][index],endTime))

		# i=list(assignment[k])[-1]
		# endTime=startTimes[k][-1]+inst.procList[i]+inst.backSU[i][list(assignment[k])[0]]
		# print(' -Task {}: [{},{}]\n'.format(i,startTimes[k][-1],endTime))

	# pdb.set_trace()

	if correct:
		print('! ~~~ Solution is correct ~~~\n')
	else:
		print('! ~~~ ERROR: Solution is incorrect! ~~~\n')



# EOF #