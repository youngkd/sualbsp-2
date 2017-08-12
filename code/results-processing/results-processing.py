# Masters Research Project 2017
# Kenneth Young
# Data Processing for the SUALBSP-2

# This file contains:
# 	-Functions to process the results

# Packages
from __future__ import division
import os, sys, pdb, csv
import numpy as np

#----------------------------------------------------------------------------------------#
# IMPORTING FILES

def import_results(ResultsDir, Model, SPsolver, classNum, alpha, creator, SearchAbrv):

	ResultsFilename = '{}results_{}_SP{}_class{}_alpha{}_{}_{}.txt'.format(ResultsDir,Model,SPsolver,classNum,
																		alpha,creator,SearchAbrv)
	try:
		results = csv.DictReader(open(ResultsFilename))
	except IOError:
		results = []

	return results

def store_all_instances_dict(ResultsDir,Model,SPsolver,classNum,alphaList,CreatorList,SearchAbvrList):

	ResultsArray = {alpha: {creator: {SearchAbrv: [] 
												for SearchAbrv in SearchAbvrList} 
												for creator in CreatorList} 
												for alpha in alphaList}

	for alpha in alphaList:
		for creator in CreatorList:
			for SearchAbrv in SearchAbvrList:
				results = import_results(ResultsDir, Model, SPsolver,classNum, alpha, creator, SearchAbrv)
				results = list(results)

				ResultsArray[alpha][creator][SearchAbrv].append(results)

	return ResultsArray


#----------------------------------------------------------------------------------------#
# SUPPORT FUNCTIONALITY
def define_lowercase_SP_solver(SPsolver, CPmodeltype):

	if SPsolver == 'MIP':
		lcSPsolver = 'mip'
	elif SPsolver == 'CP2' or SPsolver == 'CP3':
		if CPmodeltype == 2:
			lcSPsolver = 'cp2'
		elif CPmodeltype == 3:
			lcSPsolver = 'cp3'
	else:
		lcSPsolver = 'na'

	return lcSPsolver

def define_creator_list(classNum):

	if classNum == 1:
		CreatorList = ['mertens', 'bowman8', 'jaeschke', 'jackson', 'mansoor', 'mitchell']
	elif classNum == 2:
		CreatorList = ['roszieg', 'heskia', 'buxey', 'sawyer30']
	elif classNum == 3:
		CreatorList = ['lutz1', 'gunther', 'kilbrid', 'hahn', 'warnecke']
	elif classNum == 4:
		CreatorList = ['tonge', 'wee-mag', 'arc83']
	elif classNum == 5:
		CreatorList = ['lutz2', 'lutz3', 'mukherje']
	elif classNum == 6:
		CreatorList = ['arc111', 'barthold', 'barthol2']
	elif classNum == 7:
		CreatorList = ['scholl']

	return CreatorList

def write_tex_file_preamble(TableFile):

	TableFile.write("""
	\\documentclass[12pt]{article}
	{}
	\\usepackage{latexsym,amssymb,amsmath,epsfig,amsfonts,graphicx,url,fancyhdr,todonotes,verbatim}
	\\usepackage{float}
	\\usepackage{chngpage}
	\\usepackage{booktabs}

	\\setlength{\\topmargin}{-10pt}
	\\setlength{\headsep}{20pt}
	\\setlength{\headheight}{14.5pt}
	\\setlength{\\textheight}{600pt}
	\\setlength{\oddsidemargin}{0pt}
	\\setlength{\evensidemargin}{0pt}
	\\setlength{\\textwidth}{460pt}
	\\setlength{\parskip}{.30cm}
	\\parskip=10pt

	\\interdisplaylinepenalty=1000


	\\pagestyle{fancy}
	\\fancyhf{}
	\\rhead{Kenneth Young}
	\\lhead{Master Research Project - SUALBSP-2}
	\\cfoot{\\thepage}

	\\begin{document}

	""")

def write_tex_file_post(TableFile):

	TableFile.write("""

	\\end{document}""")

def write_tex_table_benders_preamble(TableFile, SPsolver, CPmodeltype, cutType):

	TableFile.write("""
	\\begin{table}[tpb]
		\\scriptsize
		\\begin{adjustwidth}{-.9in}{-.9in}
		""")

	if SPsolver != 'MIP':
		TableFile.write("\\caption{spsolver=%s-%s, cut=%s}" %(SPsolver, CPmodeltype, cutType))
	else:
		TableFile.write("\\caption{spsolver=%s, cut=%s}" %(SPsolver, cutType))

	TableFile.write("""
		\\centering
		\\vspace{2mm}
		\\begin{tabular}{ccl|rrrrrrrrrrr}
			\\toprule
	""")

def write_tex_table_benders_CP_preamble(TableFile, SPsolver, CPmodeltype, cutType):

	TableFile.write("""
	\\begin{table}[tpb]
		\\scriptsize
		\\begin{adjustwidth}{-.9in}{-.9in}
		""")

	if SPsolver != 'MIP':
		TableFile.write("\\caption{spsolver=%s-%s, cut=%s}" %(SPsolver, CPmodeltype, cutType))
	else:
		TableFile.write("\\caption{spsolver=%s, cut=%s}" %(SPsolver, cutType))

	TableFile.write("""
		\\centering
		\\vspace{2mm}
		\\begin{tabular}{cc|rrrrrrrrrrr}
			\\toprule
	""")

def write_tex_table_benders_post(TableFile):

	TableFile.write("""
			\\bottomrule
		\\end{tabular}
		\\end{adjustwidth}
	\\end{table}

	""")

def write_tex_table_benders_header_row_simple(TableFile):
	TableFile.write("""
			class & alpha & creator & \#nodes & iters & \#cuts & cuts & \\%gap & \#no sol & \#opt & \\%opt & RMP.rt & SP.rt & rt(s) \\\\\\midrule
	""")

def write_tex_table_benders_CP_header_row_simple(TableFile):
	TableFile.write("""
			search & class & \#CPnodes & iters & \#cuts & cuts & \\%gap & \#no sol & \#opt & \\%opt & RMP.rt & SP.rt & rt(s) \\\\\\midrule\n""")

def write_tex_table_mip_preamble(TableFile, mipModel):

	TableFile.write("""
	\\begin{table}[tpb]
		\\scriptsize
		\\begin{adjustwidth}{-.9in}{-.9in}
		""")
	
	TableFile.write("\\caption{MIP-model=%s}" %(mipModel))

	TableFile.write("""
		\\centering
		\\vspace{2mm}
		\\begin{tabular}{ccl|rrrrrr}
			\\toprule
	""")

def write_tex_table_mip_header_row_simple(TableFile):
	TableFile.write("""
			class & alpha & creator & \#nodes & \\%gap & \#no sol & \#opt & \\%opt & rt(s) \\\\\\midrule
	""")

def write_tex_table_mip_post(TableFile):

	TableFile.write("""
			\\bottomrule
		\\end{tabular}
		\\end{adjustwidth}
	\\end{table}

	""")

#----------------------------------------------------------------------------------------#
# SUMMARIZING RESULTS

# class storing values for the row output
class RowTotal():
	def __init__(self):
		self.TotalNodes = 0
		self.TotalSPNodes = 0
		self.TotalIters = 0
		self.TotalCuts = 0
		self.TCutList = [0,0,0,0,0,0]
		self.TotalPercGap = 0
		self.TotalFeas = 0
		self.TotalSubOptimal = 0
		self.TotalOpt = 0
		self.TotalRMPruntime = 0
		self.TotalSPruntime = 0
		self.TotalRuntime = 0

		# rSum specific attributes
		self.counter = 0
		self.TotalNoSol = 0
		self.TotalPercOpt = 0
		self.TotalNumInst = 0

# functions to write types of rows
def write_benders_row_simple_averages(TableFile,ResultsList,rSum,SPsolver):

	# initialise totals for output
	r = RowTotal()

	numInst = len(ResultsList)

	for instanceResult in ResultsList:
		# nodes
		r.TotalNodes += int(instanceResult['nodes'])

		# SP nodes
		r.TotalSPNodes += int(instanceResult['SPnodes'])

		# iters
		r.TotalIters += int(instanceResult['iters'])

		# cuts total
		r.TotalCuts += int(instanceResult['cuts'])

		# each cut total
		if instanceResult['numNG'] != '-':
			r.TCutList[0] += int(instanceResult['numNG'])
		if instanceResult['numGUB'] != '-':
			r.TCutList[1] += int(instanceResult['numGUB'])
		if instanceResult['numIC'] != '-':
			r.TCutList[2] += int(instanceResult['numIC'])
		if instanceResult['numIC2'] != '-':
			r.TCutList[3] += int(instanceResult['numIC2'])
		if instanceResult['numIC3'] != '-':
			r.TCutList[4] += int(instanceResult['numIC3'])
		if instanceResult['numLC'] != '-':
			r.TCutList[5] += int(instanceResult['numLC'])

		# feasibility
		feasible = int(instanceResult['feasible'])
		r.TotalFeas += feasible

		# optimality
		optimal = int(instanceResult['optimal'])
		r.TotalOpt += optimal

		# gap percentage
		# only add gap to running total if instance is sub-optimal
		if feasible == 1 and optimal != 1:
			r.TotalSubOptimal += 1
			r.TotalPercGap += float(instanceResult['gap'])

		# RMP runtime
		r.TotalRMPruntime += float(instanceResult['RMPtime'])

		# SP runtime
		r.TotalSPruntime += float(instanceResult['SPtime'])

		# total runtime
		r.TotalRuntime += float(instanceResult['runtime'])

	AvgNodes = r.TotalNodes/numInst
	AvgIters = r.TotalIters/numInst
	AvgCuts = r.TotalCuts/numInst
	AvgCutList = [ r.TCutList[i]/numInst for i in range(6) ]
	# AvgPercGap = r.TotalPercGap/numInst
	if r.TotalSubOptimal == 0:
		AvgPercGap = 0
	else:
		AvgPercGap = r.TotalPercGap/r.TotalSubOptimal
	r.TotalNoSol = numInst - r.TotalFeas
	r.TotalPercOpt = float(r.TotalOpt/numInst)*100
	AvgRMPruntime = r.TotalRMPruntime/numInst
	AvgSPruntime = r.TotalSPruntime/numInst
	AvgRuntime = r.TotalRuntime/numInst

	# if SPsolver == "MIP":
	TableFile.write('{:,.0f} & {:,.0f} & {:,.0f} & [{:,.0f},\\:{:,.0f},\\:{:,.0f},\\:{:,.0f},\\:{:,.0f},\\:{:,.0f}] & {:.2f} & {} & {}/{} & {:.2f} & {:.2f} & {:.2f} & {:.2f} \\\\\n'.format(AvgNodes,
						AvgIters,AvgCuts,
						AvgCutList[0],AvgCutList[1],AvgCutList[2],AvgCutList[3],AvgCutList[4],AvgCutList[5],
						AvgPercGap,r.TotalNoSol,r.TotalOpt,numInst,r.TotalPercOpt,
						AvgRMPruntime,AvgSPruntime,AvgRuntime))

	# keep track of running totals for this class
	rSum.counter += 1
	rSum.TotalNumInst +=	numInst
	rSum.TotalNodes += 		r.TotalNodes
	rSum.TotalSPNodes += 	r.TotalSPNodes
	rSum.TotalIters += 		r.TotalIters
	rSum.TotalCuts += 		r.TotalCuts
	rSum.TCutList[0] += 	r.TCutList[0]
	rSum.TCutList[1] += 	r.TCutList[1]
	rSum.TCutList[2] += 	r.TCutList[2]
	rSum.TCutList[3] += 	r.TCutList[3]
	rSum.TCutList[4] += 	r.TCutList[4]
	rSum.TCutList[5] += 	r.TCutList[5]
	rSum.TotalPercGap += 	r.TotalPercGap
	rSum.TotalNoSol += 		r.TotalNoSol
	rSum.TotalSubOptimal += r.TotalSubOptimal
	rSum.TotalOpt += 		r.TotalOpt
	rSum.TotalPercOpt += 	r.TotalPercOpt/100
	rSum.TotalRMPruntime += r.TotalRMPruntime
	rSum.TotalSPruntime +=	r.TotalSPruntime
	rSum.TotalRuntime += 	r.TotalRuntime

	return rSum

def write_benders_row_simple_summary(TableFile,rSum,SPsolver):

	rSum.TotalPercOpt = (rSum.TotalOpt/rSum.TotalNumInst)*100
	if rSum.TotalNumInst > rSum.TotalOpt:
		# totalSubOptimal = rSum.TotalNumInst - rSum.TotalOpt - rSum.TotalNoSol
		rSum.TotalPercGap = rSum.TotalPercGap/rSum.TotalSubOptimal
		# pdb.set_trace()
	else:
		rSum.TotalPercGap = 0

	rSum.TotalNodes = 		rSum.TotalNodes/rSum.TotalNumInst
	rSum.TotalSPNodes = 	rSum.TotalSPNodes/rSum.TotalNumInst
	rSum.TotalIters = 		rSum.TotalIters/rSum.TotalNumInst
	rSum.TotalCuts = 		rSum.TotalCuts/rSum.TotalNumInst
	rSum.TCutList[0] = 		rSum.TCutList[0]/rSum.TotalNumInst
	rSum.TCutList[1] = 		rSum.TCutList[1]/rSum.TotalNumInst
	rSum.TCutList[2] = 		rSum.TCutList[2]/rSum.TotalNumInst
	rSum.TCutList[3] = 		rSum.TCutList[3]/rSum.TotalNumInst
	rSum.TCutList[4] = 		rSum.TCutList[4]/rSum.TotalNumInst
	rSum.TCutList[5] = 		rSum.TCutList[5]/rSum.TotalNumInst
	# rSum.TotalPercGap = 	rSum.TotalPercGap/rSum.TotalNumInst
	# rSum.TotalPercOpt = 	rSum.TotalPercOpt/100
	rSum.TotalRMPruntime =	rSum.TotalRMPruntime/rSum.TotalNumInst
	rSum.TotalSPruntime =	rSum.TotalSPruntime/rSum.TotalNumInst
	rSum.TotalRuntime = 	rSum.TotalRuntime/rSum.TotalNumInst

	# if SPsolver == "MIP":
	# 	TableFile.write('\t\t\t\t\\midrule & & Overall\t& {:,.0f} & {:,.0f} & {:,.0f} & [{:,.0f},\\:{:,.0f},\\:{:,.0f},\\:{:,.0f},\\:{:,.0f},\\:{:,.0f}] & {:.2f} & {} & {}/{} & {:.2f} & {:.2f} & {:.2f} & {:.2f} \\\\\\midrule\n'.format(rSum.TotalNodes,
	# 		rSum.TotalIters,rSum.TotalCuts,
	# 		rSum.TCutList[0],rSum.TCutList[1],rSum.TCutList[2],rSum.TCutList[3],rSum.TCutList[4],rSum.TCutList[5],
	# 		rSum.TotalPercGap,rSum.TotalNoSol,rSum.TotalOpt,rSum.TotalNumInst,rSum.TotalPercOpt,
	# 		rSum.TotalRMPruntime,rSum.TotalSPruntime,rSum.TotalRuntime))
	# else:
	TableFile.write('\t\t\t\t\\midrule & & Overall\t& {:,.0f} & {:,.0f} & {:,.0f} & [{:,.0f},\\:{:,.0f},\\:{:,.0f},\\:{:,.0f},\\:{:,.0f},\\:{:,.0f}] & {:.2f} & {} & {}/{} & {:.2f} & {:.2f} & {:.2f} & {:.2f} \\\\\n'.format(rSum.TotalSPNodes,
		rSum.TotalIters,rSum.TotalCuts,
		rSum.TCutList[0],rSum.TCutList[1],rSum.TCutList[2],rSum.TCutList[3],rSum.TCutList[4],rSum.TCutList[5],
		rSum.TotalPercGap,rSum.TotalNoSol,rSum.TotalOpt,rSum.TotalNumInst,rSum.TotalPercOpt,
		rSum.TotalRMPruntime,rSum.TotalSPruntime,rSum.TotalRuntime))

def write_mip_row_simple_averages(TableFile,ResultsList,rSum):

	# initialise totals for output
	r = RowTotal()

	numInst = len(ResultsList)
	# pdb.set_trace()
	for instanceResult in ResultsList:
		# nodes
		r.TotalNodes += int(instanceResult['nodes'])

		# gap percentage
		if int(instanceResult['feasible']) == 1:
			r.TotalPercGap += float(instanceResult['gap'])

		# feasibility
		r.TotalFeas += int(instanceResult['feasible'])

		# optimality
		r.TotalOpt += int(instanceResult['optimal'])

		# total runtime
		r.TotalRuntime += float(instanceResult['runtime'])

	AvgNodes = r.TotalNodes/numInst
	AvgPercGap = r.TotalPercGap/(numInst-r.TotalNoSol)
	r.TotalNoSol = numInst - r.TotalFeas
	r.TotalPercOpt = float(r.TotalOpt/numInst)*100
	AvgRuntime = r.TotalRuntime/numInst

	TableFile.write('{:,.0f} & {:.2f} & {} & {}/{} & {:.2f} & {:.2f} \\\\\n'.format(AvgNodes,
				AvgPercGap,r.TotalNoSol,r.TotalOpt,numInst,r.TotalPercOpt,AvgRuntime))

	# keep track of running totals for this class
	rSum.counter += 1
	rSum.TotalNumInst +=	numInst
	rSum.TotalNodes += 		r.TotalNodes
	rSum.TotalPercGap += 	AvgPercGap
	rSum.TotalNoSol += 		r.TotalNoSol
	rSum.TotalOpt += 		r.TotalOpt
	rSum.TotalPercOpt += 	r.TotalPercOpt/100
	rSum.TotalRuntime += 	r.TotalRuntime

	return rSum

def write_mip_row_simple_summary(TableFile,rSum):

	rSum.TotalPercOpt = (rSum.TotalOpt/rSum.TotalNumInst)*100
	rSum.TotalPercGap = rSum.TotalPercGap/(rSum.counter)

	rSum.TotalNodes = 		rSum.TotalNodes/rSum.TotalNumInst
	rSum.TotalIters = 		rSum.TotalIters/rSum.TotalNumInst
	# rSum.TotalPercGap = 	rSum.TotalPercGap/rSum.TotalNumInst
	# rSum.TotalPercOpt = 	rSum.TotalPercOpt/100
	rSum.TotalRuntime = 	rSum.TotalRuntime/rSum.TotalNumInst

	TableFile.write('\t\t\t\t\\midrule & & Overall\t& {:,.0f} & {:.2f} & {} & {}/{} & {:.2f} & {:.2f} \\\\\\midrule\n'.format(rSum.TotalNodes,
			rSum.TotalPercGap,rSum.TotalNoSol,
			rSum.TotalOpt,rSum.TotalNumInst,rSum.TotalPercOpt,rSum.TotalRuntime))



#----------------------------------------------------------------------------------------#
# DATA PROCESSING DRIVERS

def main_basic_data_processing():

	debugging = 0

	# ~~~~~~~~~~~~~~~~~~~~~~
	# Input Parameters

	TIMELIMIT = 1800
	Model = "BD"
	# SPsolver = "MIP"
	SPsolver = "CP3"

	# ClassNumList = [1 ,2, 3]
	ClassNumList = [2]
	CPmodeltype = 3

	alphaList = ["1.00", "0.75", "0.50", "0.25"]
	# alphaList=["0.25"]

	if SPsolver == 'MIP':
		SearchList=["default_s"]
		SearchAbvrList = ['na']
	elif SPsolver == 'CP2' or SPsolver == 'CP3':
		# SearchList=["default_s", "priority_smallest"]
		# SearchAbvrList = ["def", "pris"]
		SearchList=["priority_first_fail"]
		SearchAbvrList = ["priff"]
	else:
		SearchList=["default_s"]
		SearchAbvrList = ['na']

	# Define Benders cutting procedure
	# CuttingList = ["001000", "000100", "000010", "011000", "010100", "010010", "010101", "011001"]
	CuttingList = ["010010"]

	# ~~~~~~~~~~~~~~~~~~~~~~
	# Output Directories
	TableDir = "tables-" + Model + "/"

	# initialise output instance file
	filename = TableDir + "table_{}_{}.tex".format(Model,SPsolver)
	TableFile = open(filename, 'w')

	# define the lowercase version of current solver for the data processing
	lcSPsolver = define_lowercase_SP_solver(SPsolver,CPmodeltype)

	write_tex_file_preamble(TableFile)

	for cutType in CuttingList:

		# Find directory containing desired redults
		if Model == "BD":
			ResultsDir = Model + "/" + SPsolver + "-" + cutType + "/results/"
		else:
			ResultsDir = Model + "/results/"

		# write start current tex table
		
		if Model == 'BD':
			write_tex_table_benders_preamble(TableFile, SPsolver, CPmodeltype, cutType)
			write_tex_table_benders_header_row_simple(TableFile)
		else:
			write_tex_table_mip_preamble(TableFile,Model)
			write_tex_table_mip_header_row_simple(TableFile)

		for classNum in ClassNumList:

			# find the right creators
			CreatorList = define_creator_list(classNum)

			# store all results for this class
			ResultsDict = store_all_instances_dict(ResultsDir,Model,lcSPsolver,classNum,alphaList,CreatorList,SearchAbvrList)

			# store the overall result for this class
			rSum = RowTotal()

			for alpha in alphaList:

				for creator in CreatorList:

					if alphaList.index(alpha) == 0:
						if CreatorList.index(creator) == 0:
							TableFile.write("\\midrule\t%s\t& %s\t& {\\tt %s}\t& " %(classNum, alpha, creator))
						else:
							TableFile.write("\t\t\t\t&\t\t& {\\tt %s}\t& " %(creator))
					else:
						if CreatorList.index(creator) == 0:
							TableFile.write("\t\t\t\t& %s\t& {\\tt %s}\t& " %(alpha, creator))
						else:
							TableFile.write("\t\t\t\t&\t\t& {\\tt %s}\t& " %(creator))
					# if (alpha=='0.50' and creator=='heskia') or (alpha=='0.25' and creator=='heskia'):
					# 	TableFile.write("\\\\\n")
					# 	continue
					for SearchAbrv in SearchAbvrList:

						ResultsList = ResultsDict[alpha][creator][SearchAbrv][0]

						try:
							# write content current tex table
							if Model == 'BD':
								rSum = write_benders_row_simple_averages(TableFile,ResultsList,rSum,SPsolver)
							else:
								rSum = write_mip_row_simple_averages(TableFile,ResultsList,rSum)
						except ZeroDivisionError:
							pdb.set_trace()


			# write overall averages for this class
			if Model == 'BD':
				write_benders_row_simple_summary(TableFile,rSum,SPsolver)
			else:
				write_mip_row_simple_summary(TableFile,rSum)

		# write end current tex table
		if Model == 'BD':
			write_tex_table_benders_post(TableFile)
		else:
			write_tex_table_mip_post(TableFile)
		
	# close the tex file
	write_tex_file_post(TableFile)
	TableFile.close

def main_basic_data_processing_CP():

	debugging = 0

	# ~~~~~~~~~~~~~~~~~~~~~~
	# Input Parameters

	TIMELIMIT = 1800
	Model = "BD"
	# SPsolver = "na"
	SPsolver = "CP"

	# ClassNumList = [1 ,2, 3]
	ClassNumList = [2]
	CPmodeltype = 3

	alphaList = ["1.00", "0.75", "0.50", "0.25"]
	# alphaList=["1.00"]

	if SPsolver == 'MIP':
		SearchList=["default_s"]
		SearchAbvrList = ['def']
	elif SPsolver == 'CP':
		# SearchList=["default_s", "start_s"]
		# SearchAbvrList = ["def", "s"]
		# SearchList=["start_s", "start_Then_spair", "priority_input_order", "priority_smallest", "priority_smallest_largest", "priority_first_fail"]
		# SearchAbvrList = ["s", "sTsp", "priio", "pris", "prisl", "priff"]
		SearchList=["priority_first_fail"]
		SearchAbvrList = ["priff"]
	else:
		SearchList=["default_s"]
		SearchAbvrList = ['na']

	# Define Benders cutting procedure
	# CuttingList = ["001000", "000100", "000010", "011000", "010100", "010010", "010101", "011001"]
	CuttingList = ["010010"]

	# ~~~~~~~~~~~~~~~~~~~~~~
	# Output Directories
	TableDir = "tables-" + Model + "/"

	# initialise output instance file
	# filename = TableDir + "table_{}_{}.tex".format(Model,SPsolver)
	# filename = TableDir + "table_{}_{}{}.tex".format(Model,SPsolver,CPmodeltype)
	filename = TableDir + "table_BD_CP3_round2.tex"
	TableFile = open(filename, 'w')

	# define the lowercase version of current solver for the data processing
	lcSPsolver = define_lowercase_SP_solver(SPsolver,CPmodeltype)

	write_tex_file_preamble(TableFile)

	for cutType in CuttingList:

		# Find directory containing desired redults
		if Model == "BD":
			if SPsolver == "MIP":
				ResultsDir = Model + "/" + SPsolver + "-" + cutType + "/results/"
			else:
				# ResultsDir = Model + "/CP-010010-cp23" + "/results/"
				ResultsDir = Model + "/CP3-010010" + "/results/"
		else:
			ResultsDir = Model + "/results/"

		# write start current tex table
		
		if Model == 'BD':
			if SPsolver == "MIP":
				write_tex_table_benders_preamble(TableFile, SPsolver, CPmodeltype, cutType)
				write_tex_table_benders_header_row_simple(TableFile)
			else:
				write_tex_table_benders_CP_preamble(TableFile, SPsolver, CPmodeltype, cutType)
				write_tex_table_benders_CP_header_row_simple(TableFile)
		else:
			write_tex_table_mip_preamble(TableFile,Model)
			write_tex_table_mip_header_row_simple(TableFile)


		for SearchAbrv in SearchAbvrList:
			TableFile.write("\t\t\t%s " %(SearchAbrv))

			for classNum in ClassNumList:
				if ClassNumList.index(classNum) == 0:
					TableFile.write("& %s & " %(classNum))
				else:
					TableFile.write("\t\t\t & %s & " %(classNum))

				# find the right creators
				CreatorList = define_creator_list(classNum)

				# store all results for this class
				ResultsDict = store_all_instances_dict(ResultsDir,Model,lcSPsolver,classNum,alphaList,CreatorList,SearchAbvrList)

				# store the overall result for this class
				rSum = RowTotal()

				for alpha in alphaList:

					for creator in CreatorList:

						# if alphaList.index(alpha) == 0:
						# 	if CreatorList.index(creator) == 0:
						# 		TableFile.write("\\midrule\t%s\t& %s\t& {\\tt %s}\t& " %(classNum, alpha, creator))
						# 	else:
						# 		TableFile.write("\t\t\t\t&\t\t& {\\tt %s}\t& " %(creator))
						# else:
						# 	if CreatorList.index(creator) == 0:
						# 		TableFile.write("\t\t\t\t& %s\t& {\\tt %s}\t& " %(alpha, creator))
						# 	else:
						# 		TableFile.write("\t\t\t\t&\t\t& {\\tt %s}\t& " %(creator))

						ResultsList = ResultsDict[alpha][creator][SearchAbrv][0]

						try:
							# write content current tex table
							if Model == 'BD':
								rSum = write_benders_row_simple_averages(TableFile,ResultsList,rSum,SPsolver)
							else:
								rSum = write_mip_row_simple_averages(TableFile,ResultsList,rSum)
						except ZeroDivisionError:
							pdb.set_trace()


				# write overall averages for this class
				if Model == 'BD':
					write_benders_row_simple_summary(TableFile,rSum,SPsolver)
				else:
					write_mip_row_simple_summary(TableFile,rSum)

		# write end current tex table
		if Model == 'BD':
			write_tex_table_benders_post(TableFile)
		else:
			write_tex_table_mip_post(TableFile)
		
	# close the tex file
	write_tex_file_post(TableFile)
	TableFile.close



if __name__ == "__main__":
	main_basic_data_processing()
	# main_basic_data_processing_CP()