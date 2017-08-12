#!/bin/bash

# if [[ -z "$1" ]]; then
# 	token=0
# else
# 	token="$1"
# fi
token=0

debugging=0

# ~~~~~~~~~~~~~~~~~~
# Input parameters

TIMELIMIT="1800"
# ModelList=("sualbsp2_fsbf" "sualbsp2_ssbf" "sualbsp2_scbf" "sualbsp2_benders")
SetNumList=("2")
# SPsolverList=("mip" "cp2" "cp3") # "tsp")
# ClassNumList=("1" "2" "3")
# TIMELIMIT="10"
ModelList=("sualbsp2_benders")
ClassNumList=("2")
SPType="opt" # use optimality (opt) or feasibility (feas) sub-problems for Benders
# SPsolverList=("mip") # only used if doing Benders
SPsolverList=("cp3") # only used if doing Benders
# CP sub-problems 02
# SearchList=("default_s" "start_s")
# CP sub-problems 03
# SearchList=("defaut_s" "start_s" "startpair_s" "start_Then_startpair" "startpair_Then_start" "io_Then_startpair" "s_Then_startpair" "sl_Then_startpair" "ff_Then_startpair")
# SearchList=("s_Then_startpair")
# SearchList=("defaut_s" "start_s" "startpair_s" "start_Then_startpair" "startpair_Then_start" "priority_input_order" "priority_smallest" "priority_smallest_largest" "priority_first_fail")
SearchList=("priority_first_fail")
# SearchList=("priority_input_order" "priority_smallest" "priority_smallest_largest" "priority_first_fail")


# if [[ $token -eq 1 ]]; then
# 	alphaList=("1.00")
# elif [[ $token -eq 2 ]]; then
# 	alphaList=("0.75")
# elif [[ $token -eq 3 ]]; then
# 	alphaList=("0.50")
# elif [[ $token -eq 4 ]]; then
# 	alphaList=("0.25")
# fi

# alphaList=("1.00" "0.75" "0.50" "0.25")
alphaList=("0.50" "0.25")

# Benders cutting properties
nc="0"
gb="1"
ic="0"
ic2="0"
ic3="1"
lc="0"

# ~~~~~~~~~~~~~~~~~~
# Iterate over input

# iterate of each chosen model
for model in "${ModelList[@]}"
do
	# define the model abbreviation
	# also if not using Benders then don't iterateover the list of sub-problem solvers
	if [[ $model == "sualbsp2_fsbf" ]]; then
		modelAbrv="FSBF"
		SPsolverList=("na")
	elif [[ $model == "sualbsp2_ssbf" ]]; then
		modelAbrv="SSBF"
		SPsolverList=("na")
	elif [[ $model == "sualbsp2_scbf" ]]; then
		modelAbrv="SCBF"
		SPsolverList=("na")
	elif [[ $model == "sualbsp2_benders" ]]; then
		modelAbrv="BD"
	fi

	# iterate of each specified sub-problem solver
	for SPsolver in "${SPsolverList[@]}"
	do
		if [[ $SPsolver == "na" ]]; then
			SearchList=("na")
		elif [[ $SPsolver == "mip" ]]; then
			SearchList=("na")
			searchAbrv="na"
		fi
		# iterate over the CP search strategies
		for Search in "${SearchList[@]}"
		do

			if [[ $SPsolver != "na" ]]; then
				if [[ $Search == "default_s" ]]; then
					searchAbrv="def"
				elif [[ $Search == "start_s" ]]; then
					searchAbrv="s"
				elif [[ $Search == "startpair_s" ]]; then
					searchAbrv="sp"
				elif [[ $Search == "start_Then_startpair" ]]; then
					searchAbrv="sTsp"
				elif [[ $Search == "startpair_Then_start" ]]; then
					searchAbrv="spTs"
				elif [[ $Search == "priority_input_order" ]]; then
					searchAbrv="priio"
				elif [[ $Search == "priority_smallest" ]]; then
					searchAbrv="pris"
				elif [[ $Search == "priority_smallest_largest" ]]; then
					searchAbrv="prisl"
				elif [[ $Search == "priority_first_fail" ]]; then
					searchAbrv="priff"
				elif [[ $Search == "io_Then_startpair" ]]; then
					searchAbrv="ioSeq"
				elif [[ $Search == "s_Then_startpair" ]]; then
					searchAbrv="sSeq"
				elif [[ $Search == "sl_Then_startpair" ]]; then
					searchAbrv="slSeq"
				elif [[ $Search == "ff_Then_startpair" ]]; then
					searchAbrv="ffSeq"
				fi
			else
				searchAbrv="na"
			fi

			# iterate over all specified datasets
			for SetNum in "${SetNumList[@]}"
			do
				# iterate over classes
				for ClassNum in "${ClassNumList[@]}"
				do
					# define list of creator's names for each class of data
					if [[ $ClassNum -eq 1 ]]; then
						# n in [0,21]
						creatorNameList=("bowman8" "jackson" "jaeschke" "mansoor" "mertens" "mitchell")
						# creatorNameList=("jackson" "jaeschke")
					elif [[ $ClassNum -eq 2 ]]; then
						# n in [25,30]
						# creatorNameList=("buxey" "roszieg" "sawyer30" "heskia")
						creatorNameList=("heskia")
					elif [[ $ClassNum -eq 3 ]]; then
						# n in [32,58]
						creatorNameList=("lutz1" "gunther" "kilbrid" "hahn" "warnecke")
						# creatorNameList=("hahn")
					elif [[ $ClassNum -eq 4 ]]; then
						# n in [70,83]
						creatorNameList=("tonge" "wee-mag" "arc83")
					elif [[ $ClassNum -eq 5 ]]; then
						# n in [89,94]
						creatorNameList=("lutz2" "lutz3" "mukherje")
					elif [[ $ClassNum -eq 6 ]]; then
						# n in [111,148]
						creatorNameList=("arc111" "barthold" "barthol2")
					elif [[ $ClassNum -eq 7 ]]; then
						# n = 297
						creatorNameList=("scholl")
					fi

					# iterate over alpha values
					for alpha in "${alphaList[@]}"
					do
						# iterate over creator's names
						for creator in "${creatorNameList[@]}"
						do
							# find the number of instances with these parameter values
							DataDir="Type-2-SBF${SetNum}/Type-2-alpha${alpha}/"
							NumInstances=$(find $DataDir -name "$creator*.alb" | wc -l)

							if [ $NumInstances -ge 1 ]; then
								# solve all instances with the given parameters
								./solve_instance.sh $token $model $SPsolver $SetNum $ClassNum $alpha $creator $NumInstances $TIMELIMIT $nc $gb $ic $ic2 $ic3 $lc $Search $SPType $debugging

								# print the number of instances where optimality was found
								results="./results/results_${modelAbrv}_SP${SPsolver}_class${ClassNum}_alpha${alpha}_${creator}_${searchAbrv}.txt"
								awk -F "," '{x+=$2}END{printf "-> " x}' "$results"
								printf "/%s optimal\n\n" ${NumInstances}
							fi
							
						done
					done
				done
			done
		done
	done
done

# EOF #