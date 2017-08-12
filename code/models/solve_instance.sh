#!/bin/bash

SOURCE="$(pwd)"

if [[ -z "$1" ]] || [[ -z "$2" ]] || [[ -z "$3" ]] || [[ -z "$4" ]] || [[ -z "$5" ]] || [[ -z "$6" ]] || [[ -z "$7" ]] || [[ -z "$8" ]] || [[ -z "$9" ]] || [[ -z "$10" ]] || [[ -z "$11" ]] || [[ -z "$12" ]] || [[ -z "$13" ]] || [[ -z "$14" ]] || [[ -z "$15" ]] || [[ -z "$16" ]] || [[ -z "$17" ]] || [[ -z "$18" ]]; then
	echo 'Usage: ./solve_instance.sh <token> <model> <SPsolver> <SetNum> <ClassNum> <alpha> <creator> <NumInstances> <TIMELIMIT> <nc> <gb> <ic> <ic2> <ic3> <lc> <Search> <SPType> <debugging>'
	exit 1
else
	token="$1"
	model="$2"
	SPsolver="$3"
	SetNum="$4"
	ClassNum="$5"
	alpha="$6"
	creator="$7"
	NumInstances="$8"
	TIMELIMIT="$9"
	nc="${10}"
	gb="${11}"
	ic="${12}"
	ic2="${13}"
	ic3="${14}"
	lc="${15}"
	Search="${16}"
	SPType="${17}"
	debugging="${18}"
fi

# define directory of input instance
DataDir="Type-2-SBF${SetNum}/Type-2-alpha${alpha}/"

# define the model abbreviation
if [[ $model == "sualbsp2_fsbf" ]]; then
	modelAbrv="FSBF"
elif [[ $model == "sualbsp2_ssbf" ]]; then
	modelAbrv="SSBF"
elif [[ $model == "sualbsp2_scbf" ]]; then
	modelAbrv="SCBF"
elif [[ $model == "sualbsp2_benders" ]]; then
	modelAbrv="BD"
fi

# define the flags to call the python script with
flags="-H -s -q -t ${TIMELIMIT} -et ${token} -spt ${SPType}"
if [[ $model == "sualbsp2_benders" ]]; then
	if [[ $nc == "1" ]]; then
		flags="${flags} -nc"
	fi
	if [[ $gb == "1" ]]; then
		flags="${flags} -gb"
	fi
	if [[ $ic == "1" ]]; then
		flags="${flags} -ic"
	fi
	if [[ $ic2 == "1" ]]; then
		flags="${flags} -ic2"
	fi
	if [[ $ic3 == "1" ]]; then
		flags="${flags} -ic3"
	fi
	if [[ $lc == "1" ]]; then
		flags="${flags} -lc"
	fi
	if [[ $SPsolver == "mip" ]]; then
		flags="${flags} -sps mip"
		Search="na"
		searchAbrv="na"
	elif [[ $SPsolver == "cp2" ]]; then
		flags="${flags} -sps cp2 -cps ${Search}"
	elif [[ $SPsolver == "cp3" ]]; then
		flags="${flags} -sps cp3 -cps ${Search}"
	fi
fi

if [[ $SPsolver != "na" ]]; then
	if [[ $SPsolver != "mip" ]]; then
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
	fi
else
	Search="na"
	searchAbrv="na"
fi

# print current parameter details
printf '%s\n' " model |  SP | class | alpha |  creator |  # | search | flags "
printf '%s\n' "-------|-----|-------|-------|----------|----|--------|------------------"
printf " %5s | %3s | %5d | %5s | %8s | %2d | %6s | %5s \n" "${modelAbrv}" "${SPsolver}" ${ClassNum} "${alpha}" "${creator}" ${NumInstances} "${searchAbrv}" "${flags}"

# define destinations for the output files
RESULTS="./results/results_${modelAbrv}_SP${SPsolver}_class${ClassNum}_alpha${alpha}_${creator}_$searchAbrv.txt"

# clear the files if they already exist
printf '' > "$RESULTS"

# Get all instances
INSTANCES="$(find $DataDir -name "${creator}*.alb" | sort)"
instNum=0	# initilise instance counter

for inst in $INSTANCES
do
	# create and clear full-output file
	instNumLeadingZeros=$(printf "%02d" ${instNum})
	FULLOUTPUT="./full-output/full-output_${modelAbrv}_SP${SPsolver}_class${ClassNum}_alpha${alpha}_${creator}_${searchAbrv}_${instNumLeadingZeros}.txt"
	printf '' > "$FULLOUTPUT"

	# call the solver with the specified flags on the curent instance and store the full output
	python3 "${model}.py" $flags $inst > "$FULLOUTPUT"

	# create header for results file if on first instances
	if [ $instNum -eq 0 ]; then
		if [[ $model == "sualbsp2_benders" ]]; then
			if [ $instNum -eq 0 ]; then
				printf 'feasible,optimal,cycle,gap,runtime,inittime,RMPtime,SPtime,SPsolvetime,SPoverheadtime,iters,nodes,RMPnodes,SPnodes,cuts,numNG,numGLB,numGUB,numIC,numIC2,numIC3,numLC\n' >> "$RESULTS"
			else
				true
			fi
		else
			printf 'feasible,optimal,cycle,gap,runtime,nodes\n' >> "$RESULTS"
		fi
	fi

	# append results to results file
	cat "summary_results_${token}.txt" >> "$RESULTS"

	# find the status of the current instance (actually that's pointless... commented out)
	# if [[ $OSTYPE == *"darwin"* ]]; then
	status="$(cat summary_results_${token}.txt | awk -F ',' '{print $2}')"
	# elif [[ $OSTYPE == *"linux"* ]]; then
	# 	status="$(cat summary_results.txt | grep -oP '^[^0-9]*\K[-0-9]+')"
	# fi

	# print a 1 if found optimal, otherwise 0
	if [ $debugging -eq 1 ]; then
		echo ""
	else
		if [ "${status}" == "1" ]; then
			printf "1"
		else
			printf "0"
		fi
	fi
	instNum=$((instNum+1))

done

# EOF #