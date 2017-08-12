#!/bin/bash

# SOURCE="$(pwd)"

# data files to process
SBFSetNum=1
alphaList=("1.00" "0.75" "0.50" "0.25")

# iterate over all specified datasets
for alpha in "${alphaList[@]}"
do
	# define output directory
	if [ $SBFSetNum -eq 1 ]; then
		outputDir="Type-2-data-sets/Type-2-SBF1/Type-2-alpha${alpha}"
		dataDir="./SBF-data-sets/DataSet-SBF1/SBF1-${alpha}"
	elif [ $SBFSetNum -eq 2 ]; then
		outputDir="Type-2-data-sets/Type-2-SBF2/Type-2-alpha${alpha}"
		dataDir="./SBF-data-sets/DataSet-SBF2/SBF2-${alpha}"
	fi

	INSTANCES="$(find ${dataDir} -name "*.alb")"
	for filename in $INSTANCES
	do
		#store new data file's name variables
		creator=$(echo $filename | cut -d "/" -f 5 | cut -d "_" -f 1 | tr '[:upper:]' '[:lower:]')
		n=$(cat ${filename} | grep -A1 'number of tasks' | grep -v 'number of tasks')
		k=$(cat ${filename} | grep -A1 'optimal' | grep -v 'optimal')
		instNum=$( printf %02d $(find ${outputDir} -iname "${creator}_n${n}_k${k}_*.alb" | wc -l) )
		#create new file's name
		new_file="${outputDir}/${creator}_n${n}_k${k}_${instNum}.alb"
		#create file
		echo "$(cat ${filename} | sed -e '/cycle/ { N; N; d; }')" > "$new_file"
		# echo ''
	done
done

# EOF #