#!/bin/bash

ncores=1
conf="/tracer/docker_helper_files/docker_tracer.conf"
resources="../resources"
readID="_R"
organism="Mmus"
single="false"
fragment_length=75
fragment_sd=2
function usage () {
    cat >&2 <<EOF
USAGE: SS2_tracer [options] <fastq file directory>
1. Merge multiple Lanes from the sequencing .fastq.gz output into one .fastq file. 
2. Then perform tracer assembly on SmartSeq2 data.
3. The results from tracer are saved in outdir/final, for each cell a folder is created.

-o <path/to/outdir>   : Output directory, where the merge .fastq files under outdir/merged and tracer results under outdir/final
are be saved.
-c number of cores    : Number of cores, which TRACER should use. Default: 1
-f configuration file : Configuration file for tracer. Default: container specific configuration file
-r resources for TCR reconstruction  : Path to folder for resources. Default "../resources", considering that you execute the script from scripts folder of the tracer repository
-i identifier in the fastq name, distinguishing Read 1 and Read 2. Default "_R"
-a organism: Which organism are you using. Default: Mmus
-s single end?: use this flag, if single end sequencing strategy is used. Default:false
-l if single end is used: determine fragment length. Default 75
-d if single end is used: determine fragment standard deviation: Default 2


EOF
}

function check_set() {
    value=$1
    name=$2
    flag=$3

    if [[ -z "$value" ]]
    then error_exit "$name has not been specified.  $flag flag is required"
    fi
}
function error_exit() {
    echo "ERROR: $1
    " >&2
    usage
    exit 1
}

 while getopts ":o:c:f:r:i:s:l:d:" options; do
   case $options in
     o ) output=$OPTARG;;
     c ) ncores=$OPTARG;;
     f ) conf=$OPTARG;;
     r ) resources=$OPTARG;;
     i ) readID=$OPTARG;;
     s ) single=$OPTARG;;
     l ) fragment_length=$OPTARG;;
     d ) fragment_sd=$OPTARG;;
     a ) organism=$OPTARG;;


   esac
 done
 shift $((OPTIND - 1))
 
 check_set "$output" "output" "-o"
 check_set "$ncores" "ncores" "-c"
fastq_files=$1

MAX_POOL_SIZE=${ncores}
CURRENT_POOL_SIZE=0

if [ ! -d ${output} ]
then
mkdir ${output}
fi;

if [ ! -d "${output}/merged" ] ##first merge the fastq files and unzip to a separate folder
then
mkdir "${output}/merged"
fi;

for x in `find "$fastq_files" -maxdepth 1 -name "*${readID}1*.fastq*" -type f`; do
file="$(basename -- $x)"
##splits at read identifier in name to get sample name
ID=${file%%"$readID"*} 
curdir=${fastq_files}
out="${output}/merged"
echo ${ID}
endR1=${file##*"$readID"}   ### this is the structure '_R1_001.fastq.gz'


endR1q=${endR1%%.gz}   #### output is not gunzipped '_R1_001.fastq'

file_outq_1="${out}/${ID}${readID}${endR1q}"

file_out_1="${out}/${ID}${readID}${endR1}"

###merge R1

count=$(find "$curdir" -maxdepth 1 -type f -name "${ID}*${endR1}" | wc -l)


if [[ ! -e ${file_out_1} && ! -e ${file_outq_1} && count -eq 1 ]]
then
echo "merge" ${ID}
find "${curdir}" -name "${ID}*${endR1}" -type f -exec cat {} > ${file_out_1} \; &
else
  echo "either files do not exist or only 1 fastq exists, skipping merge"
fi; 



# When a new job is created, the program updates the $CURRENT_POOL_SIZE variable before next iteration
	CURRENT_POOL_SIZE=$(jobs | wc -l)
	echo $CURRENT_POOL_SIZE
###only 50 processes simultaneously
	while [[ $CURRENT_POOL_SIZE -ge $MAX_POOL_SIZE ]]; do    
    # The above "echo" and "sleep" is for demo purposes only.
    # In a real usecase, remove those two and keep only the following line
    # It will drastically increase the performance of the script
		CURRENT_POOL_SIZE=$(jobs | wc -l)

done 
done


if [[ $single == "false" ]]
then
#wait 
#
####merge R2
for x in `find "$fastq_files" -maxdepth 1 -name "*${readID}2*.fastq*" -type f`; do
file="$(basename -- $x)"
##splits at read identifier in name to get sample name
ID=${file%%"$readID"*} 
curdir=${fastq_files}
out="${output}/merged"
echo ${ID}
endR2=${file##*"$readID"}   ### this is the structure '_R2_001.fastq.gz'


endR2q=${endR2%%.gz}   #### output is not gunzipped '_R2_001.fastq'

file_outq_2="${out}/${ID}${readID}${endR2q}"

file_out_2="${out}/${ID}${readID}${endR2}"

###merge R2

count=$(find "$curdir" -maxdepth 1 -type f -name "${ID}*${endR2}" | wc -l)


if [[ ! -e ${file_out_2} && ! -e ${file_outq_2} && count -eq 1 ]]
then
echo "merge" ${ID}
find "${curdir}" -name "${ID}*${endR2}" -type f -exec cat {} > ${file_out_2} \; &
else
  echo "either files do not exist or only 1 fastq exists, skipping merge"
fi; 


# When a new job is created, the program updates the $CURRENT_POOL_SIZE variable before next iteration
	CURRENT_POOL_SIZE=$(jobs | wc -l)
	echo $CURRENT_POOL_SIZE
###only 50 processes simultaneously
	while [[ $CURRENT_POOL_SIZE -ge $MAX_POOL_SIZE ]]; do    
    # The above "echo" and "sleep" is for demo purposes only.
    # In a real usecase, remove those two and keep only the following line
    # It will drastically increase the performance of the script
		CURRENT_POOL_SIZE=$(jobs | wc -l)

done 
done
wait 
fi;

echo "fastqs merged"
#cd ${out}
count=$(find "$curdir" -maxdepth 1 -type f -name "*.fastq.gz" | wc -l)
echo "$count"
if [[ ! count -eq 0 ]]
then
MAX_POOL_SIZE=${ncores}
CURRENT_POOL_SIZE=0

for x in `find "${out}" -name "*.fastq.gz" -type f`; do


echo "gunzip" ${x}
gunzip "$x" &


# When a new job is created, the program updates the $CURRENT_POOL_SIZE variable before next iteration
	CURRENT_POOL_SIZE=$(jobs | wc -l)
	echo $CURRENT_POOL_SIZE
###only 50 processes simultaneously
	while [[ $CURRENT_POOL_SIZE -ge $MAX_POOL_SIZE ]]; do    
    # The above "echo" and "sleep" is for demo purposes only.
    # In a real usecase, remove those two and keep only the following line
    # It will drastically increase the performance of the script
		CURRENT_POOL_SIZE=$(jobs | wc -l)

done 
done
else
echo "fastq already unzipped, skipping"

wait 
fi;

MAX_POOL_SIZE=${ncores}
CURRENT_POOL_SIZE=0


if [ ! -d "${output}/final" ]
then
mkdir "${output}/final"
fi;

for x in `find "$out" -maxdepth 1 -name "*${readID}1*.fastq*" -type f`; do
file="$(basename -- $x)"

ID=${file%%"$readID"*} 
curdir=${fastq_files}
out="${output}/merged"
echo ${ID}
#endR1q=${file##*"${readID}1"*}   ### this is the structure '_R2_001.fastq.gz'
endR=${file##*"$readID"1}   ### this is the structure '_R2_001.fastq.gz'

if [[ $single == "false" ]]
then
file_outq_2="${out}/${ID}${readID}2${endR}"
fi;
file_outq_1="${out}/${ID}${readID}1${endR}"


echo "TRACER" ${PWD}

if [ ! -d "${out}${ID}" ] ##don't redo analysis, if already done
then
echo "TRACER" ${file}
if [[  $single == "false" ]] ###if paired end
then
singularity run --bind $PWD:/scratch -W /scratch\
 docker://teichlab/tracer assemble -c ${conf} --resource_dir ${resources} --ncores ${ncores} -q kallisto -r -s $organism \
${file_outq_1} ${file_outq_2} ${ID} "${output}/final" &
else ##single end
singularity run --bind $PWD:/scratch -W /scratch\
 docker://teichlab/tracer assemble -c ${conf} --resource_dir ${resources} --ncores ${ncores} -q kallisto -r -s $organism --single_end --fragment_length ${fragment_length} \
--fragment_sd ${fragment_sd} \
${file_outq_1} ${ID} "${output}/final" &
fi;
fi; 



# When a new job is created, the program updates the $CURRENT_POOL_SIZE variable before next iteration
	CURRENT_POOL_SIZE=$(jobs | wc -l)
	echo $CURRENT_POOL_SIZE
###only 50 processes simultaneously
	while [[ $CURRENT_POOL_SIZE -ge $MAX_POOL_SIZE ]]; do    
    # The above "echo" and "sleep" is for demo purposes only.
    # In a real usecase, remove those two and keep only the following line
    # It will drastically increase the performance of the script
		CURRENT_POOL_SIZE=$(jobs | wc -l)

done
done
wait 
echo "TRACER done"

echo "Completed successfully."
