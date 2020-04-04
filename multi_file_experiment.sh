timing_file="./results/multifile_experiment_timing.txt"
race="./bin/sampleracemanyfilestaus"
temp="./unzipped"
taus="1,1.7,2.8,4.5,7.7,12.9,21.5,35.9,59.9,100"
outputs="./results/multifile_experiment_sample_1 ./multifile_experiment_sample_2"
inputlistfromrepo="./inputlist.txt"
inputlistfromtemp="../inputlist.txt"

# Make file for file list

cd ${temp}
pwd > ${inputlistfromtemp}
cd ..
ls ${temp} >> ${inputlistfromrepo}
mytime="$(time ( ${race} ${taus} PE ${inputlistfromrepo} ${outputs} --range 500000 --k 15 ) 2>&1 1>/dev/null )"
echo ${mytime} >> ${timing_file}