timing_file="/home/Users/bg31/DiversitySampling/results/experiment_timing.txt"
directorypath="/home/public_data/ibdmdb/part_"
race="/home/Users/bg31/DiversitySampling/bin/sampleracesavable"
temp="/home/Users/bg31/DiversitySampling/temp"
savefile="/home/Users/bg31/DiversitySampling/results/experimentsavefile"
taus="1,1.7,2.8,4.5,7.7,12.9,21.5,35.9,59.9,100"
outputs="/home/Users/bg31/DiversitySampling/results/experiment_sample_1 /home/Users/bg31/DiversitySampling/results/experiment_sample_2"

for dir in 0 1 2 3 4 5 6 7 8 9; do
    echo working on part ${dir}
    cd ${directorypath}${dir}
    mkdir ${temp}
    find . -maxdepth 1 -name '*.gz' - exec cp {} ${temp}/{} \;
    cd ${temp}
    gunzip *
    for f1 in *_1_reads.fq; do
        echo ${f1}
        f2=${f1:0:10}_2_reads.fq
        mytime="$(time ( ${race} ${taus} PE ${savefile} ${f1} ${f2} ${outputs} --range 500000 --k 15 ) 2>&1 1>/dev/null )"
        echo ${mytime} >> ${timing_file}
    done
    cd ..
    rm -r ${temp}
done

# needa add the chun thing