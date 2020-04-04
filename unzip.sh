timing_file="/home/Users/bg31/DiversitySampling/results/experiment_timing.txt"
directorypath="/home/public_data/ibdmdb/part_"
race="/home/Users/bg31/DiversitySampling/bin/sampleracesavable"
temp="/home/Users/bg31/DiversitySampling/unzipped"
savefile="/home/Users/bg31/DiversitySampling/results/experimentsavefile.bin"
taus="1,1.7,2.8,4.5,7.7,12.9,21.5,35.9,59.9,100"
outputs="/home/Users/bg31/DiversitySampling/results/experiment_sample_1 /home/Users/bg31/DiversitySampling/results/experiment_sample_2"

mkdir ${temp}
for dir in 0 1 2 3 4 5 6 7 8 9; do
    echo working on part ${dir}
    cd ${directorypath}${dir}
    for z1 in *_1_reads.fq.gz; do
        echo ${z1}
        cp ${z1} ${temp}/${z1}
        gunzip ${temp}/${z1} &
        z2=${z1:0:10}_2_reads.fq.gz
        echo ${z2}
        cp ${z2} ${temp}/${z2}
        gunzip ${temp}/${z2} &
    done
done

echo working on part chun
cd /home/public_data/hmp2/chunxiao_download_2_27_2020_hmp2
for z1 in *_1_reads.fq.gz; do
    echo ${z1}
    cp ${z1} ${temp}/${z1}
    gunzip ${temp}/${z1} &
    z2=${z1:0:10}_2_reads.fq.gz
    echo ${z2}
    cp ${z2} ${temp}/${z2}
    gunzip ${temp}/${z2} &
done