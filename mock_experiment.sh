timing_file="/home/bg31/firstRaceProject/DiversitySampling/results/experiment_timing.txt"
directorypath="/home/bg31/firstRaceProject/DiversitySampling/part_"
race="/home/bg31/firstRaceProject/DiversitySampling/bin/sampleracesavable"
taus="1,1.7,2.8,4.5,7.7,12.9,21.5,35.9,59.9,100"
outputs="/home/bg31/firstRaceProject/DiversitySampling/results/experiment_sample_1 /home/bg31/firstRaceProject/DiversitySampling/results/experiment_sample_2"

for dir in 0; do
  echo working on part ${dir}
  cd ${directorypath}${dir}
  mkdir temp
  find . -name '*.gz' -exec cp {} temp/{}\;
  cd temp
  gunzip *
  for f1 in *_1.fastq; do
    f2=${f1:0:13}_2.fastq
    (time ${race} ${taus} PE experimentsavefile.bin ${f1} ${f2} ${outputs} --range 10 --k 15) >> ${timing file}
  done
  cd ..
  rm -r temp
done