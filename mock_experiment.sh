timing_file="/home/bg31/firstRaceProject/DiversitySampling/results/experiment_timing.txt"
directorypath="/home/bg31/firstRaceProject/DiversitySampling/part_"
race="/home/bg31/firstRaceProject/DiversitySampling/bin/sampleracesavable"
savefile="/home/bg31/firstRaceProject/DiversitySampling/experimentsavefile.bin"
taus="1,1.7,2.8,4.5,7.7,12.9,21.5,35.9,59.9,100"
outputs="/home/bg31/firstRaceProject/DiversitySampling/results/experiment_sample_1 /home/bg31/firstRaceProject/DiversitySampling/results/experiment_sample_2"

for dir in 0; do
  echo working on part ${dir}
  cd ${directorypath}${dir}
  pwd  
  mkdir temp
  find . -maxdepth 1 -name '*.gz' -exec cp {} temp/{} \;
  cd temp
  pwd
  ls
  gunzip *
  for f1 in *_1.fastq; do
    f2=${f1:0:13}_2.fastq
    mytime="$(time ( ${race} ${taus} PE ${savefile} ${f1} ${f2} ${outputs} --range 10 --k 15 ))"
    echo ${mytime} >> ${timing_file}
  done
  cd ..
  pwd
  rm -r temp
  ls
done
