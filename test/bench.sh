#!/bin/bash

set -e


SQ_MODEL=dna-r10-prom
ONT_MODEL_PREFIX=dna_r10.4.1_e8.2_400bps
IDENTITY=/install/biorand/bin/identitydna.sh
GENOME=/genome/hg38noAlt.fa
GENOME_IDX=/genome/hg38noAlt.idx
#CUSTOM_MODEL="--kmer-model r10.4.1_260bps.nucleotide.9mer.model" #uncomment this line if you want to use a custom model


SQUIGULATOR=/home/hasindu/hasindu2008.git/squigulator/squigulator
BUTTERY_EEL_VENV=/home/hasindu/hasindu2008.git/buttery-eel/venv3-guppy-6.4.2
GUPPY_BIN=/install/ont-guppy-6.4.2/bin/
CUDA_DEVICE='cuda:0'



eel-solo(){
#from https://unix.stackexchange.com/questions/55913/whats-the-easiest-way-to-find-an-unused-local-port
PORT=$(netstat -aln | awk '
  $6 == "LISTEN" {
    if ($4 ~ "[.:][0-9]+$") {
      split($4, a, /[:.]/);
      port = a[length(a)];
      p[port] = 1
    }
  }
  END {
    for (i = 5000; i < 65000 && p[i]; i++){};
    if (i == 65000) {exit 1};
    print i
  }
  ') && source ${BUTTERY_EEL_VENV}/bin/activate && buttery-eel -g ${GUPPY_BIN} "$@" --port ${PORT} --use_tcp && deactivate
}

export -f eel-solo

RUN (){

MODEL_FSP=$1
ONT_MODEL=${ONT_MODEL_PREFIX}"_"${MODEL_FSP}".cfg"

echo -e "dwell_mean\tsample\tmean\tsstdev\tq1\tmedian\tq3\tn" > dwell_mean_vs_acc_${MODEL_FSP}.txt
for i in 5 6 7 8 9 10 11 12 13 14 15
do
echo -e -n "$i\t" >> dwell_mean_vs_acc_${MODEL_FSP}.txt;
${SQUIGULATOR}  -x ${SQ_MODEL}  ${GENOME} -o a.blow5 --seed 100 --dwell-std 0.0 --dwell-mean $i --amp-noise 1.0 -n 4000 ${CUSTOM_MODEL}
eel-solo  --config ${ONT_MODEL}  --device ${CUDA_DEVICE} -i a.blow5 -o  a.fastq
${IDENTITY} ${GENOME_IDX} a.fastq | tail -1  >> dwell_mean_vs_acc_${MODEL_FSP}.txt
done


echo -e "dwell_std\tsample\tmean\tsstdev\tq1\tmedian\tq3\tn" > dwell_std_vs_acc_${MODEL_FSP}.txt
for i in 0 1 2 3 4 5 6
do
echo -e -n "$i\t" >> dwell_std_vs_acc_${MODEL_FSP}.txt;
${SQUIGULATOR}  -x ${SQ_MODEL}  ${GENOME} -o a.blow5 --seed 100 --dwell-std $i --dwell-mean 10.0 --amp-noise 1.0 -n 4000 ${CUSTOM_MODEL}
eel-solo  --config ${ONT_MODEL}  --device ${CUDA_DEVICE} -i a.blow5 -o  a.fastq
${IDENTITY} ${GENOME_IDX} a.fastq | tail -1  >> dwell_std_vs_acc_${MODEL_FSP}.txt
done


echo -e "amp_noise\tsample\tmean\tsstdev\tq1\tmedian\tq3\tn" > amp_noise_vs_acc_${MODEL_FSP}.txt
for i in 0 0.25 0.5 0.75 1 1.5 2
do
echo -e -n "$i\t" >> amp_noise_vs_acc_${MODEL_FSP}.txt;
${SQUIGULATOR}  -x ${SQ_MODEL}  ${GENOME} -o a.blow5 --seed 100 --dwell-std 0 --dwell-mean 10.0 --amp-noise $i -n 4000 ${CUSTOM_MODEL}
eel-solo  --config ${ONT_MODEL}  --device ${CUDA_DEVICE} -i a.blow5 -o  a.fastq
${IDENTITY} ${GENOME_IDX} a.fastq | tail -1  >> amp_noise_vs_acc_${MODEL_FSP}.txt
done

}

RUN fast
RUN hac
RUN sup
