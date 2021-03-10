####cut adapt for raw reads

## open the raw reads dir
cd /home/wanghao/data1/DATA/ZXseq2/P101SC17051084-01-B39-21_data

set=*_1.fq.gz

for sample in ${set[*]}

do

base=$(basename "$sample" _1.fq.gz)

fq_1=$base"_1.fq.gz"
fq_2=$base"_2.fq.gz"

## cut the adapt from 3' and 5'
## save those clean reads in a new dir
cutadapt -u 2 -a N{12}AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -m 15 -o out.1.fastq -p out.2.fastq ${fq_1} ${fq_2}
cutadapt -u 12 -a N{2}AAAAAAAAAAAAAAAAAAAGATCGGAAGAGCGTCGTGTAGGG -m 15 -o ~/data1/DATA/ZXseq2/rawdata/$base.fq.2.gz -p ~/data1/DATA/ZXseq2/rawdata/$base.fq.1.gz out.2.fastq out.1.fastq
rm out.2.fastq out.1.fastq

done
