## the mapping strategy of this project

cd /home/wanghao/data1/DATA/ZXseq2/rawdata

set=*2.gz

for sample in ${set[*]}

do

hisat2 -p 16 --no-softclip --dta-cufflinks --rna-strandness F -x /home/wanghao/Refgenome/hg38_tran -U /home/wanghao/data1/DATA/ZXseq2/rawdata/$sample -S $sample.FWD.sam 2>$sample.FWD.log

done
