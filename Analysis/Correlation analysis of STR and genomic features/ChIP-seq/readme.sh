#!/bin/bash

WORKDIR=/home2/niuyw/project/STR/genomic_features_1Mb/ChIP-seq

ln -s ../../genomic_features/ChIP-seq/*.bed.gz .

# count number per 100-kb window
for i in `ls *.bed.gz`
do
  j=$(basename $i .bed.gz)
  zcat $i | \
    sort-bed - | \
    bedmap --echo --count --delim '\t' ../hg38.autosomes.bin_1M.rmCen.bed - > $j.bin_count
done

# count mean
ln -s ENCFF160CEI.bin_count H2AK5ac.count
ln -s ENCFF808BIC.bin_count H2BK120ac.count
ln -s ENCFF315NAV.bin_count H2BK12ac.count
ln -s ENCFF997ZDN.bin_count H2BK15ac.count
ln -s ENCFF628KOM.bin_count H2BK20ac.count
ln -s ENCFF068IZP.bin_count H2BK5ac.count
ln -s ENCFF642IBN.bin_count H3K14ac.count
ln -s ENCFF397EFT.bin_count H3K18ac.count
ln -s ENCFF851FNE.bin_count H3K23ac.count
ln -s ENCFF987CBQ.bin_count H3K23me2.count
ln -s ENCFF711LQB.bin_count H3K4ac.count
ln -s ENCFF441MSJ.bin_count H3K4me2.count
ln -s ENCFF053HKF.bin_count H3K56ac.count
ln -s ENCFF385DSV.bin_count H3K79me1.count
ln -s ENCFF202EZL.bin_count H3K79me2.count
ln -s ENCFF237OAY.bin_count H4K5ac.count
ln -s ENCFF760EFQ.bin_count H4K8ac.count
ln -s ENCFF964FVB.bin_count H4K91ac.count
ln -s ENCFF821AQO.bin_count CTCF.count
ln -s ENCFF418QVJ.bin_count POLR2AphosphoS5.count
ln -s ENCFF422HDN.bin_count POLR2A.count
ln -s ENCFF834UVX.bin_count EP300.count
ln -s ENCFF439CWL.bin_count H4K20me1.count

# H3K27ac
paste -d '\t' ENCFF162HPV.bin_count ENCFF045CUG.bin_count ENCFF897KUI.bin_count | \
  awk 'BEGIN{FS=OFS="\t"}{sum=$4+$8+$12;print $1,$2,$3,sum/3}' > H3K27ac.count

# H3K27me3
paste -d '\t' ENCFF296RYM.bin_count ENCFF156RHD.bin_count ENCFF760MIC.bin_count ENCFF084QDP.bin_count ENCFF411ESN.bin_count ENCFF596UYU.bin_count | \
  awk 'BEGIN{FS=OFS="\t"}{sum=$4+$8+$12+$16+$20+$24;print $1,$2,$3,sum/6}' > H3K27me3.count

# H3K36me3
paste -d '\t' ENCFF813VFV.bin_count ENCFF093YFN.bin_count ENCFF599IHW.bin_count ENCFF736WAN.bin_count ENCFF422MTS.bin_count ENCFF238MDY.bin_count | \
  awk 'BEGIN{FS=OFS="\t"}{sum=$4+$8+$12+$16+$20+$24;print $1,$2,$3,sum/6}' > H3K36me3.count

# H3K4me1
paste -d '\t' ENCFF238YJA.bin_count ENCFF694ENI.bin_count ENCFF711NNX.bin_count ENCFF558IKG.bin_count | \
  awk 'BEGIN{FS=OFS="\t"}{sum=$4+$8+$12+$16;print $1,$2,$3,sum/4}' > H3K4me1.count

# H3K4me3
paste -d '\t' ENCFF744ORJ.bin_count ENCFF047VRR.bin_count ENCFF277AOQ.bin_count ENCFF436DIF.bin_count ENCFF408FCY.bin_count ENCFF456NIF.bin_count | \
  awk 'BEGIN{FS=OFS="\t"}{sum=$4+$8+$12+$16+$20+$24;print $1,$2,$3,sum/6}' > H3K4me3.count

# H3K9ac
paste -d '\t' ENCFF806NNI.bin_count ENCFF343GTP.bin_count | \
  awk 'BEGIN{FS=OFS="\t"}{sum=$4+$8;print $1,$2,$3,sum/2}' > H3K9ac.count

# H3K9me3
paste -d '\t' ENCFF654ZZO.bin_count ENCFF348GGB.bin_count ENCFF250GSY.bin_count ENCFF378HMV.bin_count | \
  awk 'BEGIN{FS=OFS="\t"}{sum=$4+$8+$12+$16;print $1,$2,$3,sum/4}' > H3K9me3.count





