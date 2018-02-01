In this analysis, we will need the **.bed** files that have been called by either `MACS1` or `MACS2` programs and stitch the peaks in a window specified by ourselves, say, a 4-kb window. In my current analysis, I need to get the H3K4me3 peaks within a 4-kb window to be merged together. We use `BEDOPS` [BEDOPS site](http://bedops.readthedocs.io/en/latest/index.html) to carry some of the work.  

Since my peak files are not ended as ".bed", I just rename them by adding ".bed":

```sh
cd bedfiles

for f in *.peaks; do; mv ${f} ${f}.bed; done
``` 

Merge the peaks in a 4-kb window by `bedtools merge`(here I use a job script for submitting the job to the cluster - it can be run on local machine; it should be fast):  

`merge_broaddomain.sh`

```sh
#! /bin/bash
 
set -e 
set -u
set -o pipefail
 
root='pwd'
mkdir lsf
mkdir logs
mkdir H3K4me3_BroadPeaks_merged_4kb
 
cat files.txt | while read f
do
prefix=$(basename ${f} .bed)
 
JobString="
#BSUB -J ${prefix}_merge
#BSUB -W 12:00
#BSUB -o /rsrch2/headneck/zliu11/H3K4me3_broadomain_analysis/logs     
#BSUB -e /rsrch2/headneck/zliu11/H3K4me3_broadomain_analysis/logs
#BSUB -q medium
#BSUB -n 10
#BSUB -M 32768 
#BSUB -R rusage[mem=32768]
#BSUB -N
#BSUB -u zliu11@MDAnderson.org
 
mergeBed -i /rsrch2/headneck/zliu11/H3K4me3_broadomain_analysis/bedfiles/${f} -d 4000 > /rsrch2/headneck/zliu11/H3K4me3_broadomain_analysis/H3K4me3_BroadPeaks_merged_4kb/${prefix}.merged.bed
"
echo "$JobString" > /rsrch2/headneck/zliu11/H3K4me3_broadomain_analysis/lsf/${prefix}_mergebed.lsf
bsub < /rsrch2/headneck/zliu11/H3K4me3_broadomain_analysis/lsf/${prefix}_mergebed.lsf
done
```
Now have a look at the merged bed files:

```sh
less HNSC-U4XX_E_vs_HNSC-U4XX_G_macs2_peaks.broadPeak.merged.bed 

chr1    860400  860619
chr1    934709  937023
chr1    955109  957224
chr1    994701  995424
chr1    1093089 1093331
chr1    1098334 1099942
chr1    1208493 1209755
chr1    1241834 1243208
chr1    1259908 1261107
chr1    1283174 1283856
chr1    1310235 1310510
chr1    1333928 1335767
chr1    1341761 1343248
chr1    1355613 1356587
chr1    1509775 1510164
chr1    1535163 1535980
chr1    1550451 1551852
chr1    1709320 1710295
chr1    1821687 1822386
chr1    2125683 2126459
chr1    2158531 2161500
chr1    2343315 2345639
chr1    2456910 2457523
chr1    2478120 2478483
```

Run  

`find *bed | parallel -j 6 './add_name.sh {}' `  

in the console to feed the bed file names into {} for the following codes to extact the sample id and put the id in the 4th field of the bed file. Here the unique id is extracted from the bed file name matching the pattern of `s/_E_vs_.+//`. For example, my bed file name is like `HNSC-U17A_E_vs_HNSC-U17A_G_macs2_peaks.broadPeak.bed`. The sample id is extracted by `awk` shown in the following code. The regular expression should be changed accordingly based on the naming of the bed files.

`add_name.sh`

```sh
set -e
set -u
set -o pipefail

file=$1

fileName=$(echo $1 | sed -E 's/_E_vs_.+//')

cat $file | awk -v id="$fileName" '$3=$3"\t"id"_peak"NR'  OFS="\t" > ${fileName}_merged_H3K4me3_peaks.bed
```

```sh
less HNSC-U4XX_merged_H3K4me3_peaks.bed

chr1    860400  860619  HNSC-U4XX_peak1
chr1    934709  937023  HNSC-U4XX_peak2
chr1    955109  957224  HNSC-U4XX_peak3
chr1    994701  995424  HNSC-U4XX_peak4
chr1    1093089 1093331 HNSC-U4XX_peak5
chr1    1098334 1099942 HNSC-U4XX_peak6
chr1    1208493 1209755 HNSC-U4XX_peak7
chr1    1241834 1243208 HNSC-U4XX_peak8
chr1    1259908 1261107 HNSC-U4XX_peak9
chr1    1283174 1283856 HNSC-U4XX_peak10
chr1    1310235 1310510 HNSC-U4XX_peak11
chr1    1333928 1335767 HNSC-U4XX_peak12
chr1    1341761 1343248 HNSC-U4XX_peak13
chr1    1355613 1356587 HNSC-U4XX_peak14
chr1    1509775 1510164 HNSC-U4XX_peak15
chr1    1535163 1535980 HNSC-U4XX_peak16
chr1    1550451 1551852 HNSC-U4XX_peak17
chr1    1709320 1710295 HNSC-U4XX_peak18
chr1    1821687 1822386 HNSC-U4XX_peak19
chr1    2125683 2126459 HNSC-U4XX_peak20
chr1    2158531 2161500 HNSC-U4XX_peak21
chr1    2343315 2345639 HNSC-U4XX_peak22
chr1    2456910 2457523 HNSC-U4XX_peak23
chr1    2478120 2478483 HNSC-U4XX_peak24
chr1    2517055 2518306 HNSC-U4XX_peak25
chr1    2573801 2574474 HNSC-U4XX_peak26
```

Run the following to strip peak numbers and keep only the id numbers:

```sh
cat  *bed | sort-bed - |  bedmap --echo --echo-map-id-uniq --delim '\t' -  | gawk -v OFS="\t" '{a=gensub(/_peak[0-9]+/, "", "g",  $5); $5=a; print}' > merged_super.bed
```

Run this to count the unique sample id numbers and put it in the 6th column (require `gawk`):

```sh
while read line;do cnt=`echo "$line" | cut -f5 | tr ";" "\n" | sort -u | wc -l`;echo -e "$line\t$cnt";done < merged_super.bed > merged_super_clean.bed
```
But the above codes are slow in the linux command line. It is better to run the R script as the following:

```sh
Rscript counting_uniq_mergeBed_regions.R
```
The codes for `counting_uniq_mergeBed_regions.R`:

```r
library(stringr)
library(tidyverse)

bed<- read_tsv("merged_super.bed", col_names = F)

bed %>% mutate(X6 = unlist(lapply(str_split(X5, ";"), function(x) length(unique(x))))) %>% write.table("merged_super_clean.bed", quote =F, sep = "\t", row.names =F, col.names =F)
```

```sh
head merged_super_clean.bed

chr1	713372	714755	HNSC-U14B_peak1	HNSC-M19L;HNSC-U14B	2
chr1	713987	715265	HNSC-M19L_peak1	HNSC-M19L;HNSC-U14B	2
chr1	838705	840598	HNSC-M19L_peak2	HNSC-183X;HNSC-JH22;HNSC-JHU2;HNSC-M19L;HNSC-MD8L;HNSC-U14B	6
chr1	839204	840961	HNSC-MD8L_peak1	HNSC-183X;HNSC-JH22;HNSC-JHU2;HNSC-M19L;HNSC-MD8L;HNSC-U14B	6
chr1	839249	840848	HNSC-U14B_peak2	HNSC-183X;HNSC-JH22;HNSC-JHU2;HNSC-M19L;HNSC-MD8L;HNSC-U14B	6
chr1	839511	839837	HNSC-183X_peak1	HNSC-183X;HNSC-JH22;HNSC-M19L;HNSC-MD8L;HNSC-U14B	5
chr1	839536	840385	HNSC-JH22_peak1	HNSC-183X;HNSC-JH22;HNSC-JHU2;HNSC-M19L;HNSC-MD8L;HNSC-U14B	6
chr1	840155	840618	HNSC-JHU2_peak1	HNSC-JH22;HNSC-JHU2;HNSC-M19L;HNSC-MD8L;HNSC-U14B	5
chr1	855772	861347	HNSC-JH22_peak2	HNSC-183X;HNSC-584A;HNSC-HN4X;HNSC-JH22;HNSC-JHU2;HNSC-M19L;HNSC-MD8L;HNSC-U14A;HNSC-U14B;HNSC-U19X;HNSC-U4XX	11
chr1	856828	863163	HNSC-U14B_peak3	HNSC-183X;HNSC-584A;HNSC-HN4X;HNSC-JH22;HNSC-JHU2;HNSC-M19L;HNSC-MD8L;HNSC-U14A;HNSC-U14B;HNSC-U19X;HNSC-U4XX	11
```

Finally, run the following to select the number of samples that contain a unique merged region - the number is specified in `$6 >= `. The unique broad domain H3K4me3 peaks are then merged to give rise to a final `bed` file for a group of samples with a particular biological feature(s).

```sh
cat merged_super_clean.bed | awk '$6 >= 17' | sort-bed - | bedtools merge -i - > high_p53_H3K4me3_broad.bed
```