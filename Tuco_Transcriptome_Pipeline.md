## Tuco Transcriptome Assembly Pipeline 

### Initial Quality Check with SolexaQA++ (v3.1.4)

```
SolexaQA++ analysis tuco_kidney.1.fq tuco_kidney.2.fq
```
### Error Correct with Rcorrector (v1.01)
```
perl /home/molly/bin/Rcorrector/run_rcorrector.pl -k 31 -t 4 \
-1 ../rawdata/tuco_kidney.1.fq -2 ../rawdata/tuco_kidney.2.fq
```
### Transcriptome Assembly with Trinity (v2.1.1) 
```
Trinity --seqType fq --max_memory 40G --trimmomatic --CPU 10 --full_cleanup --output Rcorr_trinity_tucoKidney \
--left /home/molly/tucoKidney/Rcorrect/tuco_kidney.1.cor.fq \
--right /home/molly/tucoKidney/Rcorrect/tuco_kidney.2.cor.fq \
--quality_trimming_params "ILLUMINACLIP:/opt/trinityrnaseq/trinity-plugins/Trimmomatic/adapters/TruSeq3-PE-2.fa:2:40:15 LEADING:2 TRAILING:2 MINLEN:25"
``` 
### Optimized Assembly Including Quality Check with Transrate (v)
```
transrate -o assemb_name -t 10 \
-a Rcorr_trinity_tucoKidney.Trinity.fasta \
--left /home/molly/tucoKidney/trinity_v2.1.1/Rcorr_trinity_tucoKidney/ \
--right /home/molly/tucoKidney/trinity_v2.1.1/Rcorr_trinity_tucoKidney/ 
```
### Evaluate Original and Optimized Assembly Completeness with BUSCO (v)