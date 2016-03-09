## Tuco-tuco Transcriptome Assembly Pipeline 

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
### Generate Optimized Assembly Including Quality Check with Transrate (v1.0.2)
```
/opt/transrate-1.02-linux-x86_64 -o tuco_2 -t 8 \
-a /home/molly/tucoKidney/trinity_v2.1.1/Rcorr_trinity_tucoKidney.Trinity.fasta \
--left /home/molly/tucoKidney/Rcorrect/tuco_kidney.1.cor.fq \
--right /home/molly/tucoKidney/Rcorrect/tuco_kidney.2.cor.fq  
```
### Evaluate Original Assembly Completeness with BUSCO (v1.1b1)
```
python3 /opt/BUSCO_v1.1b1/BUSCO_v1.1b1.py -m trans --cpu 10 -l /opt/BUSCO_v1.1b1/vertebrata \
-o tuco_2 -in Rcorr_trinity_tucoKidney.Trinity.fasta
```
### Evaluate Optimized Assembly Completeness with BUSCO (v1.1b1)
```
python3 /opt/BUSCO_v1.1b1/BUSCO_v1.1b1.py -m trans --cpu 10 -l /opt/BUSCO_v1.1b1/vertebrata \
-o tuco_2.1 -in good.Rcorr_trinity_tucoKidney.Trinity.fasta
```
### Filter and Estimate Expression with Kallisto (v) 
```
kallisto index -i kallisto.idx Rcorr_trinity_tucoKidney.Trinity.fasta
kallisto quant -t 10 -i kallisto.idx -o kallisto_orig /home/molly/tucoKidney/Rcorrect/tuco_kidney.1.cor.fq /home/molly/tucoKidney/Rcorrect/tuco_kidney.2.cor.fq
```
### Filter and Estimate Expression with Salmon (v)
```
~/salmon-0.5.1/bin/salmon index -t Rcorr_trinity_tucoKidney.Trinity.fasta -i salmon.idx --type quasi -k 31
~/salmon-0.5.1/bin/salmon quant -p 32 -i salmon.idx -l IU -1 /home/molly/tucoKidney/Rcorrect/tuco_kidney.1.cor.fq -2 /home/molly/tucoKidney/Rcorrect/tuco_kidney.2.cor.fq -o salmon_orig
```
### Annotate with dammit!
```
mkdir /mnt/dammit/ && cd /mnt/dammit
dammit databases --install --database-dir /mnt/dammit --full --busco-group metazoa
dammit annotate assembly.fasta --busco-group metazoa --n_threads 10 --database-dir /mnt/dammit/ --full
```
