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
/opt/transrate-1.0.2-linux-x86_64/transrate -o tuco_1 -t 8 \
-a /home/molly/tucoKidney/trinity_v2.1.1/Rcorr_trinity_tucoKidney.Trinity.fasta \
--left /home/molly/tucoKidney/Rcorrect/tuco_kidney.1.cor.fq \
--right /home/molly/tucoKidney/Rcorrect/tuco_kidney.2.cor.fq  
```
### Evaluate Original Assembly Completeness with BUSCO (v1.1b1)
```
python3 /opt/BUSCO_v1.1b1/BUSCO_v1.1b1.py -m trans --cpu 10 -l /opt/BUSCO_v1.1b1/vertebrata \
-o tuco_1 -in /home/molly/tucoKidney/trinity_v2.1.1/Rcorr_trinity_tucoKidney.Trinity.fasta
```
### Evaluate Optimized Assembly Completeness with BUSCO (v1.1b1)
```
python3 /opt/BUSCO_v1.1b1/BUSCO_v1.1b1.py -m trans --cpu 10 -l /opt/BUSCO_v1.1b1/vertebrata \
-o tuco_1.1 -in /home/molly/tucoKidney/transrate_v1.0.2/tuco_1/Rcorr_trinity_tucoKidney.Trinity/good.Rcorr_trinity_tucoKidney.Trinity.fasta
```
##Continue Pipeline with Optimized Assembly
### Estimate Gene Expression with Kallisto (v0.42.4) 
```
/usr/local/bin/kallisto index -i kallisto.idx /home/molly/tucoKidney/transrate_v1.0.2/tuco_1/Rcorr_trinity_tucoKidney.Trinity/good.Rcorr_trinity_tucoKidney.Trinity.fasta
/usr/local/bin/kallisto quant -t 10 -i kallisto.idx -o kallisto_orig /home/molly/tucoKidney/Rcorrect/tuco_kidney.1.cor.fq /home/molly/tucoKidney/Rcorrect/tuco_kidney.2.cor.fq
```
### Estimate Gene Expression with Salmon (v0.3.0)
```
/opt/salmon-0.5.1/bin/salmon index -t /home/molly/tucoKidney/transrate_v1.0.2/tuco_1/Rcorr_trinity_tucoKidney.Trinity/good.Rcorr_trinity_tucoKidney.Trinity.fasta -i salmon.idx --type quasi -k 31
/opt/salmon-0.5.1/bin/salmon quant -p 32 -i salmon.idx -l IU -1 /home/molly/tucoKidney/Rcorrect/tuco_kidney.1.cor.fq -2 /home/molly/tucoKidney/Rcorrect/tuco_kidney.2.cor.fq -o salmon_orig
```
### Filter out Contigs Based on Gene Expression < TPM=1 
```
awk '1>$5{next}1' /home/molly/tucoKidney/kallisto_v0.42.4/kallisto_1/kallisto_orig/abundance.tsv | awk '{print $1}' > kallist
awk '1>$3{next}1' /home/molly/tucoKidney/salmon_v0.3.0/salmon_1/salmon_orig/quant.sf | sed  '1,10d' | awk '{print $1}' > salist
cat kallist salist | sort -u > uniq_list
sed -i ':begin;N;/[ACTGNn-]\n[ACTGNn-]/s/\n//;tbegin;P;D' /home/molly/tucoKidney/transrate_v1.0.2/tuco_1/Rcorr_trinity_tucoKidney.Trinity/good.Rcorr_trinity_tucoKidney.Trinity.fasta

for i in $(cat uniq_list);
   do grep --no-group-separator --max-count=1 -A1 -w $i /home/molly/tucoKidney/transrate_v1.0.2/tuco_1/Rcorr_trinity_tucoKidney.Trinity/good.Rcorr_trinity_tucoKidney.Trinity.fasta >> Rcorr_highexp.trinity.Trinity.fasta; done
```

### Generate 2nd Optimized Assembly Including Quality Check with Transrate (v1.0.1)
```
/opt/transrate-1.0.1-linux-x86_64/transrate -o tuco_2 -t 8 \
-a /home/molly/tucoKidney/TPM/TPM_1/Rcorr_highexp.trinity.Trinity.fasta \
--left /home/molly/tucoKidney/Rcorrect/tuco_kidney.1.cor.fq \
--right /home/molly/tucoKidney/Rcorrect/tuco_kidney.2.cor.fq  

```

### Evaluate TPM Filtration Assembly Completeness with BUSCO (v1.1b1)
```
python3 /opt/BUSCO_v1.1b1/BUSCO_v1.1b1.py -m trans --cpu 10 -l /opt/BUSCO_v1.1b1/vertebrata \
-o tuco_2 -in /home/molly/tucoKidney/TPM/TPM_1/Rcorr_highexp.trinity.Trinity.fasta 
```
### Evaluate 2nd Optimized Assembly Completeness with BUSCO (v1.1b1)
```
python3 /opt/BUSCO_v1.1b1/BUSCO_v1.1b1.py -m trans --cpu 10 -l /opt/BUSCO_v1.1b1/vertebrata \
-o tuco_2.1 -in /home/molly/tucoKidney/transrate_v1.0.2/tuco_2/Rcorr_highexp.trinity.Trinity/good.Rcorr_highexp.trinity.Trinity.fasta
```
## Continue Pipeline with Original Assembly
### Estimate Gene Expression with Kallisto (v0.42.4)
```
/usr/local/bin/kallisto index -i kallisto.idx /home/molly/tucoKidney/trinity_v2.1.1/Rcorr_trinity_tucoKidney.Trinity.fasta
/usr/local/bin/kallisto quant -t 10 -i kallisto.idx -o kallisto_orig /home/molly/tucoKidney/Rcorrect/tuco_kidney.1.cor.fq /home/molly/tucoKidney/Rcorrect/tuco_kidney.2.cor.fq
```
### Estimate Gene Expression with Salmon (v0.3.0)
```
/opt/salmon-0.5.1/bin/salmon index -t /home/molly/tucoKidney/trinity_v2.1.1/Rcorr_trinity_tucoKidney.Trinity.fasta -i salmon.idx --type quasi -k 31
/opt/salmon-0.5.1/bin/salmon quant -p 32 -i salmon.idx -l IU -1 /home/molly/tucoKidney/Rcorrect/tuco_kidney.1.cor.fq -2 /home/molly/tucoKidney/Rcorrect/tuco_kidney.2.cor.fq -o salmon_orig
```
### Filter out Contigs Based on Gene Expression < TPM=1
```
awk '1>$5{next}1' /home/molly/tucoKidney/kallisto_v0.42.4/kallisto_2/kallisto_orig/abundance.tsv | awk '{print $1}' > kallist
awk '1>$3{next}1' /home/molly/tucoKidney/salmon_v0.3.0/salmon_2/salmon_orig/quant.sf | sed  '1,10d' | awk '{print $1}' > salist
cat kallist salist | sort -u > uniq_list
sed -i ':begin;N;/[ACTGNn-]\n[ACTGNn-]/s/\n//;tbegin;P;D' /home/molly/tucoKidney/trinity_v2.1.1/Rcorr_trinity_tucoKidney.Trinity.fasta

for i in $(cat uniq_list);
   do grep --no-group-separator --max-count=1 -A1 -w $i /home/molly/tucoKidney/trinity_v2.1.1/Rcorr_trinity_tucoKidney.Trinity.fasta >> Rcorr_highexp.trinity.Trinity.fasta; done
```

### Generate Optimized Assembly Including Quality Check with Transrate (v1.0.1)
```
/opt/transrate-1.0.1-linux-x86_64/transrate -o tuco_3 -t 10 \
-a /home/molly/tucoKidney/TPM/TPM_2/Rcorr_highexp.trinity.Trinity.fasta \
--left /home/molly/tucoKidney/Rcorrect/tuco_kidney.1.cor.fq \
--right /home/molly/tucoKidney/Rcorrect/tuco_kidney.2.cor.fq 
```

### Evaluate TPM Filtration Assembly Completeness with BUSCO (v1.1b1)
```
python3 /opt/BUSCO_v1.1b1/BUSCO_v1.1b1.py -m trans --cpu 10 -l /opt/BUSCO_v1.1b1/vertebrata \
-o tuco_3 -in /home/molly/tucoKidney/TPM/TPM_2/Rcorr_highexp.trinity.Trinity.fasta 
```

### Evaluate Optimized Assembly Comleteness with BUSCO (v1.1b1)
```
python3 /opt/BUSCO_v1.1b1/BUSCO_v1.1b1.py -m trans --cpu 10 -l /opt/BUSCO_v1.1b1/vertebrata \
-o tuco_3.1 -in /home/molly/tucoKidney/transrate_v1.0.2/tuco_3/Rcorr_highexp.trinity.Trinity/good.Rcorr_highexp.trinity.Trinity.fasta
```

### Compare all Metrics and Scores to Determine Best Assembly

### Annotate with dammit (v0.3)
```
dammit annotate /home/molly/tucoKidney/transrate_v1.0.2/tuco_1/Rcorr_trinity_tucoKidney.Trinity/good.Rcorr_trinity_tucoKidney.Trinity.fasta \
--user-databases /data/dammit_databases/tcdb.fasta  \
--busco-group vertebrata \
--n_threads 12 \
--full \
--database-dir /data/dammit_databases/
```
