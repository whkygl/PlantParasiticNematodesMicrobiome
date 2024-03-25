    db=/d/nematode_microbiome_2023-12-05/pipline/amplicon/PairedEnd_sequence_merging/public
    wd=/d/nematode_microbiome_2023-12-05/pipline/amplicon/PairedEnd_sequence_merging
    export PATH=$PATH:${db}/win
    cd ${wd}
    mkdir -p seq
    ls -sh seq
    mkdir -p temp
    cat -A result/metadata.txt | head -n3
    cat -A result/metadata2.txt | head -n3
    gunzip seq/*.gz
    gunzip seq2/*.gz
    ls -sh seq/
    head -n4 seq/NO1_1.fastq
    for i in `tail -n+2 result/metadata.txt | cut -f 1`;do
      vsearch --fastq_mergepairs seq/${i}_1.fastq --reverse seq/${i}_2.fastq \
      --fastqout temp/${i}.merged.fq --relabel ${i}.
    done &
    for i in `tail -n+2 result/metadata2.txt | cut -f 1`;do
      usearch -fastx_relabel seq2/${i}.fastq -fastqout temp/${i}.merged.fq -prefix ${i}.
    done &
    cat temp/*.merged.fq > temp/all.fq
    ls -lsh temp/all.fq
    head -n 6 temp/all.fq|cut -c1-60
    vsearch --fastx_filter temp/all.fq \
      --fastq_stripleft 30 --fastq_stripright 20 \
      --fastq_maxee_rate 0.01 \
      --fastaout temp/filtered.fa
    head temp/filtered.fa