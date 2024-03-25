    db=~/amplicon2/public
    wd=~/amplicon2
    export PATH=$PATH:${db}/linux
    cd ${wd}
    head -n2 result/metadata.txt
    cat -A result/metadata.txt | head -n3
    head temp/filtered.fa
    vsearch --derep_fulllength temp/filtered.fa \
      --output temp/uniques.fa --relabel Uni --minuniquesize 10 --sizeout
    ls -lsh temp/uniques.fa
    head -n 2 temp/uniques.fa
    nohup vsearch --usearch_global temp/filtered.fa --db ${db}/usearch/reference.fasta --otutabout result/usearch/otutab.txt --id 0.97 --threads 38 &
    sed -i 's/\r//' result/raw/otutab.txt
    head -n3 result/raw/otutab.txt |cat -A

