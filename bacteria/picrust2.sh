mamba create -n picrust2 -c bioconda -c conda-forge picrust2=2.5.2
conda activate picrust2
picrust2_pipeline.py -s out.fasta -i bacteria2.txt -o picrust2_result -p 4