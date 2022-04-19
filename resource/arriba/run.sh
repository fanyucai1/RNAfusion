export PATH=/software/STAR-2.7.6a/bin/Linux_x86_64:/software/samtools-v1.11/bin:$PATH && cd /project/Arriba && sh /reference/arriba/arriba_v2.1.0/run_arriba.sh /reference/arriba/arriba_v2.1.0/database/STAR_index_GRCh38_GENCODE38 /reference/arriba/arriba_v2.1.0/database/GENCODE38.gtf /reference/arriba/arriba_v2.1.0/database/GRCh38.fa /reference/arriba/arriba_v2.1.0/database/blacklist_hg38_GRCh38_v2.1.0.tsv.gz /reference/arriba/arriba_v2.1.0/database/known_fusions_hg38_GRCh38_v2.1.0.tsv.gz /reference/arriba/arriba_v2.1.0/database/protein_domains_hg38_GRCh38_v2.1.0.gff3 16 ${1} ${2}