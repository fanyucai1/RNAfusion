# Email:yucai.fan@illumina.com
# 2021.12

import subprocess
import os
import argparse
import time
import re
from multiprocessing import Process

docker_name = "fanyucai1/rna:latest"

parse = argparse.ArgumentParser("This script will analysis RNA fusion.")
parse.add_argument("-p1", "--pe1", help="R1 fastq", required=True)
parse.add_argument("-p2", "--pe2", help="R2 fastq", required=True)
parse.add_argument("-r", "--ref", help="resource directory", required=True)
parse.add_argument("-p", "--prefix", help="prefix of output", required=True)
parse.add_argument("-d", "--dragen", help="hg38 dragen hash table", required=True)
parse.add_argument("-o", "--outdir", help="output directory", required=True)
args = parse.parse_args()

start = time.time()
###############################prepare input and output
args.pe1 = os.path.abspath(args.pe1)
args.pe2 = os.path.abspath(args.pe2)
args.ref = os.path.abspath(args.ref)
args.dragen = os.path.abspath(args.dragen)
if not os.path.exists(args.outdir):
    subprocess.check_call('mkdir -p %s' % (args.outdir), shell=True)
args.outdir = os.path.abspath(args.outdir)

if os.path.dirname(args.pe1) != os.path.dirname(args.pe2):
    print("%s and %s must be in the same directory." % (args.pe1, args.pe2))
    exit()
data_dir = os.path.dirname(args.pe1)

subprocess.check_call('mkdir -p %s/pizzy/output' % (args.outdir), shell=True)
subprocess.check_call('mkdir -p %s/star_fusion' % (args.outdir), shell=True)
subprocess.check_call('mkdir -p %s/Arriba' % (args.outdir), shell=True)
subprocess.check_call('mkdir -p %s/fusion_report' % (args.outdir), shell=True)
subprocess.check_call('mkdir -p %s/dragen' % (args.outdir), shell=True)


def shell_run(x):
    subprocess.check_call(x, shell=True)
def shell_run2(m,l,n):
    subprocess.check_call(m,shell=True)
    subprocess.check_call(l,shell=True)
    subprocess.check_call(n,shell=True)

##############################
docker_raw = "docker run -v %s:/reference/ -v %s:/project/ -v %s:/opt/ %s " % (args.ref, args.outdir, data_dir, docker_name)
################################star_fusion(https://github.com/STAR-Fusion/STAR-Fusion/wiki)
a = docker_raw + " sh /reference/STAR_fusion/STAR_fusion.sh /opt/%s /opt/%s /project/star_fusion" % (os.path.basename(args.pe1), os.path.basename(args.pe2))
################################Arriba(https://arriba.readthedocs.io/en/latest/)
b = docker_raw + " sh /reference/arriba/run.sh /opt/%s /opt/%s" % (os.path.basename(args.pe1), os.path.basename(args.pe2))
################################dragen
c = "dragen -f -r %s -1 %s -2 %s -a %s/gtf/gencode.v38.chr_patch_hapl_scaff.annotation.gtf.gz --output-dir %s/dragen/ --output-file-prefix %s " \
          " --RGID dragen_RGID --RGSM illumina --enable-rna true --enable-rna-gene-fusion true" % (args.dragen, args.pe1, args.pe2, args.ref, args.outdir, args.prefix)
################################pizzy(https://github.com/pmelsted/pizzly)
cmd1 = docker_raw + "/software/kallisto quant -i /reference/pizzly/index.idx --fusion -o /project/pizzy/output /opt/%s /opt/%s " % (
os.path.basename(args.pe1), os.path.basename(args.pe2))
cmd2 = docker_raw + "/software/pizzly-0.37.3/pizzly -k 31 --gtf /reference/pizzly/Homo_sapiens.GRCh38.104.chr_patch_hapl_scaff.gtf.gz --align-score 2 --insert-size 400 --fasta /reference/pizzly/Homo_sapiens.GRCh38.cdna.all.fa.gz --output /project/pizzy/output/%s /project/pizzy/output/fusion.txt" % (
    args.prefix)
cmd3 = docker_raw + "/software/python3/Python-v3.7.0/bin/python3 /software/pizzly-0.37.3/scripts/flatten_json.py /project/pizzy/output/%s.json /project/pizzy/output/%s_genetable.txt" % (
args.prefix, args.prefix)
##################################
p1 = Process(target=shell_run, args=(a,))
p2 = Process(target=shell_run, args=(b,))
p3 = Process(target=shell_run, args=(c,))
p4=Process(target=shell_run2,args=(cmd1,cmd2,cmd3,))
p1.start()
p2.start()
p3.start()
p4.start()
p1.join()
p2.join()
p3.join()
p4.join()

infile = open("%s/dragen/%s.fusion_candidates.final" % (args.outdir, args.prefix), "r")
outfile = open("%s/dragen/%s.fusion_candidates.final_new" % (args.outdir, args.prefix), "w")
for line in infile:
    line = line.strip()
    if not line.startswith("#") and re.search(r';', line.split("\t")[0]):
        array_1 = line.split("\t")[0].split('--')[0].split(";")
        array_2 = line.split("\t")[0].split('--')[1].split(";")
        a, b = [], []
        for i in range(0, len(array_1)):
            for j in range(0, len(array_2)):
                str_out = array_1[i] + "--" + array_2[j]
                for k in range(1, len(line.split("\t"))):
                    str_out += "\t" + line.split("\t")[k]
                outfile.write("%s\n" % (str_out))
    else:
        outfile.write("%s\n" % (line))
outfile.close()
################################fusion_report(https://github.com/matq007/fusion-report)
docker_final = docker_raw + " /software/python3/Python-v3.7.0/bin/fusion_report run %s /project/fusion_report/ " \
                            "/reference/fusion_report/db/ --arriba /project/Arriba/fusions.tsv --pizzly /project/pizzy/output/%s_genetable.txt " \
                            "--starfusion /project/star_fusion/star-fusion.fusion_predictions.tsv --dragen /project/dragen/%s.fusion_candidates.final_new " % (
               args.prefix, args.prefix, args.prefix)
subprocess.check_call(docker_final, shell=True)
subprocess.check_call(
    "mkdir -p %s/fusion_final_report && mv %s/fusion_report/index.html %s/fusion_final_report && rm -rf %s/fusion_report/" % (
    args.outdir, args.outdir, args.outdir, args.outdir), shell=True)

end = time.time()
print("Elapse time is %g seconds" % (end - start))
