import os
sample_info = open("/home/lixiang/project/sanqi-microbiome/data/sample.info.txt", "r")
sample_open = sample_info.readlines()

# 提取需要处理的数据
run_extracted = []
for line in sample_open:
    line2 = line.split(",")[6]
    if line2 == "plant":
        run_extracted.append(line.split(",")[0])

# 输出命令
sh_run = open("/home/lixiang/project/sanqi-microbiome/run.sh", "w")

for run in run_extracted:
    if run == "SRR14784230":
        # 将SRA序列拆分成两个fasta
        comm = "pfastq-dump --split-3 -t 30 -s data/download/%s -O data/fastq/" % run
        sh_run.write(comm + "\n")
        comm = "echo =================================================\necho 已完成%s的数据拆分！\necho =================================================" % run
        sh_run.write(comm + "\n")
        # fastqc质控
        comm = "fastqc data/fastq/* -t 2 -o results/0-fastqc/"
        sh_run.write(comm + "\n")
        comm = "echo =================================================\necho 已完成%s的数据质控！\necho =================================================" % run
        sh_run.write(comm + "\n")
        # 序列比对到宿主基因组上
        comm = "/usr/bin/bowtie2 --very-sensitive-local -p 15 -x data/host-genome/sanqi -1 data/fastq/%s_1.fastq -2 data/fastq/%s_2.fastq -S results/1-sambybowtie/%s.sam" % (run, run, run)
        sh_run.write(comm + "\n")
        comm = "echo =================================================\necho 已将%s比对到宿主基因组！\necho =================================================" % run
        sh_run.write(comm + "\n")
        # sam文件转换成bam文件
        comm = "samtools view -@ 30 -bS results/1-sambybowtie/%s.sam > results/2-bamfile/%s.bam" % (run, run)
        sh_run.write(comm + "\n")
        comm = "echo =================================================\necho 已完成%s的sam转bam！\necho =================================================" % run
        sh_run.write(comm + "\n")
        # 提取没有比对上的bam文件
        comm = "samtools view -@ 30 -b -f 12 -F 256 results/2-bamfile/%s.bam > results/3-unmapped.bam/%s.bam" % (run, run)
        sh_run.write(comm + "\n")
        comm = "echo =================================================\necho 已提取%s未比对上宿主的bam文件！\necho =================================================" % run
        sh_run.write(comm + "\n")
        # bam文件排序
        comm = "samtools sort -n results/3-unmapped.bam/%s.bam -o results/4-sorted.bam/%s.sorted.bam" % (run, run)
        sh_run.write(comm + "\n")
        comm = "echo =================================================\necho 已排序%s的bam文件！\necho =================================================" % run
        sh_run.write(comm + "\n")
        # 删除sam文件
        comm = "rm results/1-sambybowtie/%s.sam" % run
        sh_run.write(comm + "\n")
        comm = "echo =================================================\necho 已删除%s的sam文件！\necho =================================================" % run
        sh_run.write(comm + "\n")
        # 删除第一次比对的bam文件
        comm = "rm results/2-bamfile/%s.bam" % run
        sh_run.write(comm + "\n")
        comm = "echo =================================================\necho 已删除%s的第一次bam文件！\necho =================================================" % run
        sh_run.write(comm + "\n")
        # 删除未排序的bam文件
        comm = "rm results/3-unmapped.bam/%s.bam" % run
        sh_run.write(comm + "\n")
        comm = "echo =================================================\necho 已删除%s未排序的bam文件！\necho =================================================" % run
        sh_run.write(comm + "\n")
        # bam文件转换成fastq文件
        comm = "samtools bam2fq results/4-sorted.bam/%s.sorted.bam > results/5-new.fastq/%s.fastq" % (run, run)
        sh_run.write(comm + "\n")
        comm = "echo =================================================\necho 已完成%s的bam文件转fastq文件！\necho =================================================" % run
        sh_run.write(comm + "\n")
        # 第二次质控
        comm = "fastqc results/5-new.fastq/%s.fastq -o results/6-new.fastqc/" % run
        sh_run.write(comm + "\n")
        comm = "echo =================================================\necho 已完成%s的第二次质控！\necho =================================================" % run
        sh_run.write(comm + "\n")
        # megahit组装
        comm = "mkdir results/8-megahit.results/%s" % run
        sh_run.write(comm + "\n")
        comm = "megahit --12 results/5-new.fastq/%s.fastq --k-list 35 -o results/8-megahit.results/%s/out" % (run, run)
        sh_run.write(comm + "\n")
        comm = "echo =================================================\necho 已完成%s的megahit组装！\necho =================================================" % run
        sh_run.write(comm + "\n")
        # metaquast评估组装质量
        comm = "mkdir results/results/9-metaquast.results/%s" % run
        sh_run.write(comm + "\n")
        comm = "/usr/bin/metaquast.py -o results/results/9-metaquast.results/%s results/8-megahit.results/%s/out/final.contigs.fa" % (run, run)
        sh_run.write(comm + "\n")
        comm = "echo =================================================\necho 已完成%s的组装质量评估！\necho =================================================" % run
        sh_run.write(comm + "\n")
sh_run.close()