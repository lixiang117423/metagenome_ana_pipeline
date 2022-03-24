# metagenome_ana_pipeline

## 参考流程

　　[博客地址](https://rpubs.com/ednachiang/MetaG_Pipeline)

## 软件安装

### Bowtie2

* [下载地址](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#obtaining-bowtie-2)
* [文献链接](https://www.nature.com/articles/nmeth.1923)
* 下载完成后解压到指定目录即可
* 输入下方代码检验是否安装成功。

```shell
bowtie2 -h
```

### metaphlan2

　　[下载地址](https://huttenhower.sph.harvard.edu/metaphlan2/)（压缩文件500M+，用的谷歌云盘，需要挂梯子）

　　​[官网教程](https://github.com/biobakery/MetaPhlAn2)

　　[文献链接](https://www.nature.com/articles/nmeth.3589?report=reader)

　　直接下载解压即可。

　　报错`Error: Unable to find the mpa_pkl file at: /usr/bin/db_v20/mpa_v20_m200.pkl`，解决方法是把对应的数据库放到指定的位置，软连接即可：

```r
sudo ln -s /opt/software/metaphlan2/db_v20/ /usr/bin/
```

### MegaHit

　　[软件地址](https://github.com/voutcn/megahit)

　　[文献链接](https://www.sciencedirect.com/science/article/abs/pii/S1046202315301183)

　　直接用conda安装：

```r
conda install -c bioconda megahit
```

　　检测是否安装成功：

```r
megahit -h
```

### quast

　　[官网地址](http://quast.sourceforge.net/metaquast)

　　[文献链接](https://academic.oup.com/bioinformatics/article/32/7/1088/1743987?login=true)

　　安装流程：

```r
wget https://downloads.sourceforge.net/project/quast/quast-5.0.2.tar.gz
tar -xzf quast-5.0.2.tar.gz
cd quast-5.0.2
./setup.py install

sudo ln -s /opt/software/quast/quast.py /usr/bin/
sudo ln -s /opt/software/quast/metaquast.py /usr/bin/
```

### 拆分fastq

　　[项目地址](https://gist.github.com/nathanhaigh/3521724)

```r
#!/bin/bash
# Usage: deinterleave_fastq.sh < interleaved.fastq f.fastq r.fastq [compress]
#
# Deinterleaves a FASTQ file of paired reads into two FASTQ
# files specified on the command line. Optionally GZip compresses the output
# FASTQ files using pigz if the 3rd command line argument is the word "compress"
#
# Can deinterleave 100 million paired reads (200 million total
# reads; a 43Gbyte file), in memory (/dev/shm), in 4m15s (255s)
#
# Latest code: https://gist.github.com/3521724
# Also see my interleaving script: https://gist.github.com/4544979
#
# Inspired by Torsten Seemann's blog post:
# http://thegenomefactory.blogspot.com.au/2012/05/cool-use-of-unix-paste-with-ngs.html

# Set up some defaults
GZIP_OUTPUT=0
PIGZ_COMPRESSION_THREADS=10

# If the third argument is the word "compress" then we'll compress the output using pigz
if [[ $3 == "compress" ]]; then
  GZIP_OUTPUT=1
fi

if [[ ${GZIP_OUTPUT} == 0 ]]; then
  paste - - - - - - - -  | tee >(cut -f 1-4 | tr "\t" "\n" > $1) | cut -f 5-8 | tr "\t" "\n" > $2
else
  paste - - - - - - - -  | tee >(cut -f 1-4 | tr "\t" "\n" | pigz --best --processes ${PIGZ_COMPRESSION_THREADS} > $1) | cut -f 5-8 | tr "\t" "\n" | pigz --best --processes ${PIGZ_COMPRESSION_THREADS} > $2
fi
```

　　修改为可执行文件：

```r
sudo chmod +x /usr/bin/deinterleave_fastq.sh
```

### Prodigal

　　[项目地址](https://github.com/hyattpd/Prodigal)

　　[文献链接](https://link.springer.com/article/10.1186/1471-2105-11-119)

### MMseqs2

　　[官网地址](https://github.com/soedinglab/MMseqs2)

　　安装：

```shell
conda install -c conda-forge -c bioconda mmseqs2
```

## 创建文件夹

```r
mkdir data/host-genome
mkdir results
mkdir results/0-fastqc
mkdir results/1-sambybowtie
mkdir results/2-bamfile
mkdir results/3-unmapped.bam
mkdir results/4-sorted.bam
mkdir results/5-new.fastq
mkdir results/6-new.fastqc
mkdir results/7-metaphlan
mkdir results/8-megahit
mkdir results/9-metaquast
mkdir results/10-prodigal
mkdir results/11-mmseqs2
mkdir results/12-diamond
```

## 数据预处理

　　如果是从NCBI　SRA数据库下载的数据，那么需要先对数据进行转换和质控。

* 将SRA数据转换成fastq数据

  ```shell
  pfastq-dump --split-3 -t 30 -s data/download/SRR14784195 -O data/fastq
  ```

* 数据质控

  ```shell
  fastqc data/fastq/* -o results/0-fastqc
  cd results/0-fastqc
  multiqc .
  ```

## 去除宿主基因组

* 构建宿主基因组索引

  ```shell
  cd data/host-genome
  /usr/bin/bowtie2-build --seed 707 --threads 30 /home/publicdata/refgenome/panaxnotoginseng/yangshengchao/genome/Panax_notoginseng.genome.2.fasta data/host-genome/sanqi
  ```

* 将数据mapping到宿主基因组上

  ```shell
  /usr/bin/bowtie2 --very-sensitive-local -p 15 -x data/host-genome/sanqi -1 data/fastq/SRR14784195_1.fastq -2 data/fastq/SRR14784195_2.fastq -S results/1-sambybowtie/SRR14784159.sam
  ```

* 将sam文件转换成bam文件

  ```shell
  samtools view -@ 30 -bS results/1-sambybowtie/SRR14784159.sam > results/2-bamfile/SRR14784159.bam
  ```

* 提取没有mapping上宿主基因组的序列

  ```shell
  samtools view -@ 30 -b -f 12 -F 256 results/2-bamfile/SRR14784159.bam > results/3-unmapped.bam/SRR14784159.bam
  ```

* bam文件排序

  ```shell
  samtools sort -n results/3-unmapped.bam/SRR14784159.bam -o results/4-sorted.bam/SRR14784159.sorted.bam
  ```

* bam文件转换成fastq文件

  ```shell
  samtools bam2fq results/4-sorted.bam/SRR14784159.sorted.bam > results/5-new.fastq/SRR14784159.fastq
  ```

* 重新质控新的结果

  ```shell
  fastqc results/5-new.fastq/SRR14784159.fastq -o results/6-new.fastqc/
  ```

* 质控完成后将所有的质控结果整合成一个

  ```r
  cd 
  ```

## MetaPhlan2分类

```r
/usr/bin/metaphlan2.py ../10-idba.results/fq2fa/SRR14784159.fa --bowtie2out metagenome.bowtie2.2.bz2 --nproc 10 --input_type fasta > ./SRR14784159.fa.txt
```

## MegaHit组装

```r
megahit --12 results/5-new.fastq/SRR14784159.fastq --k-list 35 -o results/11-megahit.results/out
```

　　组装完后成生成的文件：

　　![image-20220323194850-5f3zb5y](https://xiang-1257290193.cos.ap-guangzhou.myqcloud.com/Typora/image-20220323194850-5f3zb5y.png)

## Metaquast评估组装质量

```r
/usr/bin/quast.py -o results/13-metaquast.results results/11-megahit.results/out/final.contigs.fa
```

　　![image-20220323195214-hmdug9k](https://xiang-1257290193.cos.ap-guangzhou.myqcloud.com/Typora/image-20220323195214-hmdug9k.png)

## 分箱binning

```r
mkdir -p results/binning
mkdir -p results/binning/deinterleave_fastq
mkdir -p results/binning/index
mkdir -p results/binning/samfile
```

* 先将先前合并的fastq文件拆成两个。输入文件为用于组装的fastq文件。

```r
/usr/bin/deinterleave_fastq.sh < results/5-new.fastq/SRR14784159.fastq results/binning/deinterleave_fastq/R1.fastq results/binning/deinterleave_fastq/R2.fastq
```

* 构建组装索引。输入文件是组装好的fasta文件。

```r
bowtie2-build --seed 707 results/11-megahit.results/out/final.contigs.fa results/binning/index/index
```

* 将拆分的两个fastq文件比对到构建的基因组上去。

```r
bowtie2 --sensitive-local -p 10 --seed 707 -x results/binning/index/index -1 results/binning/deinterleave_fastq/R1.fastq -2 results/binning/deinterleave_fastq/R2.fastq -S results/binning/samfile/sample.sam
```

## Prodigal预测基因

　　输入MegaHit组装的结果，直接预测基因。输出蛋白序列、CDS序列及gff文件等。

```r
prodigal -i results/11-megahit.results/out/final.contigs.fa -a results/prodigal.results/res.pep -d results/prodigal.results/res.cds -f gff -g 11 -o results/prodigal.results/res.gff -p single -s results/prodigal.results/res.stat
```

　　筛选特定长度的基因：

```python
import os
import re
import argparse

parser = argparse.ArgumentParser(description="You should add those parameters")
parser.add_argument("-i","--input", type = str, help = "The input fasta file")
parser.add_argument("-min","--min_length", type = int, help = "The minimal length of sequences")
parser.add_argument("-max","--max_length", type = int, default = 0, help = "The maximal length of sequences")
parser.add_argument("-o","--output", type = str, help = "The output fasta file")
#args = parser.parse_args()
args =parser.parse_known_args()[0]

# 输入的fasta文件转换成字典
file_open = open(args.input,"r")
file_read = file_open.readlines()
res_dict = {}
for line in file_read:
    if re.match(">",line):
        #print(line)
        res_dict[line] = ""
        flag = line
    else:
        res_dict[flag] = res_dict[flag] + line

# 根据长度进行筛选
# 根据长度进行筛选
final_results = []
for key, value in res_dict.items():
    if args.max_length == 0:
        if len(value) >= args.min_length:
            final_results.append(key)
            final_results.append(value)
    else:
        if args.min_length <= len(value) <= args.max_length:
            final_results.append(key)
            final_results.append(value)

# 保存文件
file_new = open(args.output,"w")
for i in final_results:
    file_new.write(i + "\n")

file_new.close()
```

```shell
python /usr/bin/get_gene_from_prodigal_with_length.py -i results/prodigal.results/res.cds -min 100 -o results/prodigal.results/res.cds.new.fa
```

　　删除生成的`gff`文件中的注释行：

```shell
sed '/#/d' res.gff > res.new.gff
```

## MMseq2构建基因集

```shell
mkdir results/mmseq2.results
```

```shell
mmseqs easy-cluster results/prodigal.results/res.cds.new.fa results/mmseq2.results/clusterRes results/mmseq2.results/tmp --min-seq-id 0.5 -c 0.8 --cov-mode 1
```

* `results/mmseq2.results/clusterRes`是输出文件
* `results/mmseq2.results/tmp`是临时文件

　　![image-20220324112938-f8c21nt](https://xiang-1257290193.cos.ap-guangzhou.myqcloud.com/Typora/image-20220324112938-f8c21nt.png)

## 蛋白序列比对到NCBI NR库

　　[NCBI NR库下载地址](https://nmdc.cn/datadownload)

