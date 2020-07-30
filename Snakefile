import pandas as pd
import sys,os
##隐式通配符：id,scaffold
configfile: "config.yaml"
##一些前处理工作
def get_filePath_fileName_fileExt(fileUrl):
    """
    获取文件路径， 文件名， 后缀名
    :param fileUrl:
    :return:
    """
    filepath, tmpfilename = os.path.split(fileUrl)
    shotname, extension = os.path.splitext(tmpfilename)
    return filepath, shotname, extension

smk_script_dir = sys.path[0] #获取snakefile所在目录

##组合新的文件名
ref_dict = os.path.join(get_filePath_fileName_fileExt(config["reference_genome"])[0],
             get_filePath_fileName_fileExt(config["reference_genome"])[1] + ".dict")
sys.stderr.write(f"[INFO] ref_dict:{ref_dict}\n")
#判断工作路径下用于软件输出临时文件的tmp目录是否存在，否则创建
tmp_dir = "tmp_dir"
if not os.path.exists(tmp_dir):
      os.makedirs(tmp_dir)

#读取包含samples id及样本对应路径的csv文件,包含两列ID，file
samples_pd = pd.read_csv(config["samples_csv_file"])

#去除所有单元格字符串两边的空格
samples_pd['ID'] = samples_pd['ID'].map(str.strip)
samples_pd['file'] = samples_pd['file'].map(str.strip)
samples_pd.dropna(inplace=True)#删除不完整的行
samples_pd.drop_duplicates(subset=['ID'], keep='first', inplace=True)#去重
samples_pd = samples_pd.set_index("ID", drop=False) #设置行索引
Individuals = samples_pd.index
sys.stderr.write(f"[INFO] Individuals:{Individuals}\n")

#读取包含scaffolds列表的文件
scaffolds_list_file_fd = open(config["scaffolds_list_file"])
scaffolds = []
for aline in scaffolds_list_file_fd:
    scaffolds.append(aline.split()[0].strip()) #行中有空格时只取第一列，并去除所有字符串两边的空格或其它空字符
scaffolds = list(set(scaffolds)) #列表元素去重
if '' in scaffolds:
    scaffolds.remove("")#去除列表中的空字符串项
sys.stderr.write(f"[INFO] scaffolds:{scaffolds}\n")

rule all:
    input:
        quals = "results/07_plots/quals.svg",
        stat = expand("results/02_mapped_bam/{a}/{a}_stat.txt",a=Individuals),
        stats = expand("results/02_mapped_bam/{a}/{a}_stats.txt",a=Individuals)

rule genome_index:
    input:
        config["reference_genome"]
    output:
        amb = config["reference_genome"] + ".amb",
        ann = config["reference_genome"] + ".ann",
        fai = config["reference_genome"] + ".fai",
        dict = ref_dict
    shell:
        "bwa index {input} && "
        "samtools faidx {input} && "
        "picard CreateSequenceDictionary R= {input} O= {output.dict}"

rule map_and_sort:
    input:
        ref = config["reference_genome"],
        index = config["reference_genome"] + ".amb",
        fq = lambda wildcards: samples_pd.loc[wildcards.id,'file']
    output:
        bam = "results/02_mapped_bam/{id}/{id}_sorted.bam",
        bai = "results/02_mapped_bam/{id}/{id}_sorted.bai"
    threads: 2
    params:
        rg = r"@RG\tID:{id}\tPU:{id}\tLB:{id}\tPL:ILLUMINA\tSM:{id}" #前面加r代表不转义
    #log:
     #   "logs/bwa_mem/{id}.log"
    benchmark:
        "benchmarks/{id}.bwa.benchmark.txt"
    shell:
        #在mapping时顺便增加RG头信息
        #bwa mem 加入-M参数以兼容Picard
        #-v 1:只输出err信息；-v 2: 输出WARNING和err
        #"bwa mem -v 1 -t {threads} -R " + r"'@RG\tID:{wildcards.id}\tPU:{wildcards.id}'"
        "bwa mem -v 1 -t {threads} -R '{params.rg}' -M {input.ref} {input.fq} | "
            "picard -Xmx4g -XX:ParallelGCThreads={threads} -Djava.io.tmpdir={tmp_dir} SortSam "
                "I=/dev/stdin "
                "O={output.bam} "
                "SORT_ORDER=coordinate "
                "CREATE_INDEX=true "
                "VERBOSITY=WARNING" #注意{params.rg}两边的引号;注意每一行结尾的空格

rule mark_duplication:
    input:
        bam = "results/02_mapped_bam/{id}/{id}_sorted.bam",
        bai = "results/02_mapped_bam/{id}/{id}_sorted.bai"
    output:
        bam = "results/02_mapped_bam/{id}/{id}_markdup.bam",
        bai = "results/02_mapped_bam/{id}/{id}_markdup.bai",
        metrics_file = "results/02_mapped_bam/{id}/{id}_metrics.txt"
    threads: 2
    #conda:
    #    "envs/samtools.yaml"
    shell:
        "picard -Xmx4g -XX:ParallelGCThreads={threads} -Djava.io.tmpdir={tmp_dir} MarkDuplicates \
            INPUT={input.bam} \
            OUTPUT={output.bam} \
            METRICS_FILE={output.metrics_file} \
            CREATE_INDEX=true \
            VERBOSITY=WARNING"

rule bam_mapping_stat:
    input:
        bam = "results/02_mapped_bam/{id}/{id}_markdup.bam",
        bai = "results/02_mapped_bam/{id}/{id}_markdup.bai"
    output:
        stat = "results/02_mapped_bam/{id}/{id}_stat.txt",
        stats = "results/02_mapped_bam/{id}/{id}_stats.txt"
    shell:
        "samtools flagstat {input.bam} > {output.stat}; "
        "samtools stats {input.bam} > {output.stats}" #注意第一行后面的分号;

rule HaplotypeCaller:
    input:
        ref = config["reference_genome"],
        bam = "results/02_mapped_bam/{id}/{id}_markdup.bam",
        bai = "results/02_mapped_bam/{id}/{id}_markdup.bai",
        fai = config["reference_genome"] + ".fai",
        dict = ref_dict
    output:
        gvcf = "results/03_gvcfs_per_sample/{id}/{id}.g.vcf.gz",
        tbi = "results/03_gvcfs_per_sample/{id}/{id}.g.vcf.gz.tbi"
    threads:2
    shell:
        "gatk --java-options '-Xmx8g -XX:ParallelGCThreads={threads}' HaplotypeCaller \
        --verbosity ERROR \
        --TMP_DIR {tmp_dir} \
        -R {input.ref} \
        -I {input.bam} \
        -O {output.gvcf} \
        -ERC GVCF"

rule gvcf_sample_map:
    input:
        gvcf = expand("results/03_gvcfs_per_sample/{a}/{a}.g.vcf.gz",a=Individuals)
    output:
        map = "results/04_DB/samples_map.txt"
    params:
        map_list = expand("{a}\tresults/03_gvcfs_per_sample/{a}/{a}.g.vcf.gz",a=Individuals)
    run:
        samples_map_fd = open(output.map,'w')
        for a_map in params.map_list:
            samples_map_fd.write(a_map + '\n') #注意加回车换行
        samples_map_fd.close()

rule GenomicsDBImport:
    input:
        samples_map = "results/04_DB/samples_map.txt"
    output:
        db = "results/04_DB/{scaffold}.db" #伪需求文件
    shell:
        """
        if [ -d "04_DB/{wildcards.scaffold}" ]
        then
            rm -R "04_DB/{wildcards.scaffold}"
        fi
        gatk --java-options '-Xmx8g -XX:ParallelGCThreads=1' \
         GenomicsDBImport \
       --TMP_DIR {tmp_dir} \
       --genomicsdb-workspace-path results/04_DB/{wildcards.scaffold} \
       -L {wildcards.scaffold} \
       --sample-name-map {input.samples_map} && \
       touch {output.db}
       """

rule GenotypeGVCFs:
    input:
        db = "results/04_DB/{scaffold}.db",#伪需求文件
        ref = config["reference_genome"]
    output:
        gvcf = "results/05_VCFs/{scaffold}/{scaffold}.vcf.gz"
    threads:4
    shell:
         #GenotypeGVCFs #使用-new-qual计算方法将大大加速大样本量的计算速度，该选项在4.10.0版本后为默认选项
        "gatk --java-options '-Xmx16g -XX:ParallelGCThreads={threads}' GenotypeGVCFs \
        -R {input.ref} \
        --TMP_DIR {tmp_dir} \
        -new-qual \
        -V gendb://results/04_DB/{wildcards.scaffold} \
        -O {output.gvcf}"

rule concat_and_filter:
    input:
        gvcfs = expand("results/05_VCFs/{a}/{a}.vcf.gz",a=scaffolds)
    output:
        vcf = "results/06_VCF/all.clean.vcf.gz",
        index = "results/06_VCF/all.clean.vcf.gz.tbi"
    threads:2
    params:
        filter = "AVG(FMT/DP)>5&AVG(FMT/DP)<40&TYPE=\"snp\"&QUAL>20&N_ALT=1&AVG(FMT/GQ)>20&MAF>0.01&F_MISSING<0.5"
        #注意"snp"引号前的转义符
    shell:
        "bcftools concat \
        -O u\
        --threads {threads} \
        {input.gvcfs} \
    | bcftools filter -i \
        '{params.filter}' \
        -Oz --threads {threads} -o {output.vcf} - \
    && tabix -p vcf {output.vcf}" #注意{params.filter}两边的引号

rule plot_quals:
    input:
        vcf = "results/06_VCF/all.clean.vcf.gz"
    output:
        "results/07_plots/quals.svg"
    script:
        "scripts/plot-quals.py"
