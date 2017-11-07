configfile: "bin/config.yaml"
import glob
SAMPLE_LIST      = config['sample']
SAMPLE_REP       = config['duplication']
# SAMPLE_LIST = ['B12-1','B12-2','B12-3']
genome            = config['datapath']['genome_path']
gtf               = config['datapath']['gtf_path']
trimmomatic       = config['trimmomatic']

Hisat2_splicesite = config['Hisat2_splicesite']
index             = config['datapath']['index_path']
Max_intron        = config['Max_intron']
Min_intron        = config['Min_intron']
Strand_specific   = config['Strand_specific']
cuff_lib_type     = config['cuff_lib_type']

hisat2_path       = config['outputpath']['hisat2']
stringtie_path    = config['outputpath']['stringtie']
gffcompare_path   = config['outputpath']['gffcompare']
ballgown_path     = config['outputpath']['ballgown']
cuffquant_path    = config['outputpath']['cuffquant']
cuffdiff_path     = config['outputpath']['cuffdiff']
lncRNA_path       = config['outputpath']['predict']
LABLES_string = []
for i in SAMPLE_LIST:
    tmp = []
    for x in SAMPLE_REP:
        a = "data/cuffquant/{0}-{1}/abundances.cxb".format(i,x)
        tmp.append(a)
    LABLES_string.append(','.join(tmp))
# print(LABLES_string)

rule all:
    input:#expand(ballgown_path  + "{sample}-{rep}/{sample}-{rep}.gtf",sample = SAMPLE_LIST,rep = SAMPLE_REP),
        # expand(cuffquant_path + "{sample}-{rep}/abundances.cxb",sample = SAMPLE_LIST,rep = SAMPLE_REP),
        # cuffdiff_path + "gene_exp.diff",
        # expand(gffcompare_path+ "{sample}-{rep}/{sample}-{rep}.filter.gtf",sample = SAMPLE_LIST,rep = SAMPLE_REP),
        # expand(gffcompare_path + "{sample}.filter.txt",sample = SAMPLE_LIST),
        # gffcompare_path + "final.txt",
        # cuffdiff_path + "gene_exp.diff",
        # gffcompare_path + "merged.annotated.filter.newname.fasta",
        # dynamic(lncRNA_path + 'split/split-{num}/split.fasta'),
        # lncRNA_path+ "final.newname.fasta"
        # dynamic(lncRNA_path + 'split/split-{num}/CPC/CPC.LncRNA.txt'),
        # dynamic(lncRNA_path + 'split/split-{num}/CNCI/CNCI.index'),
        # dynamic(lncRNA_path + 'split/split-{num}/pfam/Pfam_result.txt'),
        lncRNA_path + 'candiate.info'
    

rule hisat2:
    input:
        'raw_data/{sample}-{rep}_1.paired.fq', 'raw_data/{sample}-{rep}_2.paired.fq',
    output:
        hisat2_path + '{sample}-{rep}.sam'
    log:
        hisat2_path + '{sample}-{rep}.hisat2.log'
    shell:
        'hisat2 -x {index} -p 8 --min-intronlen {Min_intron} \
        --max-intronlen {Max_intron} --dta --rna-strandness {Strand_specific} \
        --time -1 {input[0]} -2 {input[1]} -S {output} 2> {log}'
rule sam2sortbam:
    input:
        hisat2_path + '{sample}-{rep}.sam'
    output:
        hisat2_path + '{sample}-{rep}.sort.bam'
    shell:
        "samtools sort -@ 4 -o {output} {input}"
        
rule stringtie:
    input:
        hisat2_path + '{sample}-{rep}.sort.bam'
    output:
        stringtie_path + '{sample}-{rep}.gtf'
    shell:
        "stringtie -p 8 -G {gtf} -o {output} {input}"

# rule filter:
#     input:
#         stringtie_path + '{sample}.gtf'
#     output:
#         stringtie_path + '{sample}.filter.gtf'
#     shell:
#         'python bin/filter_gtf.py -g {input[0]} -o {output[0]}'

rule mergelist:
    input:
        expand(stringtie_path + '{sample}-{rep}.gtf', sample = SAMPLE_LIST,rep = SAMPLE_REP)
    output:
        stringtie_path + "mergelist.txt"
    run:
         OUT = open(output[0],'w')
         OUT.write('\n'.join(input))
         OUT.close()

rule strintie_merged:
    input: 
        stringtie_path + "mergelist.txt"
    output:
        stringtie_path + "merge.gtf",
        stringtie_path + "merged.merge.gtf.tmap"
    shell:
        "stringtie --merge -p 8 -G {gtf} -o {output} {input}"

rule stringtie_for_Ballgown:
    input: 
        stringtie_path + "merge.gtf",hisat2_path + "{sample}-{rep}.sort.bam"
    output:
        ballgown_path  + "{sample}-{rep}/{sample}-{rep}.gtf"
    shell:
        "stringtie -e -B -p 12 -G {input[0]} -o {output} {input[1]}"

rule cuffquant:
    input:
        stringtie_path + "merge.gtf", hisat2_path + "{sample}-{rep}.sort.bam"
    params:
        cuffquant_path + "{sample}-{rep}"
    output:
        cuffquant_path + "{sample}-{rep}/abundances.cxb"
    shell:
        "cuffquant -o {params} -p 8 -b {genome} {cuff_lib_type} -u {input[0]} {input[1]} "

rule cuffdif:
    input:
        expand(cuffquant_path + "{sample}-{rep}/abundances.cxb",sample = SAMPLE_LIST,rep = SAMPLE_REP)
    output:
        cuffdiff_path + "gene_exp.diff"
    params:
        cuffdiff_path , stringtie_path + "merge.gtf", ','.join(SAMPLE_LIST), expand("{lables}",lables=LABLES_string)
    shell:
        "cuffdiff -L {params[2]} -o {params[0]} -p 8 -b {genome} -u {params[1]} \
        {cuff_lib_type} {params[3]}"

rule gffcompare:
    input:
        ballgown_path  + "{sample}-{rep}/{sample}-{rep}.gtf"
    output:
        gffcompare_path+ "{sample}-{rep}/{sample}-{rep}.annotated.gtf"
    params:
        gffcompare_path  + "{sample}-{rep}/{sample}-{rep}"
    shell:
        "gffcompare -G -r {gtf} -o {params} {input}"

rule filter_after_gffcompare:
    input:
        gffcompare_path+ "{sample}-{rep}/{sample}-{rep}.annotated.gtf"
    output:
        gffcompare_path+ "{sample}-{rep}/{sample}-{rep}.filter.gtf", gffcompare_path+ "{sample}-{rep}/{sample}-{rep}.filter.tracking"
    params:
        gffcompare_path+ "{sample}-{rep}/{sample}-{rep}",gffcompare_path+ "{sample}-{rep}/{sample}-{rep}.filter"
    shell:
        "python bin/filter.py -i {params[0]} -o {params[1]}"

rule deal_rep:
    input:
        expand(gffcompare_path+ "{{sample}}-{rep}/{{sample}}-{rep}.filter.tracking",rep = SAMPLE_REP)
    output:
        gffcompare_path + "{sample}.filter.txt", gffcompare_path + "{sample}.filter.3.txt"
    shell:
        "cat {input} | sed '/tracking-iso/d'|cut -f 4|sort| uniq -c | sed 's/^\s*//g' > {output[0]}; awk '{{if($1==3)print$2}}' {output[0]} > {output[1]}"

rule merge_all:
    input:
        expand(gffcompare_path + "{sample}.filter.3.txt",sample = SAMPLE_LIST)
    output:
        gffcompare_path + "final.txt"
    shell:
        "cat {input} | sort|uniq > {output}"

rule target_gtf:
    input:
        gffcompare_path + "final.txt", stringtie_path + "merge.gtf"
    output:
        lncRNA_path + "final.gtf"
    shell:
        "python bin/extract_target_gtf.py -t {input[0]} -g {input[1]} -o {output}"

#The gtf_to_fasta depend on the libz, export LD_LIBRARY_PATH=/home/jiajinbu/nas3/zhouran/zlib/lib:$LD_LIBRARY_PATH        
rule gtf2fasta:
    input:
        lncRNA_path + "final.gtf"
    output:
        lncRNA_path+ "final.fasta"
    shell:
        "gtf_to_fasta {input} {genome} {output}"

rule change_fa_name:
    input:
        lncRNA_path+ "final.fasta", stringtie_path + "merged.merge.gtf.tmap"
    output:
        lncRNA_path+ "final.newname.fasta", lncRNA_path + 'final.newname.list', lncRNA_path + 'final.info'
    run:
        import pandas as pd
        with open(input[0]) as IN, open(output[0],'w') as OUT, open(output[1],'w') as OUT1:
            OUT1.write('qry_id\tChr\n')
            for i in IN.readlines():
                if i.startswith('>'):
                    ALL = i.split(' ')
                    OUT1.write('\t'.join(ALL[1:3])+'\n')
                    OUT.write(">"+'|'.join(ALL[1:3]) + '\n')
                else:
                    OUT.write(i)
        IN1 = pd.read_csv(input[1],sep = '\t')
        IN2 = pd.read_csv(output[1],sep = '\t')
        RE  = pd.merge(IN2,IN1,how = 'inner',on = "qry_id")
        RE.to_csv(output[2],index = False,sep ='\t')



rule split_fa:
    input:
        lncRNA_path+ "final.newname.fasta"
    output:
        dynamic(lncRNA_path + 'split/split-{num}/split.fasta')
    params:
        1000, lncRNA_path + "split"
    shell:
        "python bin/split.py {input} {params[0]} {params[1]}"


rule CNCI:
    input:
        lncRNA_path + 'split/split-{num}/split.fasta'#,dynamic(lncRNA_path + 'split/split{clusterid}/split.fasta')
    output:
        lncRNA_path + 'split/split-{num}/CNCI/CNCI.index'
        # minifa + '/CNCI/CNCI.index'
    params:
        lncRNA_path + 'split/split-{num}/CNCI'
        # minifa + '/CNCI'
    shell:
        "python /home/jiajinbu/nas3/zhouran/soft/CNCI/CNCI.py -f {input} -o {params} -p 4 -m pl"

#export PATH=$PATH:/home/jiajinbu/nas3/zhouran/soft/hmmer-2.3.2/squid
rule CPC:
    input:
        lncRNA_path + 'split/split-{num}/split.fasta'
    output:
        lncRNA_path + 'split/split-{num}/CPC/CPC.LncRNA.txt', lncRNA_path + 'split/split-{num}/CPC/CPC.evidence.orf'
    params:
        lncRNA_path + 'split/split-{num}/CPC', lncRNA_path + 'split/split-{num}/CPC/CPC.evidence'
    shell:
        "/home/jiajinbu/nas3/zhouran/soft/cpc-0.9-r2/bin/run_predict_local.sh {input} {output[0]} {params[0]} {params[1]}"

rule pfam:
    input:
        lncRNA_path + 'split/split-{num}/split.fasta'
    output:
        lncRNA_path + 'split/split-{num}/pfam/Pfam_result.txt'
    shell:
        "pfam_scan.pl -translate orf -fasta {input} -dir /home/jiajinbu/nas3/zhouran/data/pfam/2017-10-9 -outfile {output} -cpu 4"

rule merge_all_prediction:
    input:
        lncRNA_path+ "final.newname.fasta",
        lncRNA_path + "final.gtf", 
        dynamic(expand(lncRNA_path + 'split/split-{num}/pfam/Pfam_result.txt',num = "{num}")),
        dynamic(expand(lncRNA_path + 'split/split-{num}/CPC/CPC.LncRNA.txt',num = "{num}")),
        dynamic(expand(lncRNA_path + 'split/split-{num}/CNCI/CNCI.index',num = "{num}"))
    output:
        lncRNA_path + 'candiate.gtf',lncRNA_path + 'candiate.list'
    params:
        lncRNA_path + 'split', lncRNA_path
    shell:
        "python bin/collapse.py -i {params[0]} -o {params[1]} -f {input[0]} -g {input[1]}"

rule get_candiate_info:
    input:
        stringtie_path + "merged.merge.gtf.tmap",
        lncRNA_path + 'candiate.list'
    output:
        lncRNA_path + 'candiate.info'
    run:
        import pandas as pd
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        IN1    = pd.read_csv(input[0],sep = '\t')
        IN2    = pd.read_csv(input[1],sep = '\t')
        out    = pd.merge(IN2,IN1,how = 'inner',on = "qry_id")
        count_ = pd.Series.value_counts(out['class_code'])
        
        plt.figure(figsize=(4,4))
        p = count_.plot(kind = 'bar',color = '#054E9F')
        plt.xticks(rotation=25)
        plt.title('LncRNA candiate type')
        plt.savefig(output[0]+'.pdf')

        out.to_csv(output[0],index = False,sep ='\t')

rule split_info:
    input:
        lncRNA_path + 'candiate.info', lncRNA_path + 'candiate.gtf'
    output:
        lncRNA_path + 'LincRNA.gtf'
    params:
        lncRNA_path
    shell:
        "python bin/split-info.py {input[0]} {input[1]} {params}"
