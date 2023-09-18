# PROJECT = "/mnt/project/Metabolomics/20210717_QE_1st"
# SAMPLES = ["QC1", "QC2", "QC3", "DN1", "DP1", "DG1", "DN2", "DP2", "DG2", "QC4", "DN3", "DP3", "DG3", "DN4", "DP4", "DG4", "QC5", "QC6"]
# THREAD = "8"
# ProgramDir = "/home/project/pipelines/Metabolomics/Metabolomics_pipe_V1.0"

SAMPLE_INFO = "sample_info.xls"
VS_GROUPS = "vs_group.xls"
samples_list = ",".join(SAMPLES)
mode = "all" # or plant
## rule all 
rule all:
    input:
        expand("Input/pos/{sample}.mzXML", sample = SAMPLES),
        expand("Input/neg/{sample}.mzXML", sample = SAMPLES),
        "Output/FeatureDetection/pos/pos_result.txt",
        "Output/FeatureDetection/neg/neg_result.txt",
        "Output/metaX/pos/report.html",
        "Output/metaX/neg/report.html",
        "Output/2.MetaboliteIdentification/identification.stat.csv",
        "Output/3.MetaboliteQuantification/Feautre2quant/pos/sample.correlation.png",
        "Output/4.MetaboliteComparison/comparison.stat.csv",
        "Output/5.MultiGroup/multigroup.significant.stat.csv",
        "Output/3.MetaboliteQuantification/neg-norm-metaboAnalystInput.csv",
        "Output/3.MetaboliteQuantification/pos-norm-metaboAnalystInput.csv",
        "Output/report.html"

rule raw2mzXML_pos:
    input:
        "raw/pos/{sample}.raw"
    output:
        "Input/pos/{sample}.mzXML"
    shell:
        "docker run -i --rm -u project:project -e WINEDEBUG=-all -v /mnt/project:/mnt/project chambm/pwiz-skyline-i-agree-to-the-vendor-licenses.fixuid wine msconvert {PROJECT}/{input} --mzXML -o {PROJECT}/Input/pos"


rule raw2mzXML_neg:
    input:
        "raw/neg/{sample}.raw"
    output:
        "Input/neg/{sample}.mzXML"
    shell:
        "docker run -i --rm -u project:project -e WINEDEBUG=-all -v /mnt/project:/mnt/project chambm/pwiz-skyline-i-agree-to-the-vendor-licenses.fixuid wine msconvert {PROJECT}/{input} --mzXML -o {PROJECT}/Input/neg"


rule xcms_pos:
    input:
        expand("Input/pos/{sample}.mzXML", sample = SAMPLES)
    output:
        "Output/FeatureDetection/pos/pos_result.txt"
    shell:
        "Rscript {ProgramDir}/bin/xcms.R {samples_list} pos Output/FeatureDetection/pos {SAMPLE_INFO}"


rule xcms_neg:
    input:
        expand("Input/neg/{sample}.mzXML", sample = SAMPLES)
    output:
        "Output/FeatureDetection/neg/neg_result.txt"
    shell:
        "Rscript {ProgramDir}/bin/xcms.R {samples_list} neg Output/FeatureDetection/neg {SAMPLE_INFO}"

rule metaX_pos:
    input:
        "raw/CD_pos.xlsx"
    output:
        "Output/metaX/pos/report.html",
        "Output/metaX/pos/data/pos-raw-metaboAnalystInput.csv",
        "Output/metaX/pos/data/pos-norm-metaboAnalystInput.csv",
        "Output/metaX/pos/data/pos-quant.txt"
    shell:
        "Rscript {ProgramDir}/bin/metaX.R {input} pos Output/metaX/pos {SAMPLE_INFO} {VS_GROUPS} {THREAD}"

rule metaX_neg:
    input:
        "raw/CD_neg.xlsx"
    output:
        "Output/metaX/neg/data/neg-raw-metaboAnalystInput.csv",
        "Output/metaX/neg/data/neg-norm-metaboAnalystInput.csv",
        "Output/metaX/neg/report.html",
        "Output/metaX/neg/data/neg-quant.txt"
    shell:
        "Rscript {ProgramDir}/bin/metaX.R {input} neg Output/metaX/neg {SAMPLE_INFO} {VS_GROUPS} {THREAD}"
rule MetaboliteIdentification:
    input:
        peak_neg = "raw/CD_neg.xlsx",
        peak_pos = "raw/CD_pos.xlsx",
        intensity_neg_raw = "Output/metaX/neg/data/neg-raw-metaboAnalystInput.csv",
       	intensity_neg_norm= "Output/metaX/neg/data/neg-norm-metaboAnalystInput.csv",
        intensity_pos_raw = "Output/metaX/pos/data/pos-raw-metaboAnalystInput.csv",
        intensity_pos_norm = "Output/metaX/pos/data/pos-norm-metaboAnalystInput.csv",
    output:
        identification_intensity_raw = "Output/2.MetaboliteIdentification/identification.raw.intensity.csv",
        identification_intensity_norm = "Output/2.MetaboliteIdentification/identification.norm.intensity.csv",
        stats = "Output/2.MetaboliteIdentification/identification.stat.csv"
    params:
        db = f"{ProgramDir}/bin/database/metabolics.db",
        kegg = f"{ProgramDir}/bin/database/pathway_tree2.json",
        compounds = f"{ProgramDir}/bin/database/compounds.json",
        mode=mode,
    threads: 8
    shell:
        """
        python {ProgramDir}/bin/script/annotate_compound.py --db {params.db} --kegg {params.kegg} --peak-neg {input.peak_neg} --peak-pos {input.peak_pos} \
                                    --intensity-pos {input.intensity_pos_norm} --intensity-neg {input.intensity_neg_norm} \
                                    --output {output.identification_intensity_norm} --stats {output.stats} --sample-file Output/metaX/neg/sampleList.txt --mode {params.mode} &&

        python {ProgramDir}/bin/script/annotate_compound.py --db {params.db} --kegg {params.kegg} --peak-neg {input.peak_neg} --peak-pos {input.peak_pos} \
                                    --intensity-pos {input.intensity_pos_raw} --intensity-neg {input.intensity_neg_raw} \
                                    --output {output.identification_intensity_raw} --stats {output.stats} --sample-file Output/metaX/neg/sampleList.txt --mode {params.mode} &&

        python {ProgramDir}/bin/script/2.MetaboliteIdentification/HMDB.py --input {output.identification_intensity_raw} --outdir Output/2.MetaboliteIdentification/HMDB &&
        
        python {ProgramDir}/bin/script/kegg_analysis.py --input {output.identification_intensity_raw} --outdir Output/2.MetaboliteIdentification/KEGG --compound {params.compounds} --mode {params.mode} &&
        
        Rscript {ProgramDir}/bin/script/2.MetaboliteIdentification/KEGG.R -l Output/2.MetaboliteIdentification/KEGG/kegg.level2.csv -k Output/2.MetaboliteIdentification/KEGG/kegg.csv -o Output/2.MetaboliteIdentification/KEGG && 

        Rscript {ProgramDir}/bin/script/color_pathway.R --kegg_result Output/2.MetaboliteIdentification/KEGG/kegg.csv --kegg_dir {ProgramDir}/bin/database/KGML --out_dir Output/2.MetaboliteIdentification/KEGG/kegg_map --thread {threads}

        """

rule MetaboliteQuantification:
    input:
        intensity_norm_neg = "Output/metaX/neg/data/neg-norm-metaboAnalystInput.csv",
        intensity_norm_pos = "Output/metaX/pos/data/pos-norm-metaboAnalystInput.csv",
        identification_intensity = "Output/2.MetaboliteIdentification/identification.norm.intensity.csv"
    output:
        neg = "Output/3.MetaboliteQuantification/Feautre2quant/neg/sample.correlation.png",
        pos = "Output/3.MetaboliteQuantification/Feautre2quant/pos/sample.correlation.png",
    shell:
        """
        Rscript {ProgramDir}/bin/script/3.MetaboliteQuantification/Feature2quant.R --intensity {input.intensity_norm_neg} --pca_rds Output/metaX/neg/data/neg-pca.rds --sample_file Output/metaX/neg/sampleList.txt -o Output/3.MetaboliteQuantification/Feautre2quant/neg --annotation {input.identification_intensity} --type neg && 

        Rscript {ProgramDir}/bin/script/3.MetaboliteQuantification/Feature2quant.R --intensity {input.intensity_norm_pos} --pca_rds Output/metaX/pos/data/pos-pca.rds --sample_file Output/metaX/pos/sampleList.txt -o Output/3.MetaboliteQuantification/Feautre2quant/pos --annotation {input.identification_intensity} --type pos &&
        Rscript {ProgramDir}/bin/script/3.MetaboliteQuantification/metabolite2quant.R --intensity  {input.intensity_norm_neg} --annotation {input.identification_intensity} --type neg --output_dir Output/3.MetaboliteQuantification/Metabolite2quant/neg &&

        Rscript {ProgramDir}/bin/script/3.MetaboliteQuantification/metabolite2quant.R --intensity  {input.intensity_norm_neg} --annotation {input.identification_intensity} --type pos --output_dir Output/3.MetaboliteQuantification/Metabolite2quant/pos 
        
        """
rule MetaboliteComparison:
    input:
        neg_quant = "Output/metaX/neg/data/neg-quant.txt",
        pos_quant = "Output/metaX/pos/data/pos-quant.txt",
        intensity = "Output/metaX/neg/data/neg-norm-metaboAnalystInput.csv",
        annotation = "Output/2.MetaboliteIdentification/identification.norm.intensity.csv",
    output:
        out="Output/4.MetaboliteComparison/comparison.stat.csv"
    params:
        compound = f"{ProgramDir}/bin/database/compounds.json",
        pvalue=0.05,
        log2ratio=1,
        pvalue_type='t.test_p.value',
        vip=1,
        mode=mode
    shell:
        """
        Rscript {ProgramDir}/bin/script/4.MetaboliteComparison/MetaboliteComparison.R --quant_neg {input.neg_quant} --quant_pos {input.pos_quant} --intensity {input.intensity} --sample_file Output/metaX/neg/sampleList.txt --annotation {input.annotation} --output_dir Output/4.MetaboliteComparison --compound {params.compound} --pvalue {params.pvalue} --log2ratio {params.log2ratio} --pvalue_type {params.pvalue_type} --vip 1 --mode {params.mode} --kgml {ProgramDir}/bin/database/KGML 
        """

rule MultiGroup:
    input:
        annotation = "Output/2.MetaboliteIdentification/identification.norm.intensity.csv",
        intensity = "Output/metaX/neg/data/neg-norm-metaboAnalystInput.csv",
    output:
        out="Output/5.MultiGroup/multigroup.significant.stat.csv"
    params:
        compound = f"{ProgramDir}/bin/database/compounds.json",
        pvalue=0.05,
        vip=1,
        mode=mode,
    shell:
        """
        Rscript {ProgramDir}/bin/script/5.MultiGroup/multigroup_plot.R --annotation {input.annotation} --intensity {input.intensity} --sample_file Output/metaX/neg/sampleList.txt --output_dir Output/5.MultiGroup --compound {params.compound} --pvalue {params.pvalue} --vip {params.vip} --mode {params.mode} --kgml {ProgramDir}/bin/database/KGML 
        """

rule moveFile:
    input:
        intensity_norm_neg = "Output/metaX/neg/data/neg-norm-metaboAnalystInput.csv",
        intensity_norm_pos = "Output/metaX/pos/data/pos-norm-metaboAnalystInput.csv",
        pos= "Output/3.MetaboliteQuantification/Feautre2quant/pos/sample.correlation.png",
    output:
        out1="Output/3.MetaboliteQuantification/neg-norm-metaboAnalystInput.csv",
        out2="Output/3.MetaboliteQuantification/pos-norm-metaboAnalystInput.csv",
    params:
        source_dir="Output/metaX/",
    shell:
        """
        cp {input.intensity_norm_neg} Output/3.MetaboliteQuantification/ && 
        cp {input.intensity_norm_pos} Output/3.MetaboliteQuantification/ &&
        cp {params.source_dir}/neg/data/neg-norm-peakCV.* Output/3.MetaboliteQuantification/Feautre2quant/neg/ &&
        cp {params.source_dir}/pos/data/pos-norm-peakCV.* Output/3.MetaboliteQuantification/Feautre2quant/pos/
        """

rule report:
    input:
        "Output/5.MultiGroup/multigroup.significant.stat.csv",
        "Output/4.MetaboliteComparison/comparison.stat.csv",
        "Output/3.MetaboliteQuantification/pos-norm-metaboAnalystInput.csv",
    output:
        "Output/report.html"
    shell:
        """
        cp {ProgramDir}/bin/script/report.rmd {ProgramDir}/bin/script/style.css . &&
        R -e "rmarkdown::render('report.rmd', output_file='{output}')" &&
        rm report.rmd style.css
        """

