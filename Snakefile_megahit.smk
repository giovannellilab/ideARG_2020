configfile: 'config.yaml'


rule all:
    input:
        # expand("{wdir}/{sample}/fastp", sample=config["SAMPLES"], wdir=config["WDIR"]),
        # expand("{wdir}/{sample}/fastp/R1.fastq.gz", sample=config["SAMPLES"], wdir=config["WDIR"]),
        # expand("{wdir}/{sample}/fastp/R2.fastq.gz", sample=config["SAMPLES"], wdir=config["WDIR"]),
        # expand("{wdir}/{sample}/megahit", sample=config["SAMPLES"], wdir=config["WDIR"]),
        # expand("{wdir}/{sample}/metaquast", sample=config["SAMPLES"], wdir=config["WDIR"]),
        # expand("{wdir}/rgi_card_db", wdir=config["WDIR"]),
        # expand("{wdir}/{sample}/rgi_main", wdir=config["WDIR"], sample=config["SAMPLES"]),
        expand("{wdir}/{sample}/prodigal", sample=config["SAMPLES"], wdir=config["WDIR"]),
        expand("{wdir}/{sample}/hmms_presence", sample=config["SAMPLES"], wdir=config["WDIR"]),
        # expand("{wdir}/{sample}/amr_finderplus", wdir=config["WDIR"], sample=config["SAMPLES"])


# rule run_fastp:
#     input:
#         r1="{wdir}/{sample}/R1.fastq.gz",
#         r2="{wdir}/{sample}/R2.fastq.gz",
#     output:
#         dir=directory("{wdir}/{sample}/fastp"),
#         r1="{wdir}/{sample}/fastp/R1.fastq.gz", 
#         r2="{wdir}/{sample}/fastp/R2.fastq.gz"
#     threads: 3
#     shell:
#         """fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} \
#         --unpaired1 {output.dir}/unpaired_R1.fastq.gz \
#         --unpaired2 {output.dir}/unpaired_R2.fastq.gz \
#         --json {output.dir}/report.json \
#         --html {output.dir}/report.html \
#         --failed_out {output.dir}/failed.fastq.gz \
#         --thread {threads} -q 20 u 40"""

# rule run_megahit:
#     input:
#         r1="{wdir}/{sample}/fastp/R1.fastq.gz",
#         r2="{wdir}/{sample}/fastp/R2.fastq.gz",
#     output:
#         directory("{wdir}/{sample}/megahit")
#     threads: 30
#     benchmark: "{wdir}/benchmark/{sample}/megahit.txt"
#     params:
#         complex_metagenome="--presets meta-large"
#     run:
#         shell("megahit {params.complex_metagenome} -t {threads} --memory 0.8 -1 {input.r1} -2 {input.r2} -o {output}")
#         shell("mv {output}/final.contigs.fa {output}/contigs.fasta")


# rule run_metaquast:
#     input:
#         contig_path="{wdir}/{sample}/megahit",
#     output:
#         directory("{wdir}/{sample}/metaquast")
#     threads: 3
#     params:
#         extra="--contig-thresholds 0,1000,10000,100000,1000000"
#     shell:
#         """
#         quast --threads {threads} {params.extra} -o {output} {input.contig_path}/contigs.fasta
#         """

# rule load_RGI_db:
#     output:
#         folder=directory("{wdir}/rgi_card_db")
#     params:
#         card_data="https://card.mcmaster.ca/latest/data",
#         variants_data="https://card.mcmaster.ca/latest/variants"
#     shell:
#         """
#         mkdir -p {output.folder} \
#           && rgi clean \
#           && ( cd {output.folder} \
#           && wget {params.card_data} \
#           && tar -xvf data ./card.json \
#           && wget -O wildcard_data.tar.bz2 {params.variants_data} \
#           && mkdir -p wildcard \
#           && tar -xjf wildcard_data.tar.bz2 -C wildcard \
#           && gunzip wildcard/*.gz \
#           && rgi card_annotation -i card.json > card_annotation.log 2>&1 \
#           && temp_var=$(ls *_all.fasta | grep "_v" | head -n 1) \
#           && IFS="_" read b1 b2 version remaining <<< $temp_var \
#           && version_rgi_db=$(echo $version | cut -c2-) \
#           && rgi wildcard_annotation -i wildcard --card_json card.json -v $version_rgi_db > wildcard_annotation.log 2>&1 \
#           && rgi load \
#             --card_json card.json \
#             --debug \
#             --card_annotation card_database_v${{version_rgi_db}}.fasta \
#             --card_annotation_all_models card_database_v${{version_rgi_db}}_all.fasta \
#             --wildcard_annotation wildcard_database_v${{version_rgi_db}}.fasta \
#             --wildcard_annotation_all_models wildcard_database_v${{version_rgi_db}}_all.fasta \
#             --wildcard_index wildcard/index-for-model-sequences.txt \
#             --wildcard_version ${{version_rgi_db}} \
#             --amr_kmers wildcard/all_amr_61mers.txt \
#             --kmer_database wildcard/61_kmer_db.json \
#             --kmer_size 61
#           )
#         """

# rule run_rgi_main:
#     input:
#         contig_path="{wdir}/{sample}/megahit",
#     output:
#         folder=directory("{wdir}/{sample}/rgi_main")
#     threads: 34
#     log: "{wdir}/{sample}/rgi_main/log.out"
#     benchmark: "{wdir}/benchmark/rgi_main/{sample}_bench.txt"
#     params:
#         outfile="{wdir}/{sample}/rgi_main/{sample}"
#     run:
#         shell("mkdir -p {output.folder}")
#         shell("touch {log}")
#         shell("rgi main --input_sequence {input.contig_path}/contigs.fasta \
#             --output_file {params.outfile} \
#             --input_type contig \
#             --num_threads {threads} \
#             --alignment_tool DIAMOND \
#             --orf_finder PRODIGAL \
#             --low_quality \
#             --include_nudge \
#             --split_prodigal_jobs \
#             --clean >> {log} 2>&1")

# rule amrfinderplus:
#     input:
#         contig_path="{wdir}/{sample}/megahit"
#     output:
#         folder=directory("{wdir}/{sample}/amr_finderplus")
#     threads: 5
#     benchmark: "{wdir}/benchmark/{sample}/amr_finderplus.txt"
#     log: "{wdir}/{sample}/amr_finderplus/log.out"
#     run:
#         shell("mkdir -p {output.folder}")
#         shell("amrfinder -n {input.contig_path}/contigs.fasta \
#             --threads {threads} \
#             --plus \
#             -o {log}")


## GeoMosaic conda environment
rule run_prodigal:
    input:
        contig_path="{wdir}/{sample}/megahit",
    output:
        directory("{wdir}/{sample}/prodigal")
    params:
        extra="-p meta",
        quiet="-q"
    run:
        shell("mkdir -p {output}")
        shell("prodigal -i {input.contig_path}/contigs.fasta -o {output}/genes.out -a {output}/protein_translations.faa {params.extra} {params.quiet}")

        from geomosaic.parsing_output.prodigal_orf_mapping import parsing_prodigal_orfs

        fasta_input = f"{output}/protein_translations.faa"
        output_mapping = f"{output}/orf_contig_mapping.tsv"
        output_fasta = f"{output}/orf_predicted.faa"
        output_simple_mapping = f"{output}/simple_orf_contig_mapping.tsv"

        parsing_prodigal_orfs(fasta_input, output_mapping, output_fasta, output_simple_mapping)


rule hmms_presence:
    input:        
        prodigal_path="{wdir}/{sample}/prodigal",
    output:
        folder=directory("{wdir}/{sample}/hmms_presence")
    params:
        hmm_folder=config["hmm_folder"]
    threads: 1
    run:
        shell("mkdir -p {output.folder}")
        
        import pandas as pd

        df_mapping = pd.read_csv(os.path.join(str(input.prodigal_path), "simple_orf_contig_mapping.tsv"), sep="\t")
        local_sample = output.folder.split("/")[-2]
        
        results_filename = f"{output.folder}/presence_results.tsv"
        with open(results_filename, "wt") as fd:
            fd.write("name\tmatch\tsample\n")

            for hmm in os.listdir(params.hmm_folder):
                if not hmm.endswith('.hmm'):
                    continue
                
                filename=hmm.split(".hmm")[0]
                out_path=os.path.join(output.folder, filename)
                shell("mkdir -p {out_path}")

                hmm_file=os.path.join(params.hmm_folder, hmm)
                shell("hmmsearch --tblout {out_path}/hmmsearch_table.txt -E 0.00001 -o /dev/null --cpu {threads} --notextw {hmm_file} {input.prodigal_path}/orf_predicted.faa") 
                shell("grep -v \"^#\" {out_path}/hmmsearch_table.txt > {out_path}/results.txt || true")
                if os.stat(os.path.join(out_path, "results.txt")).st_size == 0:
                    continue

                shell("awk '{{print $1\"\t\"$3}}' {out_path}/results.txt > {out_path}/hits.txt")

                df_hits = pd.read_csv(os.path.join(out_path,"hits.txt"), sep="\t", names=["orf_id", "name"])
                df_hits.drop_duplicates(inplace=True)

                m1 = df_hits.merge(df_mapping, on="orf_id", how="left")

                subset = m1.loc[:, ["orf_id", "name", "contig"]]
                name=list(subset.name.unique())[0]
                match_on_contigs = df_hits.shape[0]

                fd.write(f"{name}\t{match_on_contigs}\t{local_sample}\n")
                
                subset.to_csv(os.path.join(out_path, "hmm_presence.tsv"), sep="\t", header=True, index=False)
