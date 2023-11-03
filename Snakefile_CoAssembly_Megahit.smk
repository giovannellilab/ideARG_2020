configfile: 'config_CoAssembly_Megahit.yaml'


rule all:
    input:
        # expand("{wdir}/rgi_card_db", wdir=config["WDIR"]),
        # expand("{wdir}/{sample}/rgi_main", wdir=config["WDIR"], sample=config["SAMPLES"]),
        # expand("{wdir}/megares_forblast", wdir=config["WDIR"]),
        # expand("{wdir}/{sample}/blastn_megares", wdir=config["WDIR"], sample=config["SAMPLES"]),
        # expand("{wdir}/{sample}/blastn_megares_culling", wdir=config["WDIR"], sample=config["SAMPLES"]),
        # expand("{wdir}/resfinder_forblast", wdir=config["WDIR"]),
        # expand("{wdir}/{sample}/blastn_resfinder", wdir=config["WDIR"], sample=config["SAMPLES"]),
        # expand("{wdir}/{sample}/blastn_resfinder_culling", wdir=config["WDIR"], sample=config["SAMPLES"]),
        # expand("{wdir}/{sample}/amr_finderplus", wdir=config["WDIR"], sample=config["SAMPLES"]),
        expand("{wdir}/{sample}/prodigal", sample=config["SAMPLES"], wdir=config["WDIR"]),
        # expand("{wdir}/deepArg_db", wdir=config["WDIR"]),
        expand("{wdir}/{sample}/deeparg", sample=config["SAMPLES"], wdir=config["WDIR"]),
        expand("{wdir}/{sample}/abricate", sample=config["SAMPLES"], wdir=config["WDIR"]),
        expand("{wdir}/abricate_summary", wdir=config["WDIR"]),
        expand("{wdir}/{sample}/hmms_presence_AMRFinderPlus", sample=config["SAMPLES"], wdir=config["WDIR"]),
        expand("{wdir}/{sample}/hmms_presence", sample=config["SAMPLES"], wdir=config["WDIR"]),


rule load_RGI_db:
    output: 
        folder=directory("{wdir}/rgi_card_db")
    params:
        card_data="https://card.mcmaster.ca/latest/data",
        variants_data="https://card.mcmaster.ca/latest/variants"
    run:
        shell("mkdir -p {output.folder}")
        shell("touch {output.folder}/rgi_auto_load.txt")
        shell("echo 'RGI AUTO LOAD' >> {output.folder}/rgi_auto_load.txt")
        shell("rgi clean")
        shell("rgi auto_load")
        

rule run_rgi_main:
    input:
        contig_path="{wdir}/{sample}/megahit",
    output:
        folder=directory("{wdir}/{sample}/rgi_main")
    threads: 17
    log: "{wdir}/{sample}/rgi_main/log.out"
    benchmark: "{wdir}/benchmark/rgi_main/{sample}_bench.txt"
    params:
        outfile="{wdir}/{sample}/rgi_main/{sample}"
    run:
        shell("mkdir -p {output.folder}")
        shell("touch {log}")
        shell("rgi main --input_sequence {input.contig_path}/contigs.fasta \
            --output_file {params.outfile} \
            --input_type contig \
            --num_threads {threads} \
            --alignment_tool DIAMOND \
            --orf_finder PRODIGAL \
            --low_quality \
            --include_nudge \
            --split_prodigal_jobs \
            --clean >> {log} 2>&1")


rule create_megares_db:
    input:
        megares_fasta="databases/megares_v3.00/megares_database_v3.00.fasta"
    output:        
        megares_blastdb=directory("{wdir}/megares_forblast")
    params:
        dbtype="nucl"
    run:
        shell("makeblastdb -in {input.megares_fasta} -dbtype {params.dbtype} -out {output.megares_blastdb}/megares_blast")


rule blastn_megares:
    input: 
        contig_path="{wdir}/{sample}/megahit",
        megares_blastdb="{wdir}/megares_forblast"
    output:
        folder=directory("{wdir}/{sample}/blastn_megares"),
    threads: 10
    params:
        outfmt="\"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp\"",
        outfile="{wdir}/{sample}/blastn_megares/{sample}_blast.out",
        evalue=1e-10,
        perc_identity=50
    run:
        shell("mkdir -p {output.folder}")
        shell("blastn \
            -query {input.contig_path}/contigs.fasta \
            -db {input.megares_blastdb}/megares_blast \
            -out {params.outfile} \
            -num_threads {threads} \
            -evalue {params.evalue} \
            -perc_identity {params.perc_identity} \
            -outfmt {params.outfmt}")


rule blastn_megares_culling:
    input: 
        contig_path="{wdir}/{sample}/megahit",
        megares_blastdb="{wdir}/megares_forblast"
    output:
        folder=directory("{wdir}/{sample}/blastn_megares_culling"),
    threads: 10
    params:
        outfmt="\"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp\"",
        outfile="{wdir}/{sample}/blastn_megares_culling/{sample}_blast.out",
        evalue=1e-10,
        perc_identity=50
    run:
        shell("mkdir -p {output.folder}")
        shell("blastn \
            -query {input.contig_path}/contigs.fasta \
            -db {input.megares_blastdb}/megares_blast \
            -out {params.outfile} \
            -num_threads {threads} \
            -evalue {params.evalue} \
            -culling_limit 1 \
            -perc_identity {params.perc_identity} \
            -outfmt {params.outfmt}")


rule create_resfinder_db:
    input:
        resfinder_fasta="databases/all_resfinder_db.fasta"
    output:        
        resfinder_blastdb=directory("{wdir}/resfinder_forblast")
    params:
        dbtype="nucl"
    run:
        shell("makeblastdb -in {input.resfinder_fasta} -dbtype {params.dbtype} -out {output.resfinder_blastdb}/resfinder_blast")


rule blastn_resfinder:
    input: 
        contig_path="{wdir}/{sample}/megahit",
        resfinder_blastdb="{wdir}/resfinder_forblast"
    output:
        folder=directory("{wdir}/{sample}/blastn_resfinder"),
    threads: 10
    params:
        outfmt="\"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp\"",
        outfile="{wdir}/{sample}/blastn_resfinder/{sample}_blast.out",
        evalue=1e-10,
        perc_identity=50
    run:
        shell("mkdir -p {output.folder}")
        shell("blastn \
            -query {input.contig_path}/contigs.fasta \
            -db {input.resfinder_blastdb}/resfinder_blast \
            -out {params.outfile} \
            -num_threads {threads} \
            -evalue {params.evalue} \
            -perc_identity {params.perc_identity} \
            -outfmt {params.outfmt}")


rule blastn_resfinder_culling:
    input: 
        contig_path="{wdir}/{sample}/megahit",
        resfinder_blastdb="{wdir}/resfinder_forblast"
    output:
        folder=directory("{wdir}/{sample}/blastn_resfinder_culling"),
    threads: 10
    params:
        outfmt="\"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp\"",
        outfile="{wdir}/{sample}/blastn_resfinder_culling/{sample}_blast.out",
        evalue=1e-10,
        perc_identity=50
    run:
        shell("mkdir -p {output.folder}")
        shell("blastn \
            -query {input.contig_path}/contigs.fasta \
            -db {input.resfinder_blastdb}/resfinder_blast \
            -out {params.outfile} \
            -num_threads {threads} \
            -evalue {params.evalue} \
            -culling_limit 1 \
            -perc_identity {params.perc_identity} \
            -outfmt {params.outfmt}")


rule amrfinderplus:
    input:
        contig_path="{wdir}/{sample}/megahit"
    output:
        folder=directory("{wdir}/{sample}/amr_finderplus")
    threads: 10
    benchmark: "{wdir}/benchmark/{sample}/amr_finderplus.txt"
    log: "{wdir}/{sample}/amr_finderplus/log.out"
    run:
        shell("mkdir -p {output.folder}")
        shell("amrfinder -n {input.contig_path}/contigs.fasta \
            --threads {threads} \
            --plus \
            -o {log}")


# ## GeoMosaic conda environment
rule run_prodigal:
    input:
        contig_path="{wdir}/{sample}/megahit",
    output:
        directory("{wdir}/{sample}/prodigal")
    params:
        extra="-p meta",
        quiet="-q"
    threads: 1
    run:
        shell("mkdir -p {output}")
        shell("prodigal -i {input.contig_path}/contigs.fasta -o {output}/genes.out -a {output}/protein_translations.faa {params.extra} {params.quiet}")

        from geomosaic.parsing_output.prodigal_orf_mapping import parsing_prodigal_orfs

        fasta_input = f"{output}/protein_translations.faa"
        output_mapping = f"{output}/orf_contig_mapping.tsv"
        output_fasta = f"{output}/orf_predicted.faa"
        output_simple_mapping = f"{output}/simple_orf_contig_mapping.tsv"

        parsing_prodigal_orfs(fasta_input, output_mapping, output_fasta, output_simple_mapping)


# rule deepArg_DB:
#     output:
#         folder=directory("{wdir}/deepArg_db")
#     params:
#         simg="{wdir}/../deeparg.simg"
#     run:
#         shell("mkdir -p {output.folder}")
#         shell("singularity exec --bind {output.folder}:/deepargDB/ {params.simg} deeparg download_data -o /deepargDB")


rule deepArg:
    input:
        prodigal_path="{wdir}/{sample}/prodigal",
        deeparg_db="/mnt/data/bigdata/ideARG_2020/results_metaspades/deepArg_db",
    output:
        folder=directory("{wdir}/{sample}/deeparg")
    threads: 30
    params:
        simg="{wdir}/../deeparg.simg",
        sample_output="{sample}_output"
    run:
        shell("mkdir -p {output.folder}")
        shell("singularity exec \
            --bind {output.folder}:/deeparg_output/ \
            --bind {input.deeparg_db}:/deepargDB/ \
            --bind {input.prodigal_path}:/prodigal_input/ \
            {params.simg} deeparg predict \
            --model LS \
            -i /prodigal_input/orf_predicted.faa \
            -o /deeparg_output/{params.sample_output} \
            -d /deepargDB \
            --type prot \
            --min-prob 0.8 \
            --arg-alignment-identity 50 \
            --arg-alignment-evalue 1e-10 \
            --arg-num-alignments-per-entry 1000")


rule hmms_presence:
    input:        
        prodigal_path="{wdir}/{sample}/prodigal",
    output:
        folder=directory("{wdir}/{sample}/hmms_presence")
    params:
        hmm_folder=config["hmm_folder"]
    threads: 5
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


# ## GeoMosaic conda environment
rule run_abricate:
    input:
        contig_path="{wdir}/{sample}/megahit",
    output:
        directory("{wdir}/{sample}/abricate")
    params:
        sample_name="{sample}"
    threads: 10
    conda: "abricate"
    shell:
        """mkdir -p {output}

        for db in ncbi vfdb plasmidfinder argannot card ecoli_vf megares ecoh resfinder; do
            abricate --db $db --threads {threads} {input.contig_path}/contigs.fasta > {output}/{params.sample_name}_$db.out
        done;    
        """


rule run_abricate_summary:
    input:
        wdir_path="{wdir}",
    output:
        directory("{wdir}/abricate_summary")
    conda: "abricate"
    shell:
        """mkdir -p {output}

        for db in ncbi vfdb plasmidfinder argannot card ecoli_vf megares ecoh resfinder; do
            all_db_files=$(ls {input}/**/abricate/*_$db.out | sort)

            abricate --summary $all_db_files > {output}/summary_$db.tab
        done;    
        """


rule hmms_presence_AMRFinderPlus:
    input:        
        prodigal_path="{wdir}/{sample}/prodigal",
    output:
        folder=directory("{wdir}/{sample}/hmms_presence_AMRFinderPlus")
    params:
        hmm_folder=config["hmm_folder_AMRFinderPlus"]
    threads: 5
    run:
        shell("mkdir -p {output.folder}")
        
        import pandas as pd

        df_mapping = pd.read_csv(os.path.join(str(input.prodigal_path), "simple_orf_contig_mapping.tsv"), sep="\t")
        local_sample = output.folder.split("/")[-2]
        
        results_filename = f"{output.folder}/presence_results.tsv"
        with open(results_filename, "wt") as fd:
            fd.write("name\tmatch\tsample\n")

            for hmm in os.listdir(params.hmm_folder):
                if not hmm.endswith('.HMM'):
                    continue
                
                filename=hmm.split(".HMM")[0]
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
