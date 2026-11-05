process PREP_REFS {
    publishDir "${params.results_dir}/build_db", mode: 'copy', pattern: 'refs.tsv'

    input:
    path(refs, stageAs: 'refs_in/*')

    output:
    path('refs.prepared.fasta'), emit: fasta
    path('refs.tsv'),            emit: table

    script:
    """
    clust_assess.py prep-refs --in ${refs} --out refs.prepared.fasta --table refs.tsv \
        ${params.refs_circular ? '' : '--no-circular'}
    """

    stub:
    """
    cat ${refs} > refs.prepared.fasta
    printf 'ref_id\\tinput_id\\tsource_file\\tlength\\n' > refs.tsv
    """
}

process ACCIO_BUILD {
    publishDir "${params.db_outdir}", mode: 'copy'

    input:
    path(refs_fasta)

    output:
    path("${params.db_name}"), emit: db

    script:
    """
    nucmer --version 2>&1 | grep -qE '^4\\.' || { echo "nucmer is not MUMmer 4: \$(command -v nucmer)" >&2; exit 1; }
    accio create ${refs_fasta} -o ${params.db_name} -t ${task.cpus} ${params.accio_create_args} > accio_create.log 2>&1 \
        || { tail -30 accio_create.log >&2; exit 1; }
    mv accio_create.log ${params.db_name}/
    clust_assess.py db-manifest --db ${params.db_name}
    """

    stub:
    """
    mkdir -p ${params.db_name}/pling
    touch ${params.db_name}/plasmidDB.fasta ${params.db_name}/plasmidDB.fasta.msh
    printf 'plasmid_id,length,pling_type,community\\nm1,50000,community_0_subcommunity_0,community_0\\n' > ${params.db_name}/plasmidDB_info.csv
    """
}

workflow BUILD_DB {
    take:
    refs

    main:
    PREP_REFS(refs)
    ACCIO_BUILD(PREP_REFS.out.fasta)

    emit:
    db = ACCIO_BUILD.out.db
}
