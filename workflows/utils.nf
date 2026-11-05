process PREP_QUERIES {
    publishDir "${params.results_dir}/inputs", mode: 'copy', pattern: 'queries.tsv'

    input:
    path(inputs, stageAs: 'inputs_in/*')

    output:
    path('queries/*.fasta'), emit: fastas
    path('queries.tsv'),     emit: table

    script:
    """
    clust_assess.py prep-queries --in ${inputs} --outdir queries --table queries.tsv \
        ${params.queries_circular ? '' : '--no-circular'}
    """

    stub:
    """
    mkdir -p queries
    printf 'query_id\\tinput_id\\tdescription\\tsource_file\\tlength\\tstatus\\tnote\\n' > queries.tsv
    for n in \$(grep -h '^>' ${inputs} | cut -c2- | cut -d' ' -f1); do
        printf '>Q__%s\\nACGT\\n' "\$n" > "queries/Q__\$n.fasta"
        printf 'Q__%s\\t%s\\t\\tstub\\t4\\tok\\t\\n' "\$n" "\$n" >> queries.tsv
    done
    """
}

process VERSIONS {
    publishDir "${params.results_dir}/runningstat", mode: 'copy'

    output:
    path('software_versions.yml')

    script:
    """
    v() { "\$@" 2>&1 | grep -oE '[0-9]+(\\.[0-9]+)+[^ ]*' | head -1; }
    accio_commit() {
        local src=\$(python -c 'import importlib.util as u, os; print(os.path.dirname(os.path.dirname(u.find_spec("accio").origin)))' 2>/dev/null)
        git -C "\$src" rev-parse --short HEAD 2>/dev/null \\
        || git -C ${projectDir}/tools/accio rev-parse --short HEAD 2>/dev/null \\
        || echo unknown
    }
    {
        echo "clust_assess: ${params.version}"
        echo "accio: \$(python -c 'import importlib.metadata as m; print(m.version("accio"))' 2>/dev/null || echo NA) (commit \$(accio_commit))"
        echo "pling: \$(v pling --version)"
        echo "plasnet: \$(python -c 'import importlib.metadata as m; print(m.version("plasnet"))' 2>/dev/null || echo NA)"
        echo "skani: \$(v skani --version)"
        echo "mash: \$(v mash --version)"
        echo "nucmer: \$(v nucmer --version) (\$(command -v nucmer))"
        echo "blastn: \$(v blastn -version)"
        echo "mob_typer: \$(v mob_typer --version)"
        echo "amrfinder: \$(v amrfinder --version)"
        echo "nextflow: ${nextflow.version}"
    } > software_versions.yml
    """

    stub:
    """
    echo "clust_assess: ${params.version} (stub)" > software_versions.yml
    """
}
