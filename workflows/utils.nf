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

// Describes the database's own Pling run: within-subcommunity / within-community distance distributions
// (bin/distance_distribution.py) and the extra figures (bin/distribution_insights.py). Descriptive only:
// nothing downstream reads it, and a failure here does not stop the run (errorStrategy 'ignore' in nextflow.config).
process DB_DISTRIBUTIONS {
    publishDir "${params.results_dir}", mode: 'copy'

    input:
    path(pling, stageAs: 'db_pling')

    output:
    path('db_distributions'), emit: dir

    script:
    """
    distance_distribution.py --pling_dir db_pling --outdir db_distributions --plots \\
        --dcj_threshold ${params.dcj} --containment_threshold ${params.db_containment} > distributions.log 2>&1 \\
        || { tail -30 distributions.log >&2; exit 1; }
    distribution_insights.py --dist_dir db_distributions >> distributions.log 2>&1 \\
        || { tail -30 distributions.log >&2; exit 1; }
    mv distributions.log db_distributions/
    """

    stub:
    """
    mkdir -p db_distributions
    printf 'subcommunity\\tn_members\\tn_pairs\\n' > db_distributions/subcommunity_summary.tsv
    echo '{"stub": true}' > db_distributions/run_info.json
    touch db_distributions/within_pairs.tsv.gz
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
