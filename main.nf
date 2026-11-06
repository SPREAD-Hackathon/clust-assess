#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { BUILD_DB } from './workflows/build_db.nf'
include { ASSESS; HTML_REPORT } from './workflows/assess.nf'
include { VERSIONS; DB_DISTRIBUTIONS } from './workflows/utils.nf'

def require(names) {
    names.each { n ->
        if (!params[n]) {
            error "clust-assess: --${n} is required for --mode ${params.mode}"
        }
    }
}

workflow {
    if (!(params.mode in ['build_db', 'assess', 'all'])) {
        error "clust-assess: --mode must be build_db, assess or all (got '${params.mode}')"
    }
    log.info """
    clust-assess ${params.version}  mode=${params.mode}
      input       : ${params.input}
      db          : ${params.mode == 'assess' ? params.db : "${params.db_outdir}/${params.db_name} (built from ${params.refs})"}
      thresholds  : containment ${params.containment}, DCJ ${params.dcj}; k_closest ${params.k_closest}
      results_dir : ${params.results_dir}
    """.stripIndent()

    VERSIONS()

    if (params.mode == 'build_db') {
        require(['refs'])
    }
    else if (params.mode == 'assess') {
        require(['input', 'db'])
    }
    else {
        require(['refs', 'input'])
    }

    def db = params.mode == 'assess'
        ? channel.value(file(params.db, checkIfExists: true))
        : BUILD_DB(channel.fromPath(params.refs, checkIfExists: true).collect()).db
    if (params.mode != 'build_db') {
        ASSESS(channel.fromPath(params.input, checkIfExists: true).collect(), db)
    }
    // Distance distributions of the database's own Pling run (default <db>/pling), once per run
    if (params.db_distributions) {
        DB_DISTRIBUTIONS(db.map { d -> params.db_pling_dir ? file(params.db_pling_dir, checkIfExists: true) : d.resolve('pling') })
    }
    // Interactive HTML report: the last step. Its z-score background is the database's own pairs when DB_DISTRIBUTIONS ran;
    // the empty assets/NO_FILE stands in when it is off or failed (the report then uses the pairs in the results file).
    if (params.mode != 'build_db' && params.html_report) {
        def no_file = file("${projectDir}/assets/NO_FILE", checkIfExists: true)
        def within = params.db_distributions
            ? DB_DISTRIBUTIONS.out.dir.map { d -> d.resolve('within_pairs.tsv.gz') }.ifEmpty(no_file)
            : channel.value(no_file)
        HTML_REPORT(ASSESS.out.final_extended,
                    ASSESS.out.inputs_table,
                    within,
                    ASSESS.out.baselines,
                    file("${projectDir}/html_report", checkIfExists: true))
    }
}
