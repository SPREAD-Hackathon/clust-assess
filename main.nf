#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { BUILD_DB } from './workflows/build_db.nf'
include { ASSESS   } from './workflows/assess.nf'
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
}
