include { PREP_QUERIES } from './utils.nf'

process VERIFY_DB {
    tag "${db.name}"
    publishDir "${params.results_dir}/db_check", mode: 'copy', pattern: 'db_info.json'

    input:
    path(db)

    output:
    path('db_info.json'), emit: info

    script:
    """
    clust_assess.py verify-db --db ${db} --out db_info.json
    """

    stub:
    """
    echo '{"db_name": "${db.name}", "manifest": "stub"}' > db_info.json
    """
}

process ACCIO_CALL {
    tag "${qid}"
    publishDir "${params.results_dir}/accio", mode: 'copy', pattern: '*/*_plasmids_chosen.csv'

    input:
    tuple val(qid), path(fasta)
    path(db)
    path(_verified)

    output:
    tuple val(qid), path(fasta), path("${qid}.call.tsv"), emit: call
    path("${qid}"),                                       emit: raw

    script:
    """
    nucmer --version 2>&1 | grep -qE '^4\\.' || { echo "nucmer is not MUMmer 4: \$(command -v nucmer)" >&2; exit 1; }
    set +e
    accio find ${fasta} -d ${db} -o ${qid} -t ${task.cpus} --sample_name ${qid} ${params.accio_find_args} > accio_stdout.log 2>&1
    status=\$?
    set -e
    mkdir -p ${qid} && mv accio_stdout.log ${qid}/
    clust_assess.py parse-accio --query-id ${qid} --dir ${qid} --exit-status \$status --out ${qid}.call.tsv
    """

    stub:
    """
    mkdir -p ${qid}
    st=called; sc=community_0_subcommunity_0
    case "${qid}" in *nocall*) st=no_assignment; sc= ;; *other*) sc=community_1_subcommunity_0 ;; esac
    printf 'query_id\\tstatus\\treason\\tref_plasmid\\tsubcommunity\\tcommunity\\tis_hub\\tn_calls\\tfinal_score\\tmash_score\\tblast_identity\\tblast_coverage\\tsynteny_score\\tlength_score\\tref_length\\tref_rep_type\\tall_calls\\n' > ${qid}.call.tsv
    printf '%s\\t%s\\t\\tm1\\t%s\\t%s\\tFalse\\t1\\t0.99\\t0.999\\t\\t\\t\\t\\t\\t\\t\\n' "${qid}" "\$st" "\$sc" "\${sc%_subcommunity_*}" >> ${qid}.call.tsv
    """
}

process SELECT_MEMBERS {
    tag "${sc}"
    publishDir "${params.results_dir}/baseline", mode: 'copy', pattern: '*.members.tsv'

    input:
    tuple val(sc), val(qids), path(fastas)
    path(db)
    path(_verified)

    output:
    tuple val(sc), path("${sc}.members.tsv"), path('members/*.fasta'), emit: members

    script:
    """
    clust_assess.py candidates --db ${db} --subcommunity '${sc}' --mode ${params.baseline_members} \
        --out candidates.fasta --table candidates.tsv
    mash sketch -i -o candidates candidates.fasta > /dev/null 2>&1
    mash dist candidates.msh ${fastas} > mash.tsv
    clust_assess.py select-members --candidates candidates.tsv --candidates-fasta candidates.fasta --mash mash.tsv \
        --queries ${qids.join(' ')} --k ${params.k_closest} --outdir members --table ${sc}.members.tsv
    """

    stub:
    """
    mkdir -p members
    printf 'plasmid_id\\tlength\\tpling_type\\tcommunity\\tmin_mash_to_inputs\\trank\\tselected\\n' > ${sc}.members.tsv
    for m in m1 m2 m3; do
        printf '>%s\\nACGT\\n' \$m > members/\$m.fasta
        printf '%s\\t50000\\t${sc}\\tc\\t0.001\\t1\\tTrue\\n' \$m >> ${sc}.members.tsv
    done
    """
}

process PLING_MEMBERS {
    tag "${sc}"
    publishDir "${params.results_dir}/baseline", mode: 'copy', pattern: '*.member_pairs.tsv'

    input:
    tuple val(sc), path(members_tsv), path(member_fastas, stageAs: 'members/*')

    output:
    tuple path(members_tsv), path("${sc}.member_pairs.tsv"), emit: baseline

    script:
    """
    for f in members/*.fasta; do echo "\$PWD/\$f"; done > genomes.txt
    if [ \$(wc -l < genomes.txt) -ge 2 ]; then
        pling cluster align genomes.txt pling_out --containment_distance ${params.containment} --dcj ${params.dcj} \
            --cores ${task.cpus} ${params.pling_args} > pling.log 2>&1 || { tail -30 pling.log >&2; exit 1; }
    fi
    clust_assess.py pairs --pling-dir pling_out --genomes genomes.txt --out ${sc}.member_pairs.tsv
    """

    stub:
    """
    printf 'a\\tb\\tcontainment\\tdcj\\tdcj_censored\\nm1\\tm2\\t0.2\\t3\\tFalse\\nm1\\tm3\\t0.3\\t5\\tFalse\\nm2\\tm3\\t0.25\\t4\\tFalse\\n' > ${sc}.member_pairs.tsv
    """
}

process PLING_INPUTS {
    publishDir "${params.results_dir}/inputs", mode: 'copy', pattern: 'input_pairs.tsv'

    input:
    path(fastas, stageAs: 'inputs/*')

    output:
    path('input_pairs.tsv'), emit: pairs

    script:
    """
    for f in inputs/*.fasta; do echo "\$PWD/\$f"; done > genomes.txt
    if [ \$(wc -l < genomes.txt) -ge 2 ]; then
        pling cluster align genomes.txt pling_out --containment_distance ${params.containment} --dcj ${params.dcj} \
            --cores ${task.cpus} ${params.pling_args} > pling.log 2>&1 || { tail -30 pling.log >&2; exit 1; }
    fi
    clust_assess.py pairs --pling-dir pling_out --genomes genomes.txt --out input_pairs.tsv
    """

    stub:
    """
    ids=\$(for f in inputs/*.fasta; do basename "\$f" .fasta; done)
    printf 'a\\tb\\tcontainment\\tdcj\\tdcj_censored\\n' > input_pairs.tsv
    for a in \$ids; do for b in \$ids; do
        if [ "\$a" \\< "\$b" ]; then printf '%s\\t%s\\t0.1\\t1\\tFalse\\n' "\$a" "\$b" >> input_pairs.tsv; fi
    done; done
    """
}

process REPORT {
    publishDir "${params.results_dir}", mode: 'copy'

    input:
    path(queries_tsv)
    path(calls, stageAs: 'calls/*')
    path(input_pairs)
    path(baselines, stageAs: 'baseline/*')
    path(db_info)

    output:
    path('final_output.csv')
    path('final_output_extended.csv')
    path('inputs.tsv')
    path('subcommunities.tsv')
    path('clust_assess.json')

    script:
    def p = groovy.json.JsonOutput.toJson([containment: params.containment, dcj: params.dcj, k_closest: params.k_closest,
                                           n_resample: params.n_resample, alpha: params.alpha, seed: params.seed,
                                           db: params.mode == 'assess' ? params.db : "${params.db_outdir}/${params.db_name}"])
    """
    echo '${p}' > params.json
    clust_assess.py report --queries ${queries_tsv} --calls calls/* --input-pairs ${input_pairs} \
        --baselines \$(ls baseline/* 2>/dev/null) --params params.json --outdir .
    """

    stub:
    """
    echo '{}' > params.json
    clust_assess.py report --queries ${queries_tsv} --calls calls/* --input-pairs ${input_pairs} \
        --baselines \$(ls baseline/* 2>/dev/null) --params params.json --outdir .
    """
}

// Interactive HTML report (html_report/): one self-contained report.html built from final_output_extended.csv.
// The reference is rebuilt from this run's baseline/*.member_pairs.tsv (censored pairs dropped); if that
// gives nothing, the report uses the reference bundled in html_report/reference. Optional: a failure here
// never stops the run (errorStrategy 'ignore' in nextflow.config). Switch off with --html_report false.
process HTML_REPORT {
    publishDir "${params.results_dir}", mode: 'copy'

    input:
    path(final_output_ext)
    path(baselines, stageAs: 'baseline/*')
    path(report_dir, stageAs: 'html_report')

    output:
    path('report.html')

    script:
    """
    printf 'plasmid_1\\tplasmid_2\\tsubcommunity\\tdcj\\n' > reference.tsv
    for f in baseline/*.member_pairs.tsv; do
      [ -e "\$f" ] || continue
      sc=\$(basename "\$f" .member_pairs.tsv)
      awk -F'\\t' -v sc="\$sc" 'NR>1 && \$4 ~ /^[0-9]+(\\.[0-9]+)?\$/ && \$5!="True" {print \$1"\\t"\$2"\\t"sc"\\t"\$4}' "\$f" >> reference.tsv
    done
    if [ "\$(wc -l < reference.tsv)" -gt 1 ]; then
      CA_NOOPEN=1 sh html_report/make_report.sh ${final_output_ext} report.html reference.tsv
    else
      CA_NOOPEN=1 sh html_report/make_report.sh ${final_output_ext} report.html
    fi
    """

    stub:
    """
    touch report.html
    """
}

workflow ASSESS {
    take:
    inputs
    db

    main:
    VERIFY_DB(db)
    PREP_QUERIES(inputs)

    def queries = PREP_QUERIES.out.fastas.flatten().map { f -> tuple(f.baseName, f) }
    ACCIO_CALL(queries, db, VERIFY_DB.out.info)

    def by_subcommunity = ACCIO_CALL.out.call
        .map { qid, fa, tsv ->
            def row = tsv.splitCsv(header: true, sep: '\t')[0]
            tuple(row.status, row.subcommunity, qid, fa)
        }
        .filter { st, _sc, _qid, _fa -> st == 'called' }
        .map { _st, sc, qid, fa -> tuple(sc, [qid, fa]) }
        .groupTuple()
        .map { sc, items ->
            def sorted = items.sort { a, b -> a[0] <=> b[0] }
            tuple(sc, sorted.collect { x -> x[0] }, sorted.collect { x -> x[1] })
        }

    SELECT_MEMBERS(by_subcommunity, db, VERIFY_DB.out.info)
    PLING_MEMBERS(SELECT_MEMBERS.out.members)
    PLING_INPUTS(PREP_QUERIES.out.fastas)

    REPORT(PREP_QUERIES.out.table,
           ACCIO_CALL.out.call.map { _qid, _fa, tsv -> tsv }.collect(),
           PLING_INPUTS.out.pairs,
           PLING_MEMBERS.out.baseline.flatten().collect().ifEmpty([]),
           VERIFY_DB.out.info)

    if (params.html_report) {
        HTML_REPORT(REPORT.out[1],
                    PLING_MEMBERS.out.baseline.flatten().collect().ifEmpty([]),
                    file("${projectDir}/html_report", checkIfExists: true))
    }
}
