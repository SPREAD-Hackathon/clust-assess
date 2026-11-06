# clust-assess: pipeline flow

How the Nextflow workflow (`main.nf`) runs, from input FASTA files to the published results. See [README_workflow.md](README_workflow.md) for parameters, outputs and how to read the verdict.

- `--mode assess` uses an existing database (`--db`); `--mode build_db` only builds one from `--refs`; `--mode all` builds it and then assesses the inputs against it.
- `DB_DISTRIBUTIONS` and `VERSIONS` run in every mode.
- Inputs that ACCIO cannot place still reach `REPORT`, as `no_assignment`.

```mermaid
flowchart TB
    %% ---------------- inputs ----------------
    IN["--input<br/>plasmid FASTA files"]
    REFS["--refs<br/>reference plasmids"]
    DBIN[("--db<br/>existing ACCIO / Pling database")]

    %% ---------------- BUILD_DB ----------------
    subgraph BUILD["BUILD_DB (--mode build_db or all)"]
        PR["PREP_REFS<br/>mark circular"] --> AB["ACCIO_BUILD<br/>accio create: skani dedup,<br/>Pling communities and subcommunities"]
    end
    REFS --> PR
    AB --> DB[("database<br/>plasmidDB.fasta, plasmidDB_info.csv,<br/>Mash sketch, pling/")]
    DBIN -- "assess mode" --> DB

    %% ---------------- ASSESS ----------------
    subgraph ASSESS["ASSESS (--mode assess or all)"]
        PQ["PREP_QUERIES<br/>split into one FASTA per plasmid,<br/>rename to Q__id, length + circular QC"]
        VD["VERIFY_DB<br/>required files, MANIFEST sha256"]
        AC["ACCIO_CALL, per input<br/>accio find: BLAST, Mash, synteny, length<br/>= chosen DB plasmid + subcommunity"]
        GRP["keep status = called,<br/>group inputs by subcommunity"]
        SM["SELECT_MEMBERS, per subcommunity<br/>DB plasmids of the subcommunity,<br/>drop copies of any input,<br/>rank by Mash to the inputs, keep k closest"]
        PM["PLING_MEMBERS, per subcommunity<br/>Pling on baseline members + its inputs"]
        PI["PLING_INPUTS<br/>Pling on all inputs"]
        RE["REPORT<br/>medians, nearest DB plasmid,<br/>resampling p_tighter, status, verdict"]

        PQ --> AC
        VD --> AC
        AC --> GRP --> SM --> PM
        VD --> SM
        PQ -- "all inputs" --> SM
        PQ -- "all inputs" --> PI
        AC -- "all calls, incl. no_assignment" --> RE
        PM -- "member_pairs, input_db_pairs" --> RE
        PI -- "input_pairs" --> RE
        PQ -- "queries.tsv" --> RE
        VD -- "db_info.json" --> RE
    end
    IN --> PQ
    DB --> VD
    DB --> AC
    DB --> SM

    %% ---------------- every mode ----------------
    DD["DB_DISTRIBUTIONS (--db_distributions, on by default)<br/>distance_distribution.py + distribution_insights.py<br/>on db/pling or --db_pling_dir"]
    VER["VERSIONS"]
    DB --> DD

    %% ---------------- published results ----------------
    OUT["final_output.csv, final_output_extended.csv,<br/>inputs.tsv, subcommunities.tsv, clust_assess.json"]
    MID["accio/, baseline/, inputs/, db_check/<br/>ACCIO reports, members + pair tables, input pairs"]
    ODD["db_distributions/<br/>summaries, within_pairs.tsv.gz, histograms,<br/>insight figures, DCJ + containment heatmaps,<br/>insight_zscores.png, zscore_thresholds.tsv"]
    OV["runningstat/<br/>software_versions.yml, Nextflow reports"]
    RE --> OUT
    AC --> MID
    PM --> MID
    PI --> MID
    DD --> ODD
    VER --> OV

    classDef store fill:#cde2fb,stroke:#256abf,color:#0b0b0b
    classDef out fill:#fcfcfb,stroke:#898781,color:#0b0b0b
    class DB,DBIN store
    class OUT,MID,ODD,OV out
```

| Step | Runs | What it does |
|---|---|---|
| `PREP_REFS`, `ACCIO_BUILD` | once, `build_db` / `all` | Build the ACCIO database: de-duplicate with skani, cluster with Pling |
| `VERIFY_DB` | once | Check the database files and their checksums |
| `PREP_QUERIES` | once | One FASTA per input plasmid, renamed `Q__<id>` |
| `ACCIO_CALL` | per input | Closest database plasmid and its Pling subcommunity |
| `SELECT_MEMBERS` | per subcommunity | The k database plasmids closest to the inputs (Mash), without copies of the inputs |
| `PLING_MEMBERS` | per subcommunity | Pling distances among those plasmids, and from each input to them |
| `PLING_INPUTS` | once | Pling distances among all inputs |
| `REPORT` | once | Inputs against the baseline: medians, resampling test, status and verdict |
| `DB_DISTRIBUTIONS` | once, every mode | Distance distributions and figures of the database's own Pling run |
| `VERSIONS` | once, every mode | Tool versions |
