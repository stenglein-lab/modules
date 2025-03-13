#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { PREPROCESS_READS } from './subworkflows/stenglein-lab/preprocess_reads.nf'

//
// WORKFLOW: Run main analysis pipeline
//
workflow MAIN_WORKFLOW {
  PREPROCESS_READS(params.fastq_dir, params.fastq_pattern, params.collapse_duplicate_reads)
}

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    MAIN_WORKFLOW ()
}
