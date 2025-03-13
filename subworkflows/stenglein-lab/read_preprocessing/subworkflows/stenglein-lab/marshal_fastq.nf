include { CHECK_FASTQ_COMPRESSED   } from '../../modules/stenglein-lab/check_fastq_compressed/main'
include { COUNT_FASTQ              } from '../../modules/stenglein-lab/count_fastq/main'
include { SEQTK_SAMPLE             } from '../../modules/nf-core/seqtk/sample/main'

/*
 Create a channel with input fastq, possibly from multiple directories (specified as a comma-separated list)
 also returns IDs, determeind from filenames stripped of extensions and Illumina-added stuff like _001_ _S?_ etc

 Asumptions about input fastq:

 - gzip compressed
 - filenames end in .fastq.gz or .fq.gz
 - single-end or paired en OK
 - filenames matchable using param fastq_pattern
 - in one or more directories specified by param input_fastq_dir
 */
workflow MARSHAL_FASTQ {

 take:
 input_fastq_dir        // the path to a directory containing fastq file(s) or a comma-separated list of dirs
 fastq_pattern          // the regex that will be matched to identify fastq
 count_fastq            // boolean: should we count # of reads in input fastq?
 subsample_size         // null, or an integer value defining the # of reads to subsample each input fastq to
                        //       or an fractional value (0-1) defining the fraction of reads to subsample

 main:

  ch_versions         = Channel.empty()
  ch_fastq_counts     = Channel.empty()

  // User can specify multiple directories containing input fastq
  // In this case, the directories should be provided as a 
  // comma-separated list (no spaces)
  def fastq_dirs = input_fastq_dir.tokenize(',')

  // construct list of directories in which to find fastq
  fastq_dir_list = []
  for (dir in fastq_dirs){
     def file_pattern = "${dir}/${fastq_pattern}"
     fastq_dir_list.add(file_pattern)
  }

  /*
   These fastq files are the main input to this workflow
  */
  Channel
  .fromFilePairs(fastq_dir_list, size: -1, checkIfExists: true, maxDepth: 1)
  .map{ name, reads ->

         // define a new empty map named meta for each sample
         // and populate it with id and single_end values
         // for compatibility with nf-core module expected parameters
         // reads are just the list of fastq
         def meta        = [:]

         // 
         // TODO: should document how filenames are stripped to make sample IDs
         // 
         // this map gets rid of any of the following at the end of sample IDs:
         // .gz
         // .fastq
         // .fq
         // _001
         // _R[12]
         // _S\d+ 
         // E.g. strip _S1 from the end of a sample ID..  
         // This is typically sample #s from Illumina basecalling.
         // could cause an issue if sequenced the same sample with 
         // multiple barcodes so was repeated on a sample sheet. 
         meta.id         = name.replaceAll( /.gz$/ ,"")
         meta.id         = meta.id.replaceAll( /.fastq$/ ,"")
         meta.id         = meta.id.replaceAll( /.fq$/ ,"")
         meta.id         = meta.id.replaceAll( /.uniq$/ ,"")
         meta.id         = meta.id.replaceAll( /.trim$/ ,"")
         meta.id         = meta.id.replaceFirst( /_001$/ ,"")
         meta.id         = meta.id.replaceFirst( /_R[12]$/ ,"")
         meta.id         = meta.id.replaceFirst( /_S\d+$/ ,"")

         // if 2 fastq files then paired end data, so single_end is false
         meta.single_end = reads[1] ? false : true

         // this last statement in the map closure is the return value
         [ meta, reads ] }

  .set { ch_reads }

  // double check input is compressed - will error if not
  CHECK_FASTQ_COMPRESSED ( ch_reads ) 

  // optionally subsample input files to a certain # of random reads
  // (sampling without replacement)
  // COULD_DO: could count fastq first and only subsample if #reads > sample size
  if (subsample_size) {
     SEQTK_SAMPLE(ch_reads.map{ meta, reads -> [meta, reads, subsample_size] })
	  ch_reads = SEQTK_SAMPLE.out.reads
	  ch_versions = ch_versions.mix(SEQTK_SAMPLE.out.versions)
  }

  // count # of reads in each fastq file (or file pair)
  if (count_fastq) {
    COUNT_FASTQ ( ch_reads.map{ meta, reads -> [ meta, reads, "post_trimming"] } )
    ch_fastq_counts = ch_fastq_counts.mix(COUNT_FASTQ.out.count_file)
  }
  
 emit: 
  reads           = ch_reads
  versions        = ch_versions
  fastq_counts    = ch_fastq_counts

}
