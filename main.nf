/*
This bash script provides a means of generating scaled strand-specific BIGWIG files
from a BAM file containing paired reads.
from https://github.com/crickbabs/DRB_TT-seq/blob/master/bigwig.md
*/
nextflow.enable.dsl=2

params.scale_factor=1
params.bams="*.bam"

log.info """\
         bigwig - N F   P I P E L I N E
         ===================================
         bams           : ${params.bams}
         scale_factor   : ${params.scale_factor}
         """
         .stripIndent()



include { bigwig_all; bigwig_forward;bigwig_reverse } from './modules/bigwig.nf'
channel.fromPath(params.bams).set{ bam_ch}

workflow  {

  bam_ch.view()
  bigwig_all(bam_ch)
  bigwig_forward(bam_ch)
  bigwig_reverse(bam_ch)

}
