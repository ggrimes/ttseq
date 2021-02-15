/*
This bash script provides a means of generating strand-specific metagene, TSS and TES profiles from a BAM file using "ngs.plot".
Of course, this will only work if your libraries were created in a strand-specific fashion.
*/


params.genome = "hg38"
params.L = 5000
params.bams="*.bam"



log.info """\
         ngplot - N F   P I P E L I N E
         ===================================
         genome           : ${params.genome}
         bams             : ${params.bams}
         L                : ${params.L}
         """
         .stripIndent()

regions_ch = Channel.fromList(["genebody","tss","tes"]).
strands_ch = Channel.fromList(["both","same", "opposite"])
Channel
         .fromFilePairs(params.bams) {file -> file.name.replaceAll(/.bam|.bai$/,'')}
          .into{ bam_ch; bam_rev_ch; bam_for_ch}

//Restrict to first mate reads.
process first_reads {

  input:
  tuple(val(sampleID),path(bam)) from bam_for_ch

  output:

  script:
  """
  samtools view --threads ${task.cpus} -h -b -f 64 ${bam} -o $MATE1
  samtools index $MATE1

  MATE1REHEADER="${WORKDIR}WT.mate1.reheader.bam"

  samtools view --threads ${task.cpus} -H ${MATE1} |\
    sed -e 's/SN:\([0-9XY]*\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' |\
    samtools reheader - ${MATE1} > ${MATE1REHEADER}
  samtools index ${MATE1REHEADER}
  """
}

/*
Run ngs.plot.
Finally, use ngs.plot to create sense and anti-sense profiles for your regions of interest
using the correct BAM file.
Note that if your libraries were produced such that mate2 represents the forward strand the sense/anti-sense profiles will be reversed.
*/

process ngsplot {

  tag "${region} ${strand}"
  conda  "$baseDir/environment.yml"
  publishDir "results/ngsplot" , mode: 'copy'

  input:
  val(region) from regions_ch
  each(strand) from strands_ch
  tuple(val(sampleID),path(bam)) from bam_for_ch

  output:


  script:
  """
  ngs.plot.r -G ${params.genome} \
  -R ${region} \
  -C ${MATE1REHEADER} \
  -O ${sampleID} \
  -P ${task.cpus} \
  -SS ${strand} \
  -SE 1 \
  -L ${params.L} \
  -F chipseq \
  -D ensembl
  """
}
