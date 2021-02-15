//Create bigwig file for all reads.
process bigwig_all {
 label "bigwig"
 tag "${sampleID} bigwig_all"
 conda  "$baseDir/environment.yml"
 publishDir "results/bigwig" , mode: 'copy'


 input:
 tuple(val(sampleID),path(bam))

 output:
 path("${sampleID}.bigwig")

 script:
 """
 samtools index ${bam}
 bamCoverage --scaleFactor ${params.scale_factor} \
 -p ${task.cpus}  \
 -b ${bam} \
 -o ${sampleID}.bigwig
 """
}

/*
Get file for transcripts originating on the forward strand.
Include reads that are 2nd in a pair (128).
Exclude reads that are mapped to the reverse strand (16)

Exclude reads that are mapped to the reverse strand (16) and
 first in a pair (64): 64 + 16 = 80
*/
//Create bigwig file for all reads.
process bigwig_forward {
 label "bigwig"
 tag "${sampleID} bigwig_all"
 conda  "$baseDir/environment.yml"
 publishDir "results/bigwig" , mode: 'copy'


 input:
 tuple(val(sampleID),path(bam))

 output:
 path("${sampleID}_forward.bigwig")

 script:
 """
 samtools index ${bam}
 samtools view -b -f 128 -F 16 --threads ${task.cpus} ${bam} > ${sampleID}"_FOR1.bam"
 samtools view -b -f 80  --threads ${task.cpus} ${bam} > ${sampleID}"_FOR2.bam"
 samtools merge --threads ${task.cpus} -f ${sampleID}"_FOR.bam" ${sampleID}"_FOR1.bam" ${sampleID}"_FOR2.bam"
 samtools index ${sampleID}"_FOR.bam"
 bamCoverage --scaleFactor ${params.scale_factor} -p ${task.cpus} -b ${sampleID}"_FOR.bam" -o ${sampleID}"_forward.bigwig"
 """
}




/*
Get the file for transcripts that originated from the reverse strand:
Include reads that map to the reverse strand (128) and are second in a pair (16): 128 + 16 = 144
Include reads that are first in a pair (64), but exclude those ones that map to the reverse strand (16)
*/
//Create bigwig file for all reads.
process bigwig_reverse {
 label "bigwig"
 tag "${sampleID} bigwig_all"
 conda  "$baseDir/environment.yml"
 publishDir "results/bigwig" , mode: 'copy'


 input:
 tuple(val(sampleID),path(bam))

 output:
 path("${sampleID}_reverse.bigwig")

 script:
 """
 samtools index ${bam}
 samtools view -b -f 144 --threads ${task.cpus} ${bam} > ${sampleID}"_REV1.bam"
 samtools view -b -f 64 -F 16 --threads ${task.cpus} ${bam} > ${sampleID}"_REV2.bam"
 samtools merge --threads ${task.cpus} -f ${sampleID}"_REV.bam" ${sampleID}"_REV1.bam" ${sampleID}"_REV2.bam"
 samtools index ${sampleID}"_REV.bam"
 bamCoverage --scaleFactor ${params.scale_factor} -p ${task.cpus} -b ${sampleID}"_REV.bam" -o ${sampleID}"_reverse.bigwig"
 """
}
