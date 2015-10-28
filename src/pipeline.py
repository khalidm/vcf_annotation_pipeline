'''
Build the pipeline workflow by plumbing the stages together.
'''

from ruffus import Pipeline, suffix, formatter, add_inputs, output_from
from stages import Stages


def make_pipeline(state):
    '''Build the pipeline by constructing stages and connecting them together'''
    # Build an empty pipeline
    pipeline = Pipeline(name='vcf_annotation')
    # Get a list of paths to all the FASTQ files
    vcf_files = state.config.get_option('vcfs')
    # Stages are dependent on the state
    stages = Stages(state)

    # The original VCF files
    # This is a dummy stage. It is useful because it makes a node in the
    # pipeline graph, and gives the pipeline an obvious starting point.
    pipeline.originate(
        task_func=stages.original_vcf,
        name='original_vcf',
        output=vcf_file)

    # Decompose VCF using Vt
    pipeline.transform(
        task_func=stages.decompose_vcf,
        name='decompose_vcf',
        input=output_from('original_vcf'),
        # This will be the first input to the stage.
        # We assume the sample name may consist of only alphanumeric
        # characters.
        filter=formatter('.+/(?P<sample>[a-zA-Z0-9]+).vcf'),
        # Add an "extra" argument to the state (beyond the inputs and outputs)
        # which is the VCF file name (e.g. study/family name.
        # This is needed within the stage for finding out sample specific
        # configuration options
        extras=['{sample[0]}'],
        # The output file name is the sample name with a .bam extension.
        output='{path[0]}/{sample[0]}.decompose.normalize.vcf')

    # FILTER COMMON VARIANTS
    # ADD FILTER COMMON VARIANTS USING VEP

    # Annotate using VEP
    pipeline.transform(
        task_func=stages.annotate_vep,
        name='annotate_vep',
        input=output_from('decompose_vcf'),
        filter=suffix('.vcf'),
        output='.vep.vcf')

    # Annotate using SnpEff
    pipeline.transform(
        task_func=stages.annotate_snpeff,
        name='annotate_snpeff',
        input=output_from('annotate_vep'),
        filter=suffix('.vcf'),
        output='.snpeff.vcf')

    # Mark duplicates in the BAM file using Picard
    pipeline.transform(
        task_func=stages.mark_duplicates_picard,
        name='mark_duplicates_picard',
        input=output_from('sort_bam_picard'),
        filter=suffix('.sort.bam'),
        # XXX should make metricsup an extra output?
        output=['.sort.dedup.bam', '.metricsdup'])

    # Generate chromosome intervals using GATK
    pipeline.transform(
        task_func=stages.chrom_intervals_gatk,
        name='chrom_intervals_gatk',
        input=output_from('mark_duplicates_picard'),
        filter=suffix('.sort.dedup.bam'),
        output='.chr.intervals')

    # Local realignment using GATK
    (pipeline.transform(
        task_func=stages.local_realignment_gatk,
        name='local_realignment_gatk',
        input=output_from('chrom_intervals_gatk'),
        filter=formatter('.+/(?P<sample>[a-zA-Z0-9]+).chr.intervals'),
        add_inputs=add_inputs('{path[0]}/{sample[0]}.sort.dedup.bam'),
        output='{path[0]}/{sample[0]}.sort.dedup.realn.bam')
        .follows('mark_duplicates_picard'))

    # Base recalibration using GATK
    pipeline.transform(
        task_func=stages.base_recalibration_gatk,
        name='base_recalibration_gatk',
        input=output_from('local_realignment_gatk'),
        filter=suffix('.sort.dedup.realn.bam'),
        output=['.recal_data.csv', '.count_cov.log'])

    # Print reads using GATK
    (pipeline.transform(
        task_func=stages.print_reads_gatk,
        name='print_reads_gatk',
        input=output_from('base_recalibration_gatk'),
        filter=formatter('.+/(?P<sample>[a-zA-Z0-9]+).recal_data.csv'),
        add_inputs=add_inputs('{path[0]}/{sample[0]}.sort.dedup.realn.bam'),
        output='{path[0]}/{sample[0]}.sort.dedup.realn.recal.bam')
        .follows('local_realignment_gatk'))

    # Call variants using GATK
    pipeline.transform(
        task_func=stages.call_variants_gatk,
        name='call_variants_gatk',
        input=output_from('print_reads_gatk'),
        filter=suffix('.sort.dedup.realn.recal.bam'),
        output='.raw.snps.indels.g.vcf')

    # Combine G.VCF files for all samples using GATK
    pipeline.merge(
        task_func=stages.combine_gvcf_gatk,
        name='combine_gvcf_gatk',
        input=output_from('call_variants_gatk'),
        output='PCExomes.mergegvcf.vcf')

    # Genotype G.VCF files using GATK
    pipeline.transform(
        task_func=stages.genotype_gvcf_gatk,
        name='genotype_gvcf_gatk',
        input=output_from('combine_gvcf_gatk'),
        filter=suffix('.mergegvcf.vcf'),
        output='.genotyped.vcf')

    # SNP recalibration using GATK
    pipeline.transform(
        task_func=stages.snp_recalibrate_gatk,
        name='snp_recalibrate_gatk',
        input=output_from('genotype_gvcf_gatk'),
        filter=suffix('.genotyped.vcf'),
        output=['.snp_recal', '.snp_tranches', '.snp_plots.R'])

    # INDEL recalibration using GATK
    pipeline.transform(
        task_func=stages.indel_recalibrate_gatk,
        name='indel_recalibrate_gatk',
        input=output_from('genotype_gvcf_gatk'),
        filter=suffix('.genotyped.vcf'),
        output=['.indel_recal', '.indel_tranches', '.indel_plots.R'])

    # Apply SNP recalibration using GATK
    (pipeline.transform(
        task_func=stages.apply_snp_recalibrate_gatk,
        name='apply_snp_recalibrate_gatk',
        input=output_from('genotype_gvcf_gatk'),
        filter=suffix('.genotyped.vcf'),
        add_inputs=add_inputs(['PCExomes.snp_recal', 'PCExomes.snp_tranches']),
        output='.recal_SNP.vcf')
        .follows('snp_recalibrate_gatk'))

    # Apply INDEL recalibration using GATK
    (pipeline.transform(
        task_func=stages.apply_indel_recalibrate_gatk,
        name='apply_indel_recalibrate_gatk',
        input=output_from('genotype_gvcf_gatk'),
        filter=suffix('.genotyped.vcf'),
        add_inputs=add_inputs(
            ['PCExomes.indel_recal', 'PCExomes.indel_tranches']),
        output='.recal_INDEL.vcf')
        .follows('indel_recalibrate_gatk'))

    # Combine variants using GATK
    (pipeline.transform(
        task_func=stages.combine_variants_gatk,
        name='combine_variants_gatk',
        input=output_from('apply_snp_recalibrate_gatk'),
        filter=suffix('.recal_SNP.vcf'),
        add_inputs=add_inputs(['PCExomes.recal_INDEL.vcf']),
        output='.combined.vcf')
        .follows('apply_indel_recalibrate_gatk'))

    # Select variants using GATK
    pipeline.transform(
        task_func=stages.select_variants_gatk,
        name='select_variants_gatk',
        input=output_from('combine_variants_gatk'),
        filter=suffix('.combined.vcf'),
        output='.selected.vcf')

    return pipeline
