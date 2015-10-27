# vcf_annotation_pipeline

Author: Khalid Mahmood (khalid.mahmood@unimelb.edu.au)

vcf_annotation_pipeline is a software pipeline that take as input a raw vcf file and performs a series of operations to filtering and annotate the variants. The tasks are:

 1. Normalize the vcf file.
 2. Decompose vcf file.
 3. Filter common variants e.g. 1000 genome CEU with maf < 1%.
 4. Annotate using VEP.
 5. Annotate using SnpEff.
 6. Annotate using [VarPub](https://github.com/khalidm/VarPub)
 7. Remove

# Tools requirements

`vcf_annotation_pipeline` requires the following tools.

# Data requires

'vcf_annotation_pipeline' requires the following datasets.

## License

3 Clause BSD License. See LICENSE.txt in source repository.
