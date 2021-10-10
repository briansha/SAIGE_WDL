version 1.0

## Version 10-10-2021
##
## This WDL workflow runs SAIGE. - https://github.com/weizhouUMICH/SAIGE/wiki/Genetic-association-tests-using-SAIGE
##
## Special Notes:
## - In Step 1 - the genotype file (plink file) needs to be subset (using plink2) to include only the sample IDs for those contained in the .vcf.gz files for Step 2.
##             - From the documentation: "Currently SAIGE requires all samples used in glmm model are in the dosage file"
##             - Check to make sure all sample IDs in Step 1 are included in the VCF files in Step 2.
## - Step 2 - currently set to use 22(+) .vcf.gz files relating to individual chromosomal data.
## - Step 2 - The dosage file (vcf, bgen, or sav) requires .csi index files instead of .tbi index files for indexing.
##
## Other WDLs of mine that go with this:
## - vcf_split: Splits vcf.gz files into smaller pieces for Step 2 so it can execute faster.
##
## Cromwell version support - Successfully tested on v67
##
## Distributed under terms of the MIT License
## Copyright (c) 2021 Brian Sharber
## Contact <brian.sharber@vumc.org>

workflow run_saige {
    input {
        File vcf_and_chrom_file_step2                        # File containing tab-separated list of: .vcf.gz file, corresponding .vcf.gz.csi index file, and corresponding chromosome number - for all files used for Step2.
        String saige_docker = "wzhou88/saige:0.44.5"
        String ubuntu_docker = "ubuntu:20.10"
        String r_base_docker = "briansha/saige_r_base:4.1.0" # Docker image built off of Ubuntu 18 base image with R 4.1.0 installed.
    }
    Array[Array[String]] data_step2 = read_tsv(vcf_and_chrom_file_step2)

    call saige_step1_fitNULL  {
        input:
          docker = saige_docker
    }

    scatter (file in data_step2) {
      call saige_step2_SPAtests {
        input:
          vcfFile = file[0],
          vcfFileIndex = file[1],
          chrom = file[2],
          GMMATmodelFile = saige_step1_fitNULL.GMMATmodelFile,
          varianceRatioFile = saige_step1_fitNULL.varianceRatioFile,
          docker = saige_docker
      }
    }

    call combine_saige_results {
      input:
        saige_result_files = saige_step2_SPAtests.saige_output_file,
        docker = ubuntu_docker
    }

    call Plots{
      input:
        merged_saige_file = combine_saige_results.merged_saige_file,
        docker = r_base_docker
    }

    output {
        File merged_saige_file = combine_saige_results.merged_saige_file
        Array[File] output_plots = Plots.output_plots
        Array[File] output_saige = Plots.output_saige
    }

    meta {
        author : "Brian Sharber"
        email : "brian.sharber@vumc.org"
        description : "Run SAIGE"
    }
}

task saige_step1_fitNULL {
    input {
        # Genotype file
        File bed_step1     # Step 1 requires a .bed, .bim, .fam file - Path to plink file for creating the genetic relationship matrix (GRM). minMAFforGRM can be used to specify the minimum MAF of markers in he plink file to be used for constructing GRM. Genetic markers are also randomly selected from the plink file to estimate the variance ratios.
        File bim_step1
        File fam_step1

        # Phenotype file
        File phenoFile     # Should include phenotypes, covariates, and samples - Path to the phenotype file. The file can be either tab or space delimited. The phenotype file has a header and contains at least two columns. One column is for phentoype and the other column is for sample IDs. Additional columns can be included in the phenotype file for covariates in the null GLMM. Please note that covariates to be used in the NULL GLMM need to specified using the argument covarColList

        # Other parameters
        File? sparseGRMFile                                # Path to the pre-calculated sparse GRM file. If not specified and  IsSparseKin=TRUE, sparse GRM will be computed [default=NULL]
        File? sparseGRMSampleIDFile                        # Path to the sample ID file for the pre-calculated sparse GRM. No header is included. The order of sample IDs is corresponding to sample IDs in the sparse GRM [default=NULL]
        String phenoCol                                    # Column name for phenotype to be tested in the phenotype file, e.g CAD - Turn into list of phenotypes later
        String covarColList                                # Covariate headers in the phenoFile - List of covariates (comma separated)
        String sampleIDColinphenoFile                      # Sample ID headers in the phenoFile - Column name of sample IDs in the phenotype file, e.g. IID
        String traitType = "binary"                        # binary/quantitative [default=binary]
        String invNormalize = "FALSE"                      # if quantitative, asks SAIGE whether to perform the inverse normalization for the phenotype [default='FALSE']
        String LOCO = "FALSE"                              # Whether to apply the leave-one-chromosome-out (LOCO) approach. This option has not been extensively tested [default=FALSE]. TRUE works in version 0.36.5 and up
        String IsSparseKin = "FALSE"                       # Whether to use sparse kinship for association test [default='FALSE'] - TRUE must be specified for SAIGE-GENE
        String isCateVarianceRatio = "FALSE"               # Whether to estimate variance ratio based on different MAC categories. If yes, variance ratio will be estiamted for multiple MAC categories corresponding to cateVarRatioMinMACVecExclude and cateVarRatioMaxMACVecInclude. Currently, if isCateVarianceRatio=TRUE, then LOCO=FALSE [default=FALSE]
        String tauInit = "0,0"                             # Initial values for tau. [default=0,0] - see quantitative traits section in QUILT documentation.
        String? minMAFforGRM                               # Minimum MAF of markers used for GRM
        Int minCovariateCount = 1                          # If binary covariates have a count less than this, they will be excluded from the model to avoid convergence issues [default=-1] (no covariates will be excluded).
        Int numRandomMarkerforSparseKin = 2000             # Number of randomly selected markers to be used to identify related samples for sparse GRM [default=2000]
        Int numRandomMarkerforVarianceRatio = 30           # An integer greater than 0. Number of markers to be randomly selected for estimating the variance ratio. The number will be automatically added by 10 untial the coefficient of variantion (CV) for the variance ratio estimate is below ratioCVcutoff [default=30].
        Float relatednessCutoff = 0.125                    # For SAIGE-GENE - Threshold to treat two samples as unrelated if IsSparseKin is TRUE [default=0.125]
        Float tol = 0.02                                   # Tolerance for fitting the null GLMM to converge [default=0.02].
        Int maxiter = 20                                   # Maximum number of iterations used to fit the null GLMM [default=20].
        String? tolPCG                                     # Tolerance for PCG to converge [default=1e-5].
        Int maxiterPCG = 500                               # Maximum number of iterations for PCG [default=500].
        Int SPAcutoff = 2                                  # Cutoff for the deviation of score test statistics from mean in the unit of sd to perform SPA [default=2].
        String skipModelFitting = "FALSE"                  # Whether to skip model fitting and only to estimate the variance ratio. If TRUE, the file outputPrefix.rda is required [default='FALSE']
        Int memoryChunk = 2                                # Size (Gb) for each memory chunk [default=2]
        Float traceCVcutoff = 0.0025                       # Threshold for coefficient of variation (CV) for the trace estimator. Number of runs for trace estimation will be increased until the CV is below the threshold [default=0.0025].
        Float ratioCVcutoff = 0.001                        # Threshold for coefficient of variation (CV) for estimating the variance ratio. The number of randomly selected markers will be increased until the CV is below the threshold [default=0.001]
        String IsOverwriteVarianceRatioFile = "FALSE"      # Whether to overwrite the variance ratio file if the file exist.[default='FALSE']
        String isCovariateTransform = "TRUE"               # Whether use qr transformation on non-genetic covariates [default='TRUE'].
        String isDiagofKinSetAsOne = "FALSE"               # Whether to set the diagnal elements in GRM to be 1 [default='FALSE'].
        String? useSparseSigmaConditionerforPCG            # Whether to sparse GRM to speed up the PCG. Current this option is deactivated. [default='FALSE'].
        String useSparseSigmaforInitTau = "FALSE"          # Whether to use sparse Sigma to estiamte initial tau [default='FALSE'].
        String useSparseGRMtoFitNULL = "FALSE"             # Whether to use sparse GRM to fit the null model [default='FALSE'].
        String includeNonautoMarkersforVarRatio = "FALSE"  # Whether to allow for non-autosomal markers for variance ratio. [default, 'FALSE']
        String FemaleOnly = "FALSE"                        # Whether to run Step 1 for females only [default=FALSE]. if TRUE, --sexCol and --FemaleCode need to be specified
        String MaleOnly = "FALSE"                          # Whether to run Step 1 for males only [default=FALSE]. if TRUE, --sexCol and --MaleCode need to be specified
        String? sexCol                                     # Column name for sex in the phenotype file, e.g Sex
        String? FemaleCode                                 # Values in the column for sex in the phenotype file are used for females [default, '1']
        String? MaleCode                                   # Values in the column for sex in the phenotype file are used for males [default, '0']
        String? noEstFixedEff = "FALSE"                    # Whether to estimate fixed effect coefficients. [default, 'FALSE']

        String cateVarRatioMinMACVecExclude = "0.5,1.5,2.5,3.5,4.5,5.5,10.5,20.5"  # vector of float. Lower bound of MAC for MAC categories. The length equals to the number of MAC categories for variance ratio estimation. [default='0.5,1.5,2.5,3.5,4.5,5.5,10.5,20.5']
        String cateVarRatioMaxMACVecInclude = "1.5,2.5,3.5,4.5,5.5,10.5,20.5"      # vector of float. Higher bound of MAC for MAC categories. The length equals to the number of MAC categories for variance ratio estimation minus 1. [default='1.5,2.5,3.5,4.5,5.5,10.5,20.5']
        #String outputPrefix # Add this in later to allow for multiple phenotypes in the analysis.
        Boolean help = false  # Show this help message and exit

        # Runtime
        String docker
        Float memory = 16.0
        Int? disk_size_override
        Int cpu = 4
        Int preemptible = 1
    	Int maxRetries = 0
    }
    Float bed_size = size(bed_step1, "GiB")
    Int disk = select_first([disk_size_override, ceil(10.0 + 2.0 * bed_size)])

    command <<<
        set -euo pipefail
        step1_fitNULLGLMM.R     \
        ~{if defined(bed_step1) then "--plinkFile=~{sub(bed_step1,'\\.bed$','')} " else " "} \
        ~{if defined(phenoFile) then "--phenoFile=~{phenoFile} " else " "} \
        ~{if defined(sparseGRMFile) then "--sparseGRMFile=~{sparseGRMFile} " else " "} \
        ~{if defined(sparseGRMSampleIDFile) then "--sparseGRMSampleIDFile=~{sparseGRMSampleIDFile} " else " "} \
        ~{if defined(phenoCol) then "--phenoCol=~{phenoCol} " else " "} \
        ~{if defined(covarColList) then "--covarColList=~{covarColList} " else " "} \
        ~{if defined(sampleIDColinphenoFile) then "--sampleIDColinphenoFile=~{sep=', ' sampleIDColinphenoFile} " else " "} \
        ~{if defined(traitType) then "--traitType=~{traitType} " else " "} \
        ~{if defined(invNormalize) then "--invNormalize=~{invNormalize} " else " "} \
        ~{if defined(LOCO) then "--LOCO=~{LOCO} " else " "} \
        ~{if defined(IsSparseKin) then "--IsSparseKin=~{IsSparseKin} " else " "} \
        ~{if defined(isCateVarianceRatio) then "--isCateVarianceRatio=~{isCateVarianceRatio} " else " "} \
        ~{if defined(tauInit) then "--tauInit=~{tauInit} " else " "} \
        ~{if defined(minMAFforGRM) then "--minMAFforGRM=~{minMAFforGRM} " else " "} \
        ~{if defined(minCovariateCount) then "--minCovariateCount=~{minCovariateCount} " else " "} \
        ~{if defined(numRandomMarkerforSparseKin) then "--numRandomMarkerforSparseKin=~{numRandomMarkerforSparseKin} " else " "} \
        ~{if defined(numRandomMarkerforVarianceRatio) then "--numRandomMarkerforVarianceRatio=~{numRandomMarkerforVarianceRatio} " else " "} \
        ~{if defined(relatednessCutoff) then "--relatednessCutoff=~{relatednessCutoff} " else " "} \
        ~{if defined(tol) then "--tol=~{tol} " else " "} \
        ~{if defined(maxiter) then "--maxiter=~{maxiter} " else " "} \
        ~{if defined(tolPCG) then "--tolPCG=~{tolPCG} " else " "} \
        ~{if defined(maxiterPCG) then "--maxiterPCG=~{maxiterPCG} " else " "} \
        ~{if defined(SPAcutoff) then "--SPAcutoff=~{SPAcutoff} " else " "} \
        ~{if defined(skipModelFitting) then "--skipModelFitting=~{skipModelFitting} " else " "} \
        ~{if defined(memoryChunk) then "--memoryChunk=~{memoryChunk} " else " "} \
        ~{if defined(traceCVcutoff) then "--traceCVcutoff=~{traceCVcutoff} " else " "} \
        ~{if defined(ratioCVcutoff) then "--ratioCVcutoff=~{ratioCVcutoff} " else " "} \
        ~{if defined(IsOverwriteVarianceRatioFile) then "--IsOverwriteVarianceRatioFile=~{IsOverwriteVarianceRatioFile} " else " "} \
        ~{if defined(isCovariateTransform) then "--isCovariateTransform=~{isCovariateTransform} " else " "} \
        ~{if defined(isDiagofKinSetAsOne) then "--isDiagofKinSetAsOne=~{isDiagofKinSetAsOne} " else " "} \
        ~{if defined(useSparseSigmaConditionerforPCG) then "--useSparseSigmaConditionerforPCG=~{useSparseSigmaConditionerforPCG} " else " "} \
        ~{if defined(useSparseSigmaforInitTau) then "--useSparseSigmaforInitTau=~{useSparseSigmaforInitTau} " else " "} \
        ~{if defined(useSparseGRMtoFitNULL) then "--useSparseGRMtoFitNULL=~{useSparseGRMtoFitNULL} " else " "} \
        ~{if defined(includeNonautoMarkersforVarRatio) then "--includeNonautoMarkersforVarRatio=~{includeNonautoMarkersforVarRatio} " else " "} \
        ~{if defined(FemaleOnly) then "--FemaleOnly=~{FemaleOnly} " else " "} \
        ~{if defined(MaleOnly) then "--MaleOnly=~{MaleOnly} " else " "} \
        ~{if defined(sexCol) then "--sexCol=~{sexCol} " else " "} \
        ~{if defined(FemaleCode) then "--FemaleCode=~{FemaleCode} " else " "} \
        ~{if defined(MaleCode) then "--MaleCode=~{MaleCode} " else " "} \
        ~{if defined(noEstFixedEff) then "--noEstFixedEff=~{noEstFixedEff} " else " "} \
        ~{if defined(cateVarRatioMinMACVecExclude) then "--cateVarRatioMinMACVecExclude=~{cateVarRatioMinMACVecExclude} " else " "} \
        ~{if defined(cateVarRatioMaxMACVecInclude) then "--cateVarRatioMaxMACVecInclude=~{cateVarRatioMaxMACVecInclude} " else " "} \
        ~{if help then "--help " else " "} \
        --nThreads=~{cpu} \
        --outputPrefix=saige_step1_~{phenoCol}
	>>>

    runtime {
	  docker: docker
          memory: memory + " GiB"
	  disks: "local-disk " + disk + " HDD"
          cpu: cpu
	  preemptible: preemptible
          maxRetries: maxRetries
	}

    output {
	  File GMMATmodelFile = "saige_step1_" + phenoCol + ".rda"
          File varianceRatioFile = "saige_step1_" + phenoCol + ".varianceRatio.txt"
          Array[File] other_files = glob("*.txt")
    }
}

task saige_step2_SPAtests {
    input {
        # Dosage file
        File? vcfFile                # Path to the vcf file. Supports VCF and BCF files.
        File? bgenFile               # Path to bgen file. Currently version 1.2 with 8 bit compression is supported
        File? savFile                # Path to the sav file.

        # Index file
        File? vcfFileIndex           # Path to vcf index file. Indexed by tabix. Path to index for vcf file by tabix, .tbi file by tabix -p vcf file.vcf.gz
        File? bgenFileIndex          # Path to the .bgi file (index of the bgen file)
        File? savFileIndex           # Path to the .s1r file (index of the sav file).

        # Other parameters
        File? rangestoExcludeFile                                      # Path to a file containing genome regions to be excluded from the bgen file. The file contains three columns for chromosome, start, and end respectively with no header.
        File? rangestoIncludeFile                                      # Path to a file containing genome regions to be included from the bgen file. The file contains three columns for chromosome, start, and end respectively with no header.
        File? sampleFile                                               # Imputed samples file. Path to the file that contains one column for IDs of samples in the dosage file. For version >= 0.38, this file is only needed for bgen files.
        File? groupFile                                                # Path to the file containing the group information for gene-based tests. Each line is for one gene/set of variants. The first element is for gene/set name. The rest of the line is for variant ids included in this gene/set. For vcf/sav, the genetic marker ids are in the format chr:pos_ref/alt. For bgen, the genetic marker ids should match the ids in the bgen file. Each element in the line is seperated by tab.
        File? GMMATmodelFile                                           # Path to the input file containing the glmm model, which is output from previous step. Will be used by load()
        File? varianceRatioFile                                        # Path to the input file containing the variance ratio, which is output from the previous step
        File? sparseSigmaFile                                          # Path to the file containing the sparse Sigma output by step 1. The suffix of this file is .mtx
        File? idstoExcludeFile                                         # Path to a file containing variant ids to be excluded from the bgen file. The file does not have a header and each line is for a marker ID.
        File? idstoIncludeFile                                         # Path to a file containing variant ids to be included from the bgen file. The file does not have a header and each line is for a marker ID.
        File? sampleFile_male                                          # Path to the file containing one column for IDs of MALE samples in the bgen or vcf file with NO header. Order does not matter
        String? IsDropMissingDosages                                   # Whether to drop samples with missing dosages. If FALSE, the missing dosages with be mean imputed, otherwise, they will be removed before testing. This option only works for bgen, vcf, and sav input.
        String IsAccountforCasecontrolImbalanceinGroupTest = "TRUE"    # Whether to account for unbalanced case-control ratios for binary tratis in gene- or region-based tests. [default=TRUE]
        String IsSparse = "TRUE"                                       # Whether to exploit the sparsity of the genotype vector for less frequent variants to speed up the SPA tests or not for binary traits [default=TRUE].
        Int start = 1                                                  # start genome position in the vcf to be tested [default=1]
        Int end = 250000000                                            # end genome position in the vcf to be tested. If not specified, the whole genome will be tested [default=250000000]
        Float SPAcutoff = 2.0                                          # If the test statistic lies within the standard deviation cutoff of the mean, p-value based on traditional score test is returned. Default value is 2.
        String vcfField = "DS"                                         # "GT" - VCF containing genotypes. "DS" - VCF containing dosages or SAV for dosages.
        String? chrom                                                  # Chromosome in vcf to be tested. The string needs to exactly match the chromosome string in the vcf/sav file. For example, '1' does not match 'chr1'. If not specified, all markers in the vcf will be tested. If LOCO is specified, providing chrom will save computation cost.
        Int numLinesOutput = 10000                                     # Number of  markers to be output each time [default=10000]

        String IsOutputAFinCaseCtrl = "FALSE"                                     # For binary traits - TRUE can be specified to output allele frequency in cases and controls for dichotomous traits [default=FALSE]
        String IsOutputNinCaseCtrl = "FALSE"                                      # For binary traits - TRUE can be specified to output sample sizes in cases and controls for dichotomous traits [default=FALSE]
        String IsOutputHetHomCountsinCaseCtrl = "FALSE"                           # For binary traits - can be specified to output heterozygous and homozygous counts in cases and controls. By default, FALSE. If True, the columns homN_Allele2_cases, hetN_Allele2_cases, homN_Allele2_ctrls, hetN_Allele2_ctrls will be output [default=FALSE]
        String? method_to_CollapseUltraRare                                       # For binary traits - Method to collpase the ultra rare variants in the set-based association tests for BINARY traits only. This argument can be 'absence_or_presence', 'sum_geno', or ''. absence_or_presence:
                                                                                  #   For the resulted collpased marker, any individual having DosageCutoff_for_UltraRarePresence <= dosage < 1+DosageCutoff_for_UltraRarePresence for any ultra rare variant has 1 in the genotype vector, having dosage >= 1+DosageCutoff_for_UltraRarePresence for any ultra rare variant has 2 in the genotype vector, otherwise 0.
                                                                                  #   sum_geno: Ultra rare variants with MAC <=  MACCutoff_to_CollapseUltraRare will be collpased for set-based tests in the 'sum_geno' way and the resulted collpased marker's genotype equals weighted sum of the genotypes of all ultra rare variants. NOTE: this option sum_geno currently is NOT active. By default, ''
        String IsSingleVarinGroupTest = "FALSE"                                   # Whether to perform single-variant assoc tests for genetic markers included in the gene-based tests. By default, FALSE
        String IsOutputPvalueNAinGroupTestforBinary = "FALSE"                     # Whether to output p value if not account for case-control imbalance when performing group test (only for binary traits). [default=FALSE]
        String weightsIncludeinGroupFile = "FALSE"                                # Whether to specify customized weight for makers in gene- or region-based tests. If TRUE, weights are included in the group file. For vcf/sav, the genetic marker ids and weights are in the format chr:pos_ref/alt;weight. For bgen, the genetic marker ids should match the ids in the bgen filE, e.g. SNPID;weight. Each element in the line is seperated by tab. [default=FALSE]
        String LOCO = "FALSE"                                                     # Whether to apply the leave-one-chromosome-out option. This option has not been extensively tested.
        String IsOutputMAFinCaseCtrlinGroupTest = "FALSE"                         # Whether to output minor allele frequency in cases and controls in set-based tests By default, FALSE
        String IsOutputBETASEinBurdenTest = "FALSE"                               # Whether to output effect sizes for burden tests. [default=FALSE]
        String cateVarRatioMinMACVecExclude = "0.5,1.5,2.5,3.5,4.5,5.5,10.5,20.5" # vector of float. Lower bound of MAC for MAC categories. The length equals to the number of MAC categories for variance ratio estimation. [default='0.5,1.5,2.5,3.5,4.5,5.5,10.5,20.5']
        String cateVarRatioMaxMACVecInclude = "1.5,2.5,3.5,4.5,5.5,10.5,20.5"     # vector of float. Higher bound of MAC for MAC categories. The length equals to the number of MAC categories for variance ratio estimation minus 1. [default='1.5,2.5,3.5,4.5,5.5,10.5,20.5']
        Float dosageZerodCutoff = 0.2                                             # In gene- or region-based tests, for each variants with MAC <= 10, dosages <= dosageZerodCutoff with be set to 0. [default=0.2]
        Float MACCutoff_to_CollapseUltraRare = 10.0                               # MAC cutoff to collpase the ultra rare variants (<= MACCutoff_to_CollapseUltraRare) in the set-based association tests. By default, 10.
        Float DosageCutoff_for_UltraRarePresence = 0.5                            # Dosage cutoff to determine whether the ultra rare variants are absent or present in the samples. Dosage >= DosageCutoff_for_UltraRarePresence indicates the varaint in present in the sample. 0< DosageCutoff_for_UltraRarePresence <= 2. By default, 0.5
        String? X_PARregion                                                       # ranges of (pseudoautosomal) PAR region on chromosome X, which are seperated by comma and in the format start:end. By default: '60001-2699520,154931044-155260560' in the UCSC build hg19. For males, there are two X alleles in the PAR region, so PAR regions are treated the same as autosomes.
                                                                                  #   In the NON-PAR regions (outside the specified PAR regions on chromosome X), for males, there is only one X allele. If is_rewrite_XnonPAR_forMales=TRUE, genotypes/dosages of all variants in the NON-PAR regions on chromosome X will be mutliplied by 2.
        String is_rewrite_XnonPAR_forMales = "FALSE"                              # Whether to rewrite gentoypes or dosages of variants in the NON-PAR regions on chromosome X for males (multiply by 2). By default, FALSE. Note, only use is_rewrite_XnonPAR_forMales=TRUE when the specified VCF or Bgen file only has variants on chromosome X. When is_rewrite_XnonPAR_forMales=TRUE, the program does not check the chromosome value by assuming all variants are on chromosome X

        String? weights_for_G2_cond          # vector of float. weights for conditioning markers for gene- or region-based tests. The length equals to the number of conditioning markers, delimited by comma. e.g. '1,2,3'
        String? condition                    # For conditional analysis. Genetic marker ids (chr:pos_ref/alt if sav/vcf dosage input, marker id if bgen input) seperated by comma. e.g.chr3:101651171_C/T,chr3:101651186_G/A, Note that currently conditional analysis is only for bgen,vcf,sav input.
        Float minMAF = 0.0                   # Minimum minor allele frequency for markers to be tested. The higher threshold between minMAC and minMAF will be used [default=0].
        Float minMAC = 0.5                   # Minimum minor allele count for markers to be tested. The higher threshold between minMAC and minMAF will be used [default=0.5].
        Float maxMAFforGroupTest = 0.5       # Max MAF for markers tested in group test [default=0.5]
        Int minInfo = 0                      # Minimum Info for markers to be tested [default=0]
        Boolean help = false                 # Show this help message and exit

        # SKAT library
        String? method                       # Method for gene-based test p-values. Methods other than optimal.adj have not been extensively tested. More options can be seen in the SKAT library
        String? kernel                       # More options can be seen in the SKAT library
        Float? weights_beta_rare             # parameters for the beta distribution to weight genetic markers with MAF <= weightMAFcutoff in gene-based tests. More options can be seen in the SKAT library
        Float? weights_beta_common           # parameters for the beta distribution to weight genetic markers with MAF > weightMAFcutoff in gene-based tests. More options can be seen in the SKAT library. NOTE: this argument is not fully developed. currently, weights.beta.common is euqal to weights.beta.rare
        Float? weightMAFcutoff               # See document above for weights.beta.rare and weights.beta.common
        String? r_corr                       # More options can be seen in the SKAT library

        # Runtime
        String docker
	Float memory = 4.0
	Int? disk_size_override
        Int cpu = 1
        Int preemptible = 1
        Int maxRetries = 0
    }
    Float vcf_size = size(vcfFile, "GiB")
    Int disk = select_first([disk_size_override, ceil(10.0 + 2.0 * vcf_size)])

    command <<<
      set -euo pipefail
      step2_SPAtests.R \
        ~{if defined(vcfFile) then "--vcfFile=~{vcfFile} " else " "} \
        ~{if defined(bgenFile) then "--bgenFile=~{bgenFile} " else " "} \
        ~{if defined(savFile) then "--savFile=~{savFile} " else " "} \
        ~{if defined(vcfFileIndex) then "--vcfFileIndex=~{vcfFileIndex} " else " "} \
        ~{if defined(bgenFileIndex) then "--bgenFileIndex=~{bgenFileIndex} " else " "} \
        ~{if defined(savFileIndex) then "--savFileIndex=~{savFileIndex} " else " "} \
        ~{if defined(rangestoExcludeFile) then "--rangestoExcludeFile=~{rangestoExcludeFile} " else " "} \
        ~{if defined(rangestoIncludeFile) then "--rangestoIncludeFile=~{rangestoIncludeFile} " else " "} \
        ~{if defined(sampleFile) then "--sampleFile=~{sampleFile} " else " "} \
        ~{if defined(groupFile) then "--groupFile=~{groupFile} " else " "} \
        ~{if defined(GMMATmodelFile) then "--GMMATmodelFile=~{GMMATmodelFile} " else " "} \
        ~{if defined(varianceRatioFile) then "--varianceRatioFile=~{varianceRatioFile} " else " "} \
        ~{if defined(sparseSigmaFile) then "--sparseSigmaFile=~{sparseSigmaFile} " else " "} \
        ~{if defined(idstoExcludeFile) then "--idstoExcludeFile=~{idstoExcludeFile} " else " "} \
        ~{if defined(idstoIncludeFile) then "--idstoIncludeFile=~{idstoIncludeFile} " else " "} \
        ~{if defined(sampleFile_male) then "--sampleFile_male=~{sampleFile_male} " else " "} \
        ~{if defined(IsDropMissingDosages) then "--IsDropMissingDosages=~{IsDropMissingDosages} " else " "} \
        ~{if defined(IsAccountforCasecontrolImbalanceinGroupTest) then "--IsAccountforCasecontrolImbalanceinGroupTest=~{IsAccountforCasecontrolImbalanceinGroupTest} " else " "} \
        ~{if defined(IsSparse) then "--IsSparse=~{IsSparse} " else " "} \
        ~{if defined(start) then "--start=~{start} " else " "} \
        ~{if defined(end) then "--end=~{end} " else " "} \
        ~{if defined(SPAcutoff) then "--SPAcutoff=~{SPAcutoff} " else " "} \
        ~{if defined(vcfField) then "--vcfField=~{vcfField} " else " "} \
        ~{if defined(chrom) then "--chrom=~{chrom} " else " "} \
        ~{if defined(numLinesOutput) then "--numLinesOutput=~{numLinesOutput} " else " "} \
        ~{if defined(IsOutputAFinCaseCtrl) then "--IsOutputAFinCaseCtrl=~{IsOutputAFinCaseCtrl} " else " "} \
        ~{if defined(IsOutputNinCaseCtrl) then "--IsOutputNinCaseCtrl=~{IsOutputNinCaseCtrl} " else " "} \
        ~{if defined(IsOutputHetHomCountsinCaseCtrl) then "--IsOutputHetHomCountsinCaseCtrl=~{IsOutputHetHomCountsinCaseCtrl} " else " "} \
        ~{if defined(method_to_CollapseUltraRare) then "--method_to_CollapseUltraRare=~{method_to_CollapseUltraRare} " else " "} \
        ~{if defined(IsSingleVarinGroupTest) then "--IsSingleVarinGroupTest=~{IsSingleVarinGroupTest} " else " "} \
        ~{if defined(IsOutputPvalueNAinGroupTestforBinary) then "--IsOutputPvalueNAinGroupTestforBinary=~{IsOutputPvalueNAinGroupTestforBinary} " else " "} \
        ~{if defined(weightsIncludeinGroupFile) then "--weightsIncludeinGroupFile=~{weightsIncludeinGroupFile} " else " "} \
        ~{if defined(IsDropMissingDosages) then "--IsDropMissingDosages=~{IsDropMissingDosages} " else " "} \
        ~{if defined(LOCO) then "--LOCO=~{LOCO} " else " "} \
        ~{if defined(IsOutputMAFinCaseCtrlinGroupTest) then "--IsOutputMAFinCaseCtrlinGroupTest=~{IsOutputMAFinCaseCtrlinGroupTest} " else " "} \
        ~{if defined(IsOutputBETASEinBurdenTest) then "--IsOutputBETASEinBurdenTest=~{IsOutputBETASEinBurdenTest} " else " "} \
        ~{if defined(cateVarRatioMinMACVecExclude) then "--cateVarRatioMinMACVecExclude=~{cateVarRatioMinMACVecExclude} " else " "} \
        ~{if defined(cateVarRatioMaxMACVecInclude) then "--cateVarRatioMaxMACVecInclude=~{cateVarRatioMaxMACVecInclude} " else " "} \
        ~{if defined(dosageZerodCutoff) then "--dosageZerodCutoff=~{dosageZerodCutoff} " else " "} \
        ~{if defined(MACCutoff_to_CollapseUltraRare) then "--MACCutoff_to_CollapseUltraRare=~{MACCutoff_to_CollapseUltraRare} " else " "} \
        ~{if defined(DosageCutoff_for_UltraRarePresence) then "--DosageCutoff_for_UltraRarePresence=~{DosageCutoff_for_UltraRarePresence} " else " "} \
        ~{if defined(X_PARregion) then "--X_PARregion=~{X_PARregion} " else " "} \
        ~{if defined(is_rewrite_XnonPAR_forMales) then "--is_rewrite_XnonPAR_forMales=~{is_rewrite_XnonPAR_forMales} " else " "} \
        ~{if defined(weights_for_G2_cond) then "--weights_for_G2_cond=~{weights_for_G2_cond} " else " "} \
        ~{if defined(condition) then "--condition=~{condition} " else " "} \
        ~{if defined(minMAF) then "--minMAF=~{minMAF} " else " "} \
        ~{if defined(minMAC) then "--minMAC=~{minMAC} " else " "} \
        ~{if defined(maxMAFforGroupTest) then "--maxMAFforGroupTest=~{maxMAFforGroupTest} " else " "} \
        ~{if defined(minInfo) then "--minInfo=~{minInfo} " else " "} \
        ~{if defined(method) then "--method=~{method} " else " "} \
        ~{if defined(kernel) then "--kernel=~{kernel} " else " "} \
        ~{if defined(weights_beta_rare) then "--weights.beta.rare=~{weights_beta_rare} " else " "} \
        ~{if defined(weights_beta_common) then "--weights.beta.common=~{weights_beta_common} " else " "} \
        ~{if defined(weightMAFcutoff) then "--weightMAFcutoff=~{weightMAFcutoff} " else " "} \
        ~{if defined(r_corr) then "--r.corr=~{r_corr} " else " "} \
        ~{if help then "--help " else " "} \
        --SAIGEOutputFile=~{chrom}.txt

        gzip ~{chrom}.txt
	>>>

    runtime {
	docker: docker
	memory: memory + " GiB"
	disks: "local-disk " + disk + " HDD"
        cpu: cpu
	preemptible: preemptible
        maxRetries: maxRetries
    }

    output {
	File saige_output_file = chrom + ".txt.gz"
    }
}

task combine_saige_results {
    input {
        Array[File] saige_result_files
        String phenoCol

        String docker
        Float memory = 3.5
	Int? disk_size_override
        Int cpu = 1
        Int preemptible = 1
    	Int maxRetries = 0
    }
    Float saige_files_size = size(saige_result_files, "GiB")
    Int disk = select_first([disk_size_override, ceil(10.0 + 2.0 * saige_files_size)])

    command <<<
        zcat ~{sep=' ' saige_result_files} >> saige_~{phenoCol}_results_merged.txt
    >>>

    runtime {
	  docker: docker
          memory: memory + " GiB"
	  disks: "local-disk " + disk + " HDD"
          cpu: cpu
	  preemptible: preemptible
          maxRetries: maxRetries
	}

    output {
        File merged_saige_file = "saige_" + phenoCol + "_results_merged.txt"
    }
}

# QQ and Manhattan plots
task Plots {
  input {
    File merged_saige_file
    String docker
    Float memory = 16.0
    Int disk = 200
    Int cpu = 4
    Int preemptible = 1
    Int maxRetries = 0
  }
  String merged_saige_file_prefix = basename(merged_saige_file, ".txt")
  String merged_saige_file_in_current_dir = basename(merged_saige_file)

  # Plots are produced for each phenotype.
  # For each phenotype, a file containing all of the hits from Step 2 is output.
  # For each phenotype, a file containing a subset of all of the hits where "p < 0.5" from Step 2 is output.
  command <<<
    awk '$13 < 0.5' ~{merged_saige_file} >> ~{merged_saige_file_prefix}_subset.txt
    mv ~{merged_saige_file_prefix}_subset.txt .
    mv ~{merged_saige_file} .

    R --no-save --args ~{merged_saige_file_in_current_dir} <<RSCRIPT
    library(data.table)
    library(qqman)
    args <- commandArgs(trailingOnly = TRUE)
    for (file in args) {
      saige_output <- fread(file, drop=c("Allele1", "Allele2", "AC_Allele2", "AF_Allele2", "imputationInfo", "SE", "Tstat", "Is.SPA.converge", "varT", "varTstar"))
      saige_subset <-subset.data.frame(saige_output)
      saige_subset[,"CHR"] <-as.numeric(unlist(saige_subset[,"CHR"]))
      saige_subset[,"p.value"] <-as.numeric(unlist(saige_subset[,"p.value"]))
      saige_subset[,"POS"] <-as.numeric(unlist(saige_subset[,"POS"]))
      saige_complete_subset = saige_subset[complete.cases(saige_subset), ]
      qq_plot = substr(file,1,nchar(file)-4)
      qq_plot = paste(qq_plot, "qqplot.png", sep="_")
      png(qq_plot)
      print(qq(as.numeric(unlist(saige_complete_subset[,"p.value"]))))
      dev.off()
      manhattan_plot = substr(file,1,nchar(file)-4)
      manhattan_plot = paste(manhattan_plot, "manhattan.png", sep="_")
      png(manhattan_plot)
      print(manhattan(saige_complete_subset, chr="CHR", bp="POS", snp="SNPID", p="p.value", annotatePval = 1E-5))
      dev.off()
      }
    RSCRIPT
  >>>

  output {
        Array[File] output_plots = glob("*.png")
        Array[File] output_saige = glob("*.txt")
  }

  runtime {
        docker: docker
        memory: memory + " GiB"
	disks: "local-disk " + disk + " HDD"
        cpu: cpu
        preemptible: preemptible
        maxRetries: maxRetries
  }
}
