# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## [1.0.0](https://github.com/nasa/GeneLab_Metagenomics_Workflow/tree/NF_MetagenomeSeq_1.0.0/)

This is the initial release of the NF_MetagenomeSeq workflow which is an extension of the previous
[NF_MGIllumina workflow](https://github.com/nasa/GeneLab_Data_Processing/blob/master/Metagenomics/Illumina/Workflow_Documentation/NF_MGIllumina/).

### Added
- Low biomass metagenomics processing for both short-read (Illumina) and long-read (Nanopore) data
  - long-read specific pre-processing
  - long-read specific updates to the Assembly-based processing subworkflow
  - read decontamination/filtering during pre-processing for both long- and short-read data
- Long-read data support for processing standard metagenomics data
- Additional taxonomic profiling tools in the Read-based processing subworkflow
  - Kaiju taxonomic profiling
  - Kraken2 taxonomic profiling
- Downstream analysis
  - Feature filtering for all output datatypes
  - Barplots or Heatmaps generated for each output datatype
- Feature decontamination during downstream analysis for low biomass data

### Changed
- Update to the latest standard short-read pipeline version [GL-DPPD-7101-B](https://github.com/nasa/GeneLab_Data_Processing/blob/master/Metagenomics/Illumina/Pipeline_GL-DPPD-7107_Versions/GL-DPPD-7107-B.md) 
of the GeneLab Metagenomics consensus processing pipelines.
- Replace bbduk with fastp for initial read quality filtering and adapter trimming

<BR>

---

> ***Note:** All previous workflow changes were associated with the previous versions of the GeneLab Metagenomics Standard Illumina Pipeline and can be found in the main [GeneLab_Data_Processing](https://github.com/nasa/GeneLab_Data_Processing) github repository in either the [NF_MGIllumina change log](https://github.com/nasa/GeneLab_Data_Processing/blob/master/Metagenomics/Illumina/Workflow_Documentation/NF_MGIllumina/CHANGELOG.md) for pipeline version [GL-DPPD-7101-A](https://github.com/nasa/GeneLab_Data_Processing/blob/master/Metagenomics/Illumina/Pipeline_GL-DPPD-7107_Versions/GL-DPPD-7107-A.md) or the [SW_MGIllumina change log](https://github.com/nasa/GeneLab_Data_Processing/blob/master/Metagenomics/Illumina/Workflow_Documentation/SW_MGIllumina/CHANGELOG.md) for pipeline version [GL-DPPD-7101](https://github.com/nasa/GeneLab_Data_Processing/blob/master/Metagenomics/Illumina/Pipeline_GL-DPPD-7107_Versions/GL-DPPD-7107-A.md) 