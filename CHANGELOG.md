# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## [1.0.0-beta](https://github.com/nasa/GeneLab_Metagenomics_Workflow/tree/NF_MetagenomeSeq_1.0.0-beta/)

This is the initial release of the NF_MetagenomeSeq workflow which is an extension of the previous
[NF_MGIllumina workflow](https://github.com/nasa/GeneLab_Data_Processing/blob/master/Metagenomics/Illumina/Workflow_Documentation/NF_MGIllumina/).

### Added
- Add two additional taxonomic profiling tools to the Read-based processing subworkflow
  - Kaiju taxonomic profiling
  - Kraken2 taxonomic profiling
- Add downstream analysis and visualization for both Read-based and Assembly-based processing outputs
  - Feature filtering for all output datatypes
  - Barplots or Heatmaps for each output datatype
- Add low biomass metagenomics processing support for both short-read (Illumina) and long-read (Nanopore) data
  - Implement short-read low biomass pipeline [GL-DPPD-7117](https://github.com/nasa/GeneLab_Data_Processing/blob/master/Metagenomics/Low_Biomass/Pipeline_GL-DPPD-7117_Versions/GL-DPPD-7117.md)
  - Implement long-read low biomass pipeline [GL-DPPD-7116](https://github.com/nasa/GeneLab_Data_Processing/blob/master/Metagenomics/Low_Biomass/Pipeline_GL-DPPD-7116_Versions/GL-DPPD-7116.md)
    - Long-read specific pre-processing
    - Long-read specific steps in the Assembly-based processing subworkflow
  - Read decontamination/filtering during pre-processing for both short- and long-read data
  - Feature decontamination using decontam package during downstream analysis for low biomass data
- Long-read data support for processing standard metagenomics data

### Changed
- Update to the latest standard short-read pipeline version [GL-DPPD-7107-B](https://github.com/nasa/GeneLab_Data_Processing/blob/master/Metagenomics/Illumina/Pipeline_GL-DPPD-7107_Versions/GL-DPPD-7107-B.md) 
of the GeneLab Metagenomics consensus processing pipelines.
- Replace bbduk with fastp for initial read quality filtering and adapter trimming in standard Illumina workflow

<BR>

---

> ***Note:** All previous workflow changes were associated with the previous versions of the GeneLab Metagenomics Standard Illumina Pipeline and can be found in the main [GeneLab_Data_Processing](https://github.com/nasa/GeneLab_Data_Processing) github repository in either the [NF_MGIllumina change log](https://github.com/nasa/GeneLab_Data_Processing/blob/master/Metagenomics/Illumina/Workflow_Documentation/NF_MGIllumina/CHANGELOG.md) for pipeline version [GL-DPPD-7101-A](https://github.com/nasa/GeneLab_Data_Processing/blob/master/Metagenomics/Illumina/Pipeline_GL-DPPD-7107_Versions/GL-DPPD-7107-A.md) or the [SW_MGIllumina change log](https://github.com/nasa/GeneLab_Data_Processing/blob/master/Metagenomics/Illumina/Workflow_Documentation/SW_MGIllumina/CHANGELOG.md) for pipeline version [GL-DPPD-7101](https://github.com/nasa/GeneLab_Data_Processing/blob/master/Metagenomics/Illumina/Pipeline_GL-DPPD-7107_Versions/GL-DPPD-7107-A.md) 