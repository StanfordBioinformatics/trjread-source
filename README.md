# trajectoread_source

## Overview
Source contains source files of DNAnexus Apps developed to process Illumina reads data. These applets are used in the processing and delivery of sequencing reads data generated from Illumina sequencing platforms.

## Functions
* Curate data with multi-dimensional metadata
* Use the following metadata types to make data searchable
        * Details: Fixed key-value pairs used to describe data provenance
        * Properties: Mutable key-value pairs that describe data to be used by downstream analyses
        * Tags: Mutable strings that tag data with phrases useful to end-users
* Use a standard dictionary for storing metadata

## Value
* Make it easy to find relevant data for downstream analyses
* Allow users to quickly generate sophisticated profiles of large data pools
* Allow users to do fine-grained data searching/filtering
* Describe & maintaing data provenance

1. Clone the trajectoread_source repo to your local or remote machine

        $ git clone <trajectoread_source_github_URL>
        $ cd trajectoread_source
        
2. Install and configure the DNAnexus software development kit (dx-toolkit). https://wiki.dnanexus.com/Downloads

3. Build the applet using dx-toolkit.
        
        $ dx build bcl2fastq2_by_lane

4. You can also deploy this applet in stand-alone form or as part of workflows using trajectoread_builder. Instructions for using trajectoread_builder can be found here: https://github.com/StanfordBioinformatics/trajectoread_builder.
        
        $ python builder.py -e production -a bcl2fastq2_by_lane
