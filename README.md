# trajectoread_source

This repository contains the source files for DNAnexus applets developed for the Stanford Sequencing Center pipeline. These applets are used in the processing and delivery of sequencing reads data generated from Illumina sequencing platforms.

Currently, only the bcl2fastq_by_lane applet is publicly available. You can build the applet with the following steps.

1. Clone the trajectoread_source repo to your local or remote machine

        $ git clone <trajectoread_source_github_URL>
        $ cd trajectoread_source
        
2. Install and configure the DNAnexus software development kit (dx-toolkit). https://wiki.dnanexus.com/Downloads

3. Build the applet using dx-toolkit.
        
        $ dx build bcl2fastq2_by_lane

4. You can also deploy this applet in stand-alone form or as part of workflows using trajectoread_builder. Instructions for using trajectoread_builder can be found here: https://github.com/StanfordBioinformatics/trajectoread_builder.
        
        $ python builder.py -e production -a bcl2fastq2_by_lane
