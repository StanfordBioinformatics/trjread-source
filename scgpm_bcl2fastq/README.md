# SCGPM Bcl2fastq

## What does this app do?
This app converts BCL files generated from Illumina platforms into FASTQ files. It is used as part of the standard Stanford Center for Genomics and Personalized Medicine (http://med.stanford.edu/gssc.html) sequencing pipeline.
The app will process data from a single Illumina flowcell lane.

## What are typical use cases for this app?
After sequencing your library on an Illumina platform, use this app to convert the output into FASTQ files.

## What data are required for this app to run?

This app requires: 
- A tar archive with sequencing lane files.
- A tar archive with sequencing metadata files. 
- (optional) A text file with a list of barcodes.


The tarball files are generated automatically as part of the GSSC pipeline. If your data has already been sequencing through GSSC, you can use those tarballs found in the /raw_data directory of our DNAnexus project.

If not, you can generate the required tarballs with the following commands:

### Generate sequencing archive
    $ cd <sequencing_run_directory>
    $ tar -cf lane_archive.tar Data/Intensities/BaseCalls/L00N 

### Generate metdata archive
    $ cd <sequencing_run_directory>
    $ tar -cvf metadata_archive.tar runParameters.xml RunInfo.xml RTAConfiguration.xml

## What does this app output?

- An array of FASTQ files.
- A lane.html file with basic library read statistics.
- A tools used files that describes the executables run to generate this data.
- (optional) A sample sheet describing the barcodes use for demultiplexing.

## How does this app work?

This app is essentially a wrapper for Illumina's bcl2fastq program. It does the extra work of automatically generating a sample sheet and base mask patterns used for demultipexing.

For more information, consult the manual at:

https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq2_guide_15051736_v2.pdf

