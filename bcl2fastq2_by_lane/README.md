## Bcl2fastq2

### What does this applet do?

This applet converts BCL files, generated from Illumina sequencers, to Fastq files. During the process it also demultiplexes samples tagged with different barcodes. You can find more information on bcl2fastq2 here: http://support.illumina.com/downloads/bcl2fastq-conversion-software-v217.html.

### What are typical use cases for this applet?

Processing reads data produces by Illumina sequencers. It has been tested on MiSeq, HiSeq2000/2500, and HiSeq 4000 platforms.

### What data are required for this applet to run?

1. Lane data archive
This applet requires a tar archive of the BaseCalls lane directory from the sequencing run directory, i.e. Run_directory/Data/Intensitites/BaseCalls/L001. You should be able to generate the required archive by running the following tar command from the command-line console.
    
        $ tar -C run_directory_path -cf output_tar_file_path Data/Intensities/BaseCalls/L00N

2. Metadata archive
It also requires a tar archive of run metadata files. The metadata archive should include the following files: runParameters.xml, RunInfo.xml, RTAConfiguration.xml, RTAComplete.txt, Data/Intensities/s.locs, as well as the following directories: RTALogs, Recipe, Config. Example console command.
    
        $ tar -C run_directory_path -cvf output_tar_file_path \\
        runParameters.xml RunInfo.xml RTAConfiguration.xml \\
        RTALogs RTAComplete.txt \\
        Recipe Config \\
        Data/Intensities/s.locs
    
3. Barcodes file (optional)
If your lane contains multiplexed samples, you will also need to provide a barcodes file. This is a text file with each barcode sequence on a separate line of the file. If you do not provide a --use-bases-mask pattern, the applet will automatically infer it from the barcodes provided.
    
### What does this applet output?

This applet a separate set of reads files, in fastq format, for each of your samples. It also outputs a set of miscellaneous describing metadata about the bcl2fastq2 run.

### How does this applet work?

This applet works by automatically generating the sample sheet and use-bases-mask value necessary for bcl2fastq2, before running it. It has been specifically designed for the Stanford Sequencing Center pipeline.

It performs the following steps:
- Download all data file
- Create bcl2fastq2 sample sheet
- Infer use-bases-mask
- Convert bcl to fastq files using bcl2fastq2
- Upload files back to DNAnexus project
