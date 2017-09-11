#!usr/bin/env python
'''Convert one flowcell lane of illumina BCL files to demultiplexed FastQs.

'''

__author__ = 'pbilling@stanford.edu (Paul Billing-Ross)'

import os
import re
import ast
import sys
import dxpy
import glob
import time
import json
import shutil
import fnmatch
import datetime
import subprocess

from xml.etree import ElementTree

LOCAL_OUTPUT = 'output'
PROJECT_DXID = dxpy.PROJECT_CONTEXT_ID

def parse_applet_inputs(applet_inputs):
    '''Parse applet arguments into functional categories.

    Args:
        applet_inputs (dict): All applet inputs.

    Returns:
        tuple: 5-item tuple with inputs separated into functional categories.

    '''

    # Used for DNAnexus operations but not added to file properties.
    applet_keys = (
                   'project_dxid', 
                   'project_folder',
                   'lane_data_tar',
                   'metadata_tar',
                   'barcodes_file')
    
    # Sequencing library information & added to file properties.
    sample_keys = (
                   'run_name', 
                   'lane_index',
                   'library_name',
                   'project_name') 
    
    # Arguments passed to bcl2fastq2 executable & added to file properties.
    options_values = (
                      'barcode_mismatches', 
                      'tiles',
                      'use_bases_mask')
    
    # Flags passed to bcl2fastq2 executable & added to file properties.
    options_flags = (
                     'create_fastq_for_index_reads', 
                     'ignore_missing_bcls', 
                     'ignore_missing_filter', 
                     'ignore_missing_positions', 
                     'tiles',  
                     'with_failed_reads')

    applet_args = { 
                   key : value 
                   for key, value in applet_inputs.items() 
                   if key in applet_keys
                  }
    sample_args = { 
                   key : value 
                   for key, value in applet_inputs.items() 
                   if key in sample_keys
                  }
    options_dict = { 
                    key : value 
                    for key, value in applet_inputs.items() 
                    if key in options_values
                   }
    flags_dict = {
                  key : value 
                  for key, value in applet_inputs.items() 
                  if key in options_flags 
                  and value is True
                 }
    
    # Optional input to be added as file tags.
    if 'tags' in applet_inputs:
        tags = applet_inputs['tags']
    else:
        tags = []

    if 'properties' in applet_inputs.keys():
        logger.warning('Extra properties will overwrite existing values' +
                       'in sample_args')
        sample_args.update(applet_inputs['properties'])

    options_dict['output_dir'] = LOCAL_OUTPUT
    applet_args['project_dxid'] = PROJECT_DXID

    return applet_args, sample_args, options_dict, flags_dict, tags

def download_file(file_dxid):
    '''Download file from DX Object store

    Args: 
        file_dxid (str): DNAnexus ID of file to be downloaded.
    
    Returns: 
        str: Path to downloaded file.
    
    '''

    dx_file = dxpy.DXFile(file_dxid)
    filename = dx_file.describe()['name']
    dxpy.download_dxfile(dxid=dx_file.get_id(), filename=filename)
    return filename

def untar_file(filename):
    '''Extract tar archive.

    Args:
        filename (str): Name of tar archive.

    '''

    command = 'tar -xf %s --owner root --group root --no-same-owner' % filename
    create_subprocess(cmd=command, pipeStdout=False)

def get_flowcell_id(run_info_xml):
    '''Parse Flowcell ID from RunInfo.xml file.

    Args:
        run_info_xml (str): Name of local RunInfo.xml file.

    Returns:
        str: Flowcell ID.

    '''

    with open(run_info_xml, 'r') as XML:
        for line in XML:
            match = re.search(r'<Flowcell>(.+)</Flowcell>', line)
            if match:
                flowcell_id = match.group(1)
                return flowcell_id
        logger.error('Could not parse Flowcell ID from RunInfo.xml')

def truncate_flowcell_id(flowcell_id):
    '''Helper function to truncate flowcell IDs.

    Some flowcell IDs are weirdly formatted. This function gets the 5 
    character flowcell ID relevant for identification.
        
    HFNKGBBXX => HFNKG
    000000000-AMG8G => AMG8G

    Args:
        flowcell_id (str): Raw flowcell ID.

    Returns:
        str: A 5-character ID.

    '''

    elements = flowcell_id.split('-')
    if len(elements) == 2:
        trunc_flowcell_id = elements[1][:5]
    elif len(elements) == 1:
        trunc_flowcell_id = elements[0][:5]
    return trunc_flowcell_id

def configure_logger(name, file_handle=False):
    '''Configure logger object.
    
    Args: 
        file_handle (bool): Indicates whether to write log to file.

    Returns:
        logger: A formatted logger object.

    '''

    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")

    if file_handle:
        log_file = '{}.log'.format(name)
        LOG = logging.FileHandler(log_file)
        LOG.setLevel(logging.DEBUG)
        LOG.setFormatter(formatter)
        logger.addHandler(LOG)


    return logger

def create_subprocess(cmd, pipeStdout=False, checkRetcode=True):
    '''Creates a subprocess via a call to subprocess
    
    Popen with the argument 'shell=True', and pipes stdout and stderr. 
    Stderr is always  piped, but stdout can be turned off. If the argument 
    checkRetcode is True, which it is by default, then for any non-zero return 
    code, an Exception is raised that will print out the the command, stdout, 
    stderr, and the returncode when not caught. Otherwise, the Popen instance 
    will be return, in which case the caller must call the instance's 
    communicate() method (and not it's wait() method!!) in order to get the 
    return code to see if the command was a success. communicate() will return 
    a tuple containing (stdout, stderr). But at that point, you can then check 
    the return code with Popen instance's 'returncode' attribute.
    
    Args: 
        cmd (str): The command line for the subprocess wrapped in the subprocess. 
        pipeStdout (bool): True means to pipe stdout of the subprocess.
        checkRetcode (bool): See documentation in the description above for specifics.
    
    Returns: 
        A two-item tuple containing stdout and stderr, respectively.
    
    '''

    stdout = None
    if pipeStdout:
        stdout = subprocess.PIPE
        stderr = subprocess.PIPE
    popen = subprocess.Popen(cmd,shell=True,stdout=stdout,stderr=subprocess.PIPE)
    if checkRetcode:
        stdout,stderr = popen.communicate()
        if not stdout:  # Will be None if not piped
            stdout = ""
        stdout = stdout.strip()
        stderr = stderr.strip()
        retcode = popen.returncode
        if retcode:
            raise Exception("subprocess command '{cmd}' failed with returncode '{returncode}'.\n\nstdout is: '{stdout}'.\n\nstderr is: '{stderr}'.".format(cmd=cmd,returncode=retcode,stdout=stdout,stderr=stderr))
        return stdout,stderr
    else:
        return popen

def format_library_name(library_name):
    '''Remove datestamp from library name

    Remove artifact of SCGPM library naming.

    Args:
        library_name (str): Arbitrary name of sequencing library.
    
    Returns:
        formatted_name (str): The reformatted name. 
    '''

    elements = library_name.split('rcvd')
    stripped_name = elements[0].rstrip()
    formatted_name = re.sub(r"[^a-zA-Z0-9]+", "-", stripped_name)
    return formatted_name

class Bcl2fastqJob:
    '''Converts Illumina intensity files to fastqs.

    Run Illumina bcl2fastq software to generate fastq files
    for one lane of flowcell data.

    Args:
        run_name (str): Sequencing run name.
        lane_index: Index of Illumina flowcell lane (1-8).

    Attributes:
        run_name (str): Sequencing run name.
        lane_index (int): Index of Illumina flowcell lane (1-8).
        sample_sheet (str): Name of the CSV sample sheet to be generated.
    
    '''

    def __init__(self, run_name, lane_index):

        self.run_name = run_name
        self.lane_index = lane_index

    def create_sample_sheet(self, barcodes_file):
        '''Creates sample sheet for bcl2fastq.

        Creates a CSV formatted samplesheet with barcode and sample 
        information used in demultiplexing. Barcodes should be specified
        in a tab-separated text file. Each barcode should be on a separate 
        line where the barcode is the first value and the optional sample 
        name is the second.

        "
        ATCGG   SampleA 
        GCTCA   SampleB
        "

        Args:
            barcodes_file (str): n

        Returns:
            tuple: Name of sample sheet and dictionary of barcode:sample_id

        '''
        barcode_sample_dict = {}

        sample_sheet = '{}-L{}-samplesheet.csv'.format(
                                                       self.run_name, 
                                                       self.lane_index)
        SHEET = open(sample_sheet, 'w')
        SHEET.write('[Data]\n')
        SHEET.write('Lane,Sample_ID,index,index2\n')

        with open(barcodes_file, 'r') as CODES:
            for line in CODES:
                elements = line.split()

                # Parse barcodes
                barcode = elements[0].upper()
                indexes = barcode.split('-')
                i7 = indexes[0]
                if len(indexes) == 2:
                    i5 = indexes[1]
                elif len(indexes) > 2:
                    logger.error('Length of indexes is greater than 2')
                    logger.error(indexes)
                else:
                    i5 = None

                # Get sample_id if provided, use barcode if not
                if len(elements) == 2:
                    sample_id = elements[1]
                    barcode_sample_dict[barcode] = sample_id
                elif len(elements) == 1:
                    sample_id = barcode
                    barcode_sample_dict[barcode] = None

                # Write to sample sheet
                if i5:
                    SHEET.write('{},{},{},{}\n'.format(
                                                       self.lane_index, 
                                                       sample_id, 
                                                       i7, 
                                                       i5))
                else:
                    SHEET.write('{},{},{},,\n'.format(
                                                      self.lane_index, 
                                                      sample_id, 
                                                      i7))
        SHEET.close()

        return sample_sheet, barcode_sample_dict

    def run(self, tools_used_dict, options_dict, flags_dict):
        '''Run bcl2fastq2 program.

        Args:
            use_bases_mask (str): Description of barcode masking pattern
            tools_used_dict (dict): Descritpion of executables and 
                                    configurations used by App
            options_dict (dict): Qualitative app configuration data
            flags_dict (dict): Boolean app configuration data

        Returns:
            tuple: stdout,stderr information from running bcl2fastq2 executable

        '''
        
        command = 'bcl2fastq ' 

        for option, value in options_dict.items():
            option = option.replace('_','-')
            command += '--{} {} '.format(option, value)

        for flag in flags_dict.keys():
            flag = flag.replace('_','-')
            command += '--{} '.format(flag)

        logger.info('Running bcl2fastq v2 with command: {}'.format(command))

        tools_used_dict['commands'].append(command)
        stdout,stderr = create_subprocess(cmd=command, pipeStdout=True)

class GetUseBasesMaskJob:
    '''Calculates use_bases_mask value from barcodes & RunInfo.xml.

    The bases mask indicates how bases should be parsed into reads and indexes.
    This class will generate the bases mask automatically based on the 
    RunInfo.xml file and the provided barcodes. Unnecessary if 

    '''

    def __init__(self):
        pass

    def _get_use_bases_mask(self, run_info_file, index_lengths_list):
        '''Get --use-bases-mask value.

        Args:
            run_info_file (str): Filename of RunInfo.xml file ("RunInfo.xml")
            index_lengths_list (list): These lengths should all be the same

        Returns:
            str: Describes how to parse barcodes for demultiplexing

        '''

        mask_elements = {}       
        index_lengths_set = {
                             index_length 
                             for index_length in index_lengths_list 
                             if index_length != None
                            }
        unique_index_lengths = list(index_lengths_set)
        if len(unique_index_lengths) != 1:
            logger.error("Inconsistent index lengths: {}".format(unique_index_lengths))
            sys.exit()
        else:
            index_lengths = unique_index_lengths[0]


        tree = ElementTree.parse(run_info_file)
        for read_element in tree.findall(".//Read"):
            number = int(read_element.get('Number'))
            read_length = int(read_element.get('NumCycles'))
            is_indexed = read_element.get('IsIndexedRead')

            if is_indexed == 'Y':
                # RunInfo index numbers are 2 & 3. Corresponding index_length
                # indices are 0 & 1. Subtract 2 to get corresponding index length.
                index_length = int(index_lengths[number-2])
                if index_length > 0:
                    index_mask = 'I{}'.format(index_length)
                else:
                    index_mask = ''

                if read_length - index_length > 0:
                    n_mask = 'n{}'.format(read_length - index_length)
                elif read_length - index_length < 0:
                    logger.error(
                                 'Index length is longer than read length. ' + 
                                 'Read {}: {} Cycles. '.format(number, read_length) +
                                 'Index length: {}.'.format(index_length))
                    sys.exit()
                else:
                    n_mask = ''
                mask_elements[number] = '{}{}'.format(index_mask, n_mask) 
            
            elif is_indexed == 'N':
                mask_elements[number] = 'y{}'.format(read_length)
        
        sorted_mask_elements = [mask for number, mask in sorted(mask_elements.items())]
        use_bases_mask = ','.join(sorted_mask_elements)
        return use_bases_mask

    def _count_index_lengths(self, barcodes, index_lengths=[]):
        '''Recursive function to get index lengths of all barcodes.

        Args:
            barcodes (list): List of barcodes (i.e. i7_index-i5_index)
            index_lengths (list): List of index lengths (?)

        Returns:
            list: List of tuples with i7 and i5 index lengths

        '''

        try:
            barcode = barcodes.pop()
            logger.info('Getting barcode length: {}'.format(barcode))
        except:
            # All barcodes have been parsed
            logger.info('No more barcodes')
            return index_lengths

        elements = barcode.split('-')
        i7_length = len(elements[0])
        if len(elements) == 2:
            i5_length = len(elements[1])
        elif len(elements) > 2:
            logger.error('Barcode has more than two indexes: {}').format(barcode)
            sys.exit()
        else:
            i5_length = 0
        index_lengths.append((i7_length, i5_length))
        return self._count_index_lengths(barcodes, index_lengths)

    def run(self, barcodes_list, run_info_file):
        '''Calculate use-bases-mask from read & index lengths.

        Args:
            barcodes_list (list): List of barcodes.
            run_info_file (str): Path to RunInfo file.

        Returns:
            str: --use-bases-mask argument passed to bcl2fastq2.

        '''
        
        if len(barcodes_list) == 0:
            index_lengths = (0,0)
        else:
            index_lengths = self._count_index_lengths(barcodes_list)
        logger.info('Index lengths: {}'.format(index_lengths))

        use_bases_mask = self._get_use_bases_mask(run_info_file, index_lengths)
        logger.info('--use-bases-mask {}'.format(use_bases_mask))
        return use_bases_mask

class Bcl2fastqFileUploader:
    '''Upload bcl2fastq2 output files.

    Upload resulting fastq files as well as intermediate and accessory files.

    Args:
        project_dxid (str): ID of project where files will be uploaded.
        project_path (str): Folder path where files will be uploaded. 

    Attributes:
        project_dxid (str): ID of project where files will be uploaded.
        project_path (str): Folder path where files will be uploaded.

    '''

    def __init__(self, project_dxid, project_path):

        self.project_dxid = project_dxid
        self.project_path = project_path

    def upload_fastq_files(self, raw_properties, tags):
        '''Recursively find and upload all fastqs to DNAnexus object store.

        Args:
            raw_properties (dict): Properties with values of different types.
            tags (list): List of descriptive tags.

        Returns:
            list: DXLinks to all uploaded fastq files.
        
        '''

        fastq_dxlinks = []
        for fastq in self._find_fastqs():
            fastq_metadata = self._get_scgpm_fastq_name(
                                                        fastq, 
                                                        raw_properties['flowcell_id'], 
                                                        raw_properties['library_name'], 
                                                        raw_properties['lane_index'])
            scgpm_name = fastq_metadata[0] 
            barcode = fastq_metadata[1]
            read_index = fastq_metadata[2]
            
            # Convert all property values to strings
            properties = {key : str(value) for key, value in raw_properties.items()}
            properties['barcode'] = barcode
            properties['read_index'] = read_index

            project_folder = '{}/fastqs'.format(self.project_path)
            fastq_dxid = dxpy.upload_local_file(
                                                filename = fastq, 
                                                name = scgpm_name, 
                                                properties = properties, 
                                                tags = tags,
                                                project = self.project_dxid, 
                                                folder = project_folder,
                                                parents = True)
            fastq_dxlink = dxpy.dxlink(fastq_dxid)
            fastq_dxlinks.append(fastq_dxlink)
        return fastq_dxlinks

    def upload_sample_sheet(self, local_file_path, raw_properties):
        '''Upload sample sheet to DNAnexus project.

        Args:
            local_file_path (str): Local path of sample sheet.
            raw_properties (dict): Properties with values of different types.

        Returns:
            str: DXLink to sample sheet on DNAnexus object store.

        '''

        # Convert all property values to strings
        properties = {key : str(value) for key, value in raw_properties.items()}
        properties['file_type'] = 'sample_sheet'

        project_folder = '{}/miscellany'.format(self.project_path)
        sample_sheet_dxid = dxpy.upload_local_file(
                                                   filename = local_file_path,
                                                   properties = properties, 
                                                   project = self.project_dxid, 
                                                   folder = project_folder, 
                                                   parents = True)
        return dxpy.dxlink(sample_sheet_dxid)

    def upload_lane_html(self, raw_properties, tags):
        '''Upload lane.html file to DNAnexus project.
        
        Args:
            local_file_path (str): Local path of sample sheet.
            raw_properties (dict): Properties with values of different types.

        Returns:
            str: DXLink to lane.html file on DNAnexus object store.

        '''

        # Convert all property values to strings
        properties = {key : str(value) for key, value in raw_properties.items()}
        properties['file_type'] = 'lane_html'

        project_folder = '{}/miscellany'.format(self.project_path)

        local_file_path = (
                           '{}/Reports/html/'.format(LOCAL_OUTPUT) +  
                           '{}/all/all/all/lane.html'.format(properties['flowcell_id']))
        remote_file_name = '{}_L{}.lane.html'.format(
                                                     properties['run_name'], 
                                                     properties['lane_index'])
        lane_html_dxid = dxpy.upload_local_file(
                                                filename = local_file_path,
                                                name = remote_file_name,
                                                properties = properties, 
                                                tags = tags,
                                                project = self.project_dxid, 
                                                folder = project_folder, 
                                                parents = True)
        return dxpy.dxlink(lane_html_dxid)

    def upload_tools_used(self, tools_used_dict, raw_properties):
        '''Write console commands to Tools Used file & upload.

        Args:
            tools_used_dict (dict): Description of executables & configurations.
            raw_properties (dict): Properties with values of different types.

        Returns:
            str: DXLink to "tools used" file on DNAnexus object store.
        '''

        # Convert all property values to strings
        properties = {key : str(value) for key, value in raw_properties.items()}
        properties['file_type'] = 'lane_html'

        # Write file
        local_file_path = 'bcl2fastq_tools_used.json'
        with open(local_file_path, 'w') as TOOLS:
            TOOLS.write(json.dumps(tools_used_dict))

        # Upload file
        properties['file_type'] = 'tools_used'
        project_folder = '{}/miscellany'.format(self.project_path)
        tools_used_dxid = dxpy.upload_local_file(
                                                 filename = local_file_path, 
                                                 properties = properties, 
                                                 project = self.project_dxid, 
                                                 folder = project_folder, 
                                                 parents = True)
        return dxpy.dxlink(tools_used_dxid)

    def _find_fastqs(self):
        '''Recursively find all fastq files in directory path.

        Returns:
            list: Paths of found fastq files.

        '''
        
        matches = []
        # /home/dnexus is default location for dxpy.download_dxfile
        for root, dirnames, filenames in os.walk('/home/dnanexus'):
            for filename in fnmatch.filter(filenames, '*.fastq.gz'):
                matches.append(os.path.join(root, filename))
        logger.info('Found the following fastq files: {}'.format(matches))
        return matches

    def _get_scgpm_fastq_name(self, fastq_path, flowcell_id, library_name, lane_index):
        '''Get SCGPM formatted fastq name, barcode, and read index. 

        Args:
            fastq_path (str): Local path to fastq file.
            flowcell_id (str): Flowcell ID fastq was sequenced on.
            library_name (str): Arbitrary library name.
            lane_index (int): Flowcell lane index (1-8).

        Returns:
            tuple: (str) new fastq name; (str) fastq barcode; (int) read index (1/2)

        '''

        fastq_filename = os.path.basename(fastq_path)
        elements = fastq_filename.split('_')
        if len(elements) < 5 or len(elements) > 6:
            logger.error('Fastq filename has unusual number of elements: ' +
                         '{}'.format(fastq_filename))
            sys.exit()
        
        elif len(elements) == 6:
            # Dual barcode : TCTCGCGC_TCAGAGCC_S47_L001_R2_001.fastq.gz
            barcodes = (elements[0], elements[1])
            barcode = '%s-%s' % (barcodes[0], barcodes[1])
            read = elements[4]       
        
        elif len(elements) == 5:
            # No barcode : Undetermined_S1_L001_R1_001.fastq.gz
            # Single barcode: TCAGAGCC_S47_L001_R2_001.fastq.gz
            barcode = elements[0]
            read = elements[3]
            
        # Parse read index
        read_index_match = re.match(r'R(\d)', read)
        if read_index_match:
            read_index = read_index_match.group(1)
        else:
            print 'Could not determine read index: %s' % read
            
        trunc_flowcell_id = truncate_flowcell_id(flowcell_id)
        scgpm_name = 'SCGPM_%s_%s_L%d_%s_R%d.fastq.gz' % (
                                                          library_name, 
                                                          trunc_flowcell_id,
                                                          lane_index,  
                                                          barcode,
                                                          int(read_index))
        return(scgpm_name, barcode, read_index)

@dxpy.entry_point("main")
def main(**applet_inputs):
    '''Run bcl2fastq2 to generate and demultiplex fastq files.

    Use illumina bcl2fastq applet to perform demultiplex and 
    convert bcl files to fastq files. Currently handles files generated from
    RTA version 2.7.3 and earlier.

    Args:
        applet_input (dict): Input parameters specified when calling applet 
                             from DNAnexus.

    Returns:
        dict: Names of outputs and corresponding file dxids.

    '''

    # Define all variables here
    global logger
    logger = configure_logger(name = 'RunBcl2fastq2', file_handle = True)

    tools_used_dict = {'name': 'Bcl to Fastq Conversion and Demultiplexing', 'commands': []}
    output = {}

    # Parse applet inputs
    parsed_inputs = parse_applet_inputs(applet_inputs)
    applet_args = parsed_inputs[0]
    sample_args = parsed_inputs[1]
    options_dict = parsed_inputs[2]
    flags_dict = parsed_inputs[3]
    tags = parsed_inputs[4]

    # Determines whether or not to create sample sheet, use bases mask.
    if not 'barcodes_file' in applet_args.keys():
        barcodes = False
    else:
        barcodes = True

    # Download and untar lane files in /home/dnanexus
    logger.info('Downloading lane data archive: %s' % applet_args['lane_data_tar'])
    lane_tar_filename = download_file(applet_args['lane_data_tar'])
    untar_file(lane_tar_filename)
    logger.info('Downloading lane metadata archive: %s' % applet_args['metadata_tar'])
    metadata_filename = download_file(applet_args['metadata_tar'])
    untar_file(metadata_filename)
    if barcodes:
        logger.info('Downloading barcodes file: %s' % applet_args['barcodes_file'])
        barcodes_filename = download_file(applet_args['barcodes_file'])
    else:
        logger.info('No barcodes associated with this sample')

    # Create upload & bcl2fastq runner objects
    uploader = Bcl2fastqFileUploader(
                                     applet_args['project_dxid'], 
                                     applet_args['project_folder'])
    bcl_job = Bcl2fastqJob(
                           run_name = sample_args['run_name'], 
                           lane_index = sample_args['lane_index'])

    # Create sample sheet
    if barcodes:
        logger.info('Creating sample sheet')
        sample_sheet, barcode_sample_dict = bcl_job.create_sample_sheet(barcodes_filename)
        options_dict['sample_sheet'] = sample_sheet
        output['sample_sheet'] = uploader.upload_sample_sheet(sample_sheet, sample_args)
    else:
        logger.info('Skipping sample sheet generation; no barcodes.')
    
    # Parse use-bases-mask from RunInfo.xml and samplesheet barcodes
    if not 'use_bases_mask' in options_dict.keys() and barcodes:
        logger.info('Inferring use-bases-mask from barcodes')
        base_mask_job = GetUseBasesMaskJob()
        use_bases_mask = base_mask_job.run(
                                           barcodes_list = barcode_sample_dict.keys(),
                                           run_info_file = 'RunInfo.xml')
        options_dict['use_bases_mask'] = use_bases_mask
        
    else:
        logger.info('Not generating bases mask.')
    
    logger.info('Convert bcl to fastq files')
    bcl_job.run(
                tools_used_dict = tools_used_dict,
                options_dict = options_dict,
                flags_dict = flags_dict)

    # Get fastq metadata
    logger.info('Getting fastq metadata') 
    sample_args['flowcell_id'] = get_flowcell_id('RunInfo.xml')
    sample_args['trunc_flowcell_id'] = truncate_flowcell_id(sample_args['flowcell_id'])
    sample_args['library_name'] = format_library_name(sample_args['library_name'])
    
    # Combine multiple dictionaries to add to file properties
    fastq_properties = {}
    fastq_properties.update(sample_args)
    fastq_properties.update(options_dict)
    fastq_properties.update(flags_dict)

    # Call uploader to upload results files
    logger.info('Create tools used file')
    output['tools_used'] = uploader.upload_tools_used(tools_used_dict, fastq_properties)
    
    # Upload fastqs to DNAnexus
    logger.info('Uploading fastq files back to DNAnexus')  
    output['fastqs'] = uploader.upload_fastq_files(
                                                   raw_properties = fastq_properties,
                                                   tags = tags)
    output['lane_html'] = uploader.upload_lane_html(
                                                    raw_properties = fastq_properties,
                                                    tags = tags)
    return output

dxpy.run()
