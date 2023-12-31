#!/nfs/goldstein/software/python2.7.7/bin/python
# -*- coding: utf-8 -*-
"""
Create dragen alignment configuration file for single sample runs
"""

from datetime import datetime
def create_align_config(sample,RGLane,RGFCID,RGPrepID):

    raw_config="""

#================================================================================
# Dragen 2.1 Configuration File
#================================================================================
#
# Desc: 	Configuration file for Dragen 2.1 for alignment and variant
# 		calling from Casava 1.8 raw fastq files
#
#================================================================================

#================================================================================
# SAMPLE SETUP
#================================================================================

intermediate-results-dir = /staging/tmp
ref-dir = /staging/REF/b37_decoy/
fastq-file1 = {first_fastq1}
fastq-file2 = {first_fastq2}
fastq-offset = 33 		# For CASAVA1.8 samples
enable-auto-multifile = true

#================================================================================
# OUTPUT
#================================================================================

output-file-prefix = {sample_name}.{pseudo_prepid}.{RGFCID}.{RGLane}
output-format = BAM
output-directory = {output_dir}

#================================================================================
# READ GROUP INFO
#================================================================================

RGID = {RGFCID}.{RGLane}		# Read group ID
RGLB = {prepid} 			# Library
RGPL = ILLUMINA 		# Platform/technology
RGSM = {sample_name}  		# Sample Name
RGPU = {RGFCID}.{RGLane}		# Platform unit
RGCN = IGM  			# Sequencing Center
"""
    if sample.metadata['seqtime'] != '1970-1-1':
        raw_config+="RGDT = {seqtime}		# Date the run was produced.  Format:ISO_8601\n"
    raw_config+="""
#================================================================================
# ALIGNMENT
#================================================================================

enable-map-align-output = true
enable-bam-indexing = true
enable-sort = true
enable-duplicate-marking = true
remove-duplicates = false
enable-sampling = true 			# automatically detect paired-end parameters with aligner test
enable-deterministic-sort = true 	# ensure sort order is completely repeatable at cost of a small decrease in speed

[Aligner]

match-score = 1 	# Score increment for matching reference nucleotide
mismatch-pen = 4 	# Score penalty for a mismatch
gap-open-pen = 6 	# Score penalty for opening a gap (insertion or deletion)
gap-ext-pen = 1 	# Score penalty for gap extension
unclip-score = 5 	# Score bonus for reaching each edge of the read
global = 0 		# If alignment is global (N-W) rather than local (S-W)
pe-orientation = 0 	# Expected paired-end orientation: 0=FR, 1=RF, 2=FF
pe-max-penalty = 60 	# Maximum pairing score penalty, for unpaired or distant ends
mapq-max = 60 		# Ceiling on reported MAPQ
supp-aligns = 3 	# Maximum supplimentary (chimeric) alignments to report per read
sec-aligns = 0 		# Maximum secondary (suboptimal) alignments to report per read
supp-as-sec = 0 	# If supplementary alignments should be reported with secondary flag
hard-clips = 6 		# Flags for hard clipping: 0 primary, 1 supplementary, 2 secondary
unpaired-pen = 80 	# Penalty for unpaired alignments in Phred scale
dedup-min-qual = 15 	# Minimum base quality for calculating read quality metric for deduplication
no-unpaired = 0 	# If only properly paired alignments should be reported for paired reads

#================================================================================
# MAPPER
#================================================================================

[Mapper]

seed-density = 0.5 	# Requested density of seeds from reads queried in the hash table
edit-mode = 0 		# 0 = No edits, 1 = Chain len test, 2 = Paired chain len test, 3 = Edit all std seeds
edit-seed-num = 6 	# For edit-mode 1 or 2: Requested number of seeds per read to allow editing on 
edit-read-len = 100 	# For edit-mode 1 or 2: Read length in which to try edit-seed-num edited seeds
edit-chain-limit = 29 	# For edit-mode 1 or 2: Maximum seed chain length in a read to qualify for seed editing
map-orientations = 0 	# 0=Normal, 1=No Rev Comp, 2=No Forward  (paired end requires Normal)


#================================================================================
# NOTES
#================================================================================
# Defaults taken from:
# Dragen User Defaults Configuration File - Revision 01.002.060.01.01.14.12343
# Adapter masking is performed during BCL
#
# http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/human/ #Patches to GRCh37 can be found here
# http://www.ncbi.nlm.nih.gov/variation/docs/multi_assembly_support/ #Multi Assembly Support
# http://exac.broadinstitute.org/faq
#================================================================================

#================================================================================
# FIN
#================================================================================
"""
    final_config = raw_config.format(RGPrepID=RGPrepID,
                                     RGLane=RGLane,
                                     RGFCID=RGFCID,
                                     prepid=RGPrepID,
                                     first_fastq1=sample.metadata['first_fastq1'],
                                     first_fastq2=sample.metadata['first_fastq2'],
                                     sample_name=sample.metadata['sample_name'],
                                     pseudo_prepid=sample.metadata['pseudo_prepid'],
                                     output_dir=sample.metadata['output_dir'],
                                     bed_file_loc=sample.metadata['bed_file_loc'],
                                     seqtime=sample.metadata['seqtime'])
    conf_file_fh = ("{script_dir}/{sample_name}.{pseudo_prepid}.DragenAlignment.{RGPrepID}.{RGFCID}.{RGLane}.conf"
               ).format(RGPrepID=RGPrepID,
                        RGFCID=RGFCID,
                        RGLane=RGLane,
                        **sample.metadata)
    sample.set('conf_file',conf_file_fh)
    conf_file = open(conf_file_fh,'w')
    conf_file.write(final_config)
    conf_file.close()



