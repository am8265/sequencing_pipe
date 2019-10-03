# sequencing_pipe
<h1>Sequencing Pipelines</h1>
Contains both the BCL2FASTQ pipelines and Dragen alignment pipelines for the Columbia University Medical Center's Institute for Genomic 
Medicine.  

* The BCL2Fastq pipeline converts BCL data outputted from our Illumina Sequencers (NovaSeqs and Hi2XXX) to fastq.gz files, generates 
samples sheets, updates SequenceDB and archives fastqs.gz to an archive location.
* The Dragen alignment pipeline iterates over a MySQL queue on SequenceDB named dragen_queue.  Each sample is broken down into its 
constituent read groups and their corresponding fastqs are aligned and outputted as a bam.  Upon successful alignment of all read groups,
an entry in dragen_sample_metadata is either populated or updated with the dragen_sample_metadata file is_merged set to 0 and the comma 
delimited, lexicographically sorted, location of each bam read group is inserted into the component bam field.  The field is_merge == 0 is 
the trigger for downstream processing

<h2>Usage</h2>

* /nfs/goldstein/software/python3.6.1-x86_64_shared/bin/python /nfs/goldstein/software/sequencing_pipe/master/sequencing_pipe/bcl2fastq_pipeline_wrapper.py -f [FCILLUMID]

* /nfs/goldstein/software/python3.6.1-x86_64_shared/bin/python /nfs/goldstein/software/sequencing_pipe/master/sequencing_pipe/dragen_pipe.py -p [pseudo_prepid]

<h2>Custom Configurations</h2>
Edit config.ini

