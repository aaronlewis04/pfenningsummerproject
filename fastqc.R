library("fastqcr")

#connected via to the cmu cluster with ssfhs which i downloaded from a github repo following instructions from user @thallgren on this page : https://github.com/telepresenceio/telepresence/issues/1654#issuecomment-873538291
fastqc(fq.dir = "/Users/aaronlewis/sshfs/projects/pfenninggroup/rada/single_cell", # FASTQ files directory
       qc.dir = "/Users/aaronlewis/sshfs/projects/pfenninggroup/rada/single_cell", # Results directory
       threads = 4                    # Number of threads
)

#in the same folder after running this one function was an html file for each fastq sample with the qc figures for each.


