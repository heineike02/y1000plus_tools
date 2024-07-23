import os
import subprocess
#from Bio import SeqIO


##define input directory
#eg: aln_dir = os.path.normpath('/home/heineikeb')


for alignment_full in os.listdir(base_dir + os.sep + os.path.normpath('msas/structural/tm_align/cds_aln')):
    alignment  = alignment_full.split('.')[0]
    print(alignment)

    #Run Mafft on proteome with auto parameter.  

    #Output the screen as a log so that I can extract the strategy as .aln.fasta.log
    print('Obtaining mafft alignment')
    mafft_command = ['mafft', '--auto', 
                     aln_dir + os.sep + alignment + '.fasta']

    mafft_output_fname =  aln_dir + os.sep + os.path.normpath('mafft/fasta/'+ alignment + '.aln.fasta')
    mafft_log_fname = aln_dir + os.sep + os.path.normpath('mafft/fasta/'+ alignment + '.aln.fasta.log')

    with open(mafft_output_fname, 'w') as mafft_output: 
        with open(mafft_log_fname, 'w') as mafft_log: 
            subprocess.run(mafft_command, stdout=mafft_output, stderr=mafft_log)