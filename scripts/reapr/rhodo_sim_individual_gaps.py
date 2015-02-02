import sys,os,subprocess

import argparse
from genomics_tools.file_formats import bam
from mapping import align
from genomics_tools.file_formats.various_annotations import to_AGP,to_GFF

def ReadInContigseqs(contigfile,contig_filter_length):
    cont_dict = {}
    k = 0
    temp = ''
    accession = ''
    for line in contigfile:
        if line[0] == '>' and k == 0:
            accession = line[1:].strip().split()[0]
            cont_dict[accession] = ''
            k += 1
        elif line[0] == '>':
            cont_dict[accession] = temp
            temp = ''
            accession = line[1:].strip().split()[0]
        else:
            temp += line.strip()
    cont_dict[accession] = temp

    # for contig in cont_dict.keys():
    #     print 'Initial length:', len(cont_dict[contig])
    if contig_filter_length:
        singled_out = 0
        for contig in cont_dict.keys():
            if len(cont_dict[contig]) < contig_filter_length:
                del cont_dict[contig]
                singled_out += 1
    return(cont_dict)


def simulate_instance(args):
    print 'Started modyfiyng genome'
    if not os.path.exists(args.output_path):
        os.makedirs(args.output_path)

    # genome_path = args.genome
    contig_path = os.path.join(args.output_path, 'ctgs.fa')
    read1_path = args.read1 
    read2_path = args.read2 
    bam_path = os.path.join(args.output_path, 'mapped')
    gff_path = os.path.join(args.output_path, 'true_error_pos.gff')
    gff_file = open(gff_path,'w')
    genome_seqs = ReadInContigseqs(open(args.genome, 'r'),10)


    #contigs/scaffolds
    gap = args.gapsize
    error = args.error
    chunk_size = 20000 
    modified_genome = {}
    modified_chunks = []
    for acc,seq in genome_seqs.iteritems():
        if acc == 'sequence_0':
            pos = 0
            chunks = [seq[i:i+chunk_size] for i in range(0, len(seq), chunk_size)]
            i=0

            for sample in range(100):
                N_s = 'N'* max(0,(gap + error))
                cut_size = gap + max(0,-error)
                modified_chunk = chunks[i][: len(chunks[i])-(cut_size)] + N_s
            #print modified_chunk
                modified_chunks.append(modified_chunk) 
                i+=1
                #print len(modified_chunk)
    
                pos += len(modified_chunk)

                if (gap + error) > 0:
                    error_start = pos - len(N_s)  
                    error_stop = pos  # error is anywhere in the introduced gap (either contraction or expansion)
                else:
                    error_start = pos 
                    error_stop = pos + 1 # error is at a specific position where a contraction has occured
                if error < 0:
                    to_GFF(gff_file, '{0}'.format(acc), 'TRUTH','FCD', error_start, error_stop, 1, '+', '.', 'Note=Error:Contraction {0}bp'.format(abs(error)))
                else:
                    to_GFF(gff_file, '{0}'.format(acc), 'TRUTH','FCD', error_start, error_stop, 1, '+', '.', 'Note=Error:Expansion {0}bp'.format(abs(error)))

            mod_seq = ''.join(modified_chunks)
            if error < 0:
                modified_genome['scf_gap{0}_errorsize_minus{1}'.format(gap,error)] = mod_seq
            else:
                modified_genome['scf_gap{0}_errorsize{1}'.format(gap,error)] = mod_seq

        else:
            modified_genome[acc] = seq

        #print and map

        ctgs = open(contig_path,'w')
        for acc,seq in modified_genome.iteritems():
            ctgs.write('>{0}\n{1}\n'.format(acc,seq))
        ctgs.close()
        ctgs = open(contig_path,'r')

        print 'Started mapping'
        align.bwa_mem(read1_path, read2_path, contig_path, bam_path, args)


def main(args):
    successful_experiments = 0
    while successful_experiments < 1: 
        try:
            simulate_instance(args)
        except subprocess.CalledProcessError:
            print 'Error'
            continue

        successful_experiments += 1
	

if __name__ == '__main__':
    ##
    # Take care of input
    parser = argparse.ArgumentParser(description="Simulate paired end reads with a given mean and standard deviation on mean insert size. A given coverage should also be provided.")
    parser.add_argument('genome', type=str, help='Path to genome. ')
    parser.add_argument('read1', type=str, help='Path to reads1. ')
    parser.add_argument('read2', type=str, help='Path to reads2. ')
    parser.add_argument('output_path', type=str, help='path to folder output. ')
    parser.add_argument( 'gapsize', dest='gapsize', type=int, default=False, help='Gapsize' )
    parser.add_argument( 'error', dest='error', type=int, default=False, help='Error size' )

    parser.add_argument( '-sort', dest='sort', action='store_false', default=True, help='Coordinate sort the reads in the bam file' )
    parser.add_argument( '-sam', dest='sam', action='store_true', default=False, help='Output a samfile (default is bam)' )
    #parser.add_argument( '-errors', dest='errorsize', type=int, nargs='+', default=False, help='gap distance error' )
    #parser.add_argument( '-nrgaps', dest='nrgaps', type=int, default=False, help='Number of gaps' )
    parser.add_argument('--threads', type=str, dest='threads', default='8', required=False, help='Number of threads for bwa mem.')


    args = parser.parse_args()

    main(args)








