#!/usr/bin/env python
"""
wrapper for transcriptskimmer  

requirements:
    biopython :- http://biopython.org 

"""

import os 
import sys 
import numpy 
import pysam 
import shutil
import subprocess
import collections
from Bio import SeqIO 
from utils import GFFParser, helper 

from optparse import OptionParser, OptionGroup

def main():
    """
    main inputs to TranscriptSkimmer
    """

    #TODO options to the trkm program 
    parser = OptionParser(usage='usage: %prog [options] arguments')

    required = OptionGroup(parser, 'Required')

    required.add_option( "-g", "--fasta_file", dest="fasta_file", action="store", help="Genome sequence file in fasta format")
    required.add_option( "-b", "--bam_file", dest="bam_file", action="store", help="BAM file for storing sequence read alignment")
    required.add_option( "-D", "--result_dir", dest="result_dir", help="Result directory; storing output files.")

    optional = OptionGroup(parser, 'Optional')

    optional.add_option( "-A", "--gff_file", dest="gff_file", action="store", type="str", help="Assembled isoforms output file; formate is GFF", default="trsk_genes.gff")
    optional.add_option( "-R", "--bed_file", dest="bed_file", action="store", type="str", help="Assembled isoforms region file; formate is BED", default="trsk_regions.bed")
    optional.add_option( "-I", "--max_intron_len", dest="max_intron_len", action="store", type="int", help="max intron length", default=20000)
    optional.add_option( "-M", "--mm_filter", dest="mm_filter", action="store", type="int", help="mm_filter", default=2)
    optional.add_option( "-E", "--el_filter", dest="el_filter", action="store", type="int", help="el_filter", default=8)
    optional.add_option( "-X", "--max_exon_len", dest="max_exon_len", action="store", type="int", help="max exon length", default=8000)
    optional.add_option( "--min_exon_len", dest="min_exon_len", action="store", type="int", help="min exon length", default=10)
    optional.add_option( "--exon_mean", dest="exon_mean", action="store", type="int", help="exon_mean", default=5)
    optional.add_option( "--exon_term_thresh", dest="exon_term_thresh", action="store", type="int", help="exon_mean", default=3)
    optional.add_option( "--exon_drop", dest="exon_drop", action="store", type="int", help="exon drop", default=5)
    optional.add_option( "--exon_cut", dest="exon_cut", action="store", type="int", help="exon cut", default=3)
    optional.add_option( "--intron_cut", dest="intron_cut", action="store", type="int", help="exon cut", default=3)
    optional.add_option( "--intron_conf", dest="intron_conf", action="store", type="int", help="exon cut", default=1)
    optional.add_option( "--intron_dist", dest="intron_dist", action="store", type="int", help="exon cut", default=0)
    optional.add_option( "--intron_seed_conf", dest="intron_seed_conf", action="store", type="int", help="exon cut", default=3)
    #optional.add_option( "--reject_retained_introns", dest="reject_retained_introns", action="store", type="str", help="exon cut", default=0)
    optional.add_option( "--term_filter", dest="term_filter", action="store", type="float", help="term filter", default=2.0)
    #optional.add_option( "--find_orf", dest="find_orf", action="store", type="str", help="find orf", default=1)
    optional.add_option( "-O", "--min_orf_len", dest="min_orf_len", action="store", type="int", help="min orf length", default=300)
    optional.add_option( "--min_orf_sep", dest="min_orf_sep", action="store", type="float", help="min orf sep", default=0.7)
    optional.add_option( "--term_offset", dest="term_offset", action="store", type="int", help="term offset", default=100)
    optional.add_option( "--region_rel_length", dest="region_rel_length", action="store", type="float", help="region_rel_length", default=0.25)
    optional.add_option( "--min_intergenic_len", dest="min_intergenic_len", action="store", type="int", help="min_intergenic_len", default=50)
    optional.add_option( "--max_intergenic_len", dest="max_intergenic_len", action="store", type="int", help="max_intergenic_len", default=20000)
    #optional.add_option( "--intergenic_win", dest="intergenic_win", action="store", type="int", help="intergenic_win", default=100)
    optional.add_option( "-s", "--strand_specific", dest="strand_specific", action="store", type="str", help="cDNA library preparation protocol - strand specific (T) or non strand specific (F) default: F", default="F")

    parser.add_option_group(required)
    parser.add_option_group(optional)

    ( options, args ) = parser.parse_args()
    
    #TODO check condition for input variables, specifically we need genome, bam file 
    if not ( options.fasta_file and options.bam_file and \
        options.result_dir ):
        parser.print_help()
        sys.exit(-1)

    #TODO genome preparation step 
    gio_path_temp = os.path.join(options.result_dir, "temp_gio")
    make_gio(options.fasta_file, gio_path_temp)

    #options="-maxel %d -ss -reglen 0.66 -maxic %d -minic 20 -maxin %d -mm 2 -exm 3 -indt 150 -exd 20 -tf 0.5 -inscf 3 -excut 3 -toff 100 -el 15" % (max_exon_length, max_intergenic_region, max_intron_length)
    #TODO run the program 
    advance_options = "-maxin %d \
        -mm %d \
        -el %d \
        -maxel %d \
        -minel %d \
        -exm %d \
        -exts %d \
        -exd %d \
        -excut %d \
        -incut %d \
        -incf %d \
        -indt %d \
        -inscf %d \
        -tf %.2f \
        -orflen %d \
        -orfsep %.2f \
        -toff %d \
        -reglen %.2f \
        -minic %d \
        -maxic %d \
        " % ( options.max_intron_len, options.mm_filter, \
            options.el_filter, options.max_exon_len, \
            options.min_exon_len, options.exon_mean, \
            options.exon_term_thresh, options.exon_drop, \
            options.exon_cut, options.intron_cut, \
            options.intron_conf, options.intron_dist, \
            options.intron_seed_conf, options.term_filter, \
            options.min_orf_len, options.min_orf_sep, \
            options.term_offset, options.region_rel_length, \
            options.min_intergenic_len, options.max_intergenic_len ) 

    if options.strand_specific in ['T', 't', 'True', 'true']:
        advance_options = "%s -ss" % advance_options
    elif options.strand_specific in ['F', 'f', 'False', 'false']:
        advance_options = "%s -nss" % advance_options

    #print advance_options

    ## running the core transcript_skimmer engine 
    gtf_db = run_trsk(gio_path_temp, options.bam_file, options.result_dir, advance_options)
    
    
    #gtf_db 
    
    """
    #TODO clean the genome annotation based on the transcript model

    check the consistency of splice sites 

    read coverage to the promoter directions

    clean the region file based on the new gff file 
    """

    shutil.rmtree(gio_path_temp)


def run_trsk(gio_file, bam_file, res_path, options, tmp_gff_file="tmp_trsk_genes.gff", tmp_reg_file="tmp_trsk_regions.bed"):
    """
    run TransriptSkimmer for provided genome
    """

    ## indexing the in bam file 
    if not os.path.exists(bam_file + ".bai"):
        ## TranscriptSkimmer need the bam file sorted by coordinates 
            
        print "bam file index command throw errors! Looking through the bam file"
        file_prefix, ext = os.path.splitext(bam_file)
        sorted_bam = "%s_sortbyCoord" % file_prefix

        print "trying to sort based by the coordinates with output prefix as: %s" % sorted_bam
        if not os.path.isfile("%s.bam" % sorted_bam):
            print 'sorting...'
            pysam.sort(bam_file, sorted_bam)
            
        sorted_bam = "%s.bam" % sorted_bam
        print "now creating the index for %s " % sorted_bam
        if not os.path.exists(sorted_bam + ".bai"):
            print 'indexing...'
            pysam.index(sorted_bam) 

        bam_file = sorted_bam 

    print "bam file using is %s" % bam_file

    gio_file = "%s/genome.config"  % gio_file 
    if not os.path.isfile("%s" % gio_file):
        print "error: failed to fetch genome index object file from %s" % gio_file
        sys.exit(-1)
        
    cli_trsk = "infer_genes -gio %s -bam %s -gff %s -reg %s %s" % (gio_file, bam_file, tmp_gff_file, tmp_reg_file, options)  
    sys.stdout.write('\trun TranscriptSkimmer as: %s \n' % cli_trsk)
    
    try:
        os.chdir(res_path) 
        process = subprocess.Popen(cli_trsk, shell=True) 
        returncode = process.wait()
        
        if returncode !=0:
            raise Exception, "Exit status return code = %i" % returncode
        sys.stdout.write("transcript_skimmer run finished\n")

    except Exception, e:
        print 'Error running TranscriptSkimmer.\n%s' %  str( e )
        sys.exit(-1)

    return "%s/%s" % (res_path, tmp_gff_file)


def make_gio(in_file_name, gio_path):
	"""
    make_gio builds a genome information object for an input fasta file. 

	takes 2 arguments:

	@args fasta_file: is the input file in fasta format
    @type fasta_file: str 
	@args gio_path: is the directory to which the genome information object will be written to
    @type gio_path: dir 
	"""

	try:
		f_in = file(in_file_name, "r")
	except Exception, msg:
		print msg
		print "cannot open infile '" + in_file_name + "'"
		sys.exit(1)
	   
	write_dna = 1 
	flat_path = os.path.join(gio_path, "genome")
	try:
		if os.path.exists(flat_path):
			print "directory " + flat_path + " exists already."
		else:
			os.makedirs(flat_path)

	except Exception, msg:
		print msg
		print "cannot create path '" + flat_path + "'"
		sys.exit(1)

	f_out = None
	f_out_dna = None
	contig_list = []
	num_bases = 0
	
	for line in f_in:

		if line.isspace():
			print "warning: wrong format. ignoring empty line in file '" + in_file_name + "'"
			continue

		if line[0].isspace():
			print "wrong format: leading white space in file '" + in_file_name + "'"
			sys.exit(1)
	
		if line.startswith(">"):
			
			if f_out != None:
				f_out.close()

			if f_out_dna != None:
				f_out_dna.close()

			contig_list.append(line[1:-1].split()[0])
			out_name = os.path.join(flat_path, contig_list[-1] + ".flat")
			out_dna_name = os.path.join(flat_path, contig_list[-1] + ".dna")
			try:
				f_out = file(out_name, "w")
				print "creating file '" + out_name + "'"
				if write_dna==1:
					f_out_dna = file(out_dna_name, "w")
					f_out_dna.write(line)
					print "creating file '" + out_dna_name + "'"

			except Exception, msg:
				print msg
				print "cannot open file '" + out_name + "'"
				sys.exit(1)
				
		else:
			try:
				f_out.write(line[0:-1].lower())
				if write_dna==1:
					f_out_dna.write(line.lower())
			except Exception, msg:
				if f_out != None:
					print msg
					print "cannot write to file '" +out_name + "'"
					sys.exit(1)
				else:
					print "improper input format. No header in first line"
					sys.exit(1)

			num_bases += len(line)-1

	f_out.close()

	try:
		print "creating file '" + os.path.join(gio_path, "genome.config") + "'"
		f_conf = file(os.path.join(gio_path, "genome.config"),"w")
		f_conf.write("BASEDIR " +  os.path.abspath(gio_path) +"\n\n")
		f_conf.write("CONTIGS " +  str(len(contig_list)) +"\n")
		for c in contig_list:
			f_conf.write(c + "\tgenome/" + c + ".flat\tgenome/" + c + ".dna\n")
		f_conf.write("\nALPHABET acgt\n\n")
		f_conf.write("ESTFILES 0\n\n")
		f_conf.write("CDNAFILES 0\n\n")
		f_conf.write("ANNOTATIONFILES 0\n")
		f_conf.close()
	except Exception, msg:
		print msg
		print "cannot create file '" + os.path.join(gio_path, "genome.config") + "'"
		sys.exit(1)
		

if __name__=="__main__":
    main()
