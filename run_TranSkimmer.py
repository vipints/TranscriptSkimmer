#!/usr/bin/env python
"""
wrapper for transcriptskimmer  
"""

from __future__ import division
import os 
import sys 
import numpy 
import pysam 
import shutil
import subprocess
import collections
from Bio import SeqIO 
from gfftools import GFFParser, helper 

from optparse import OptionParser

def main():
    """
    main inputs to trsk 
    """

    parser = OptionParser()

    #TODO options to the trkm program 
    parser.add_option( "-g", "--genome_sequence", dest="genome_sequence", action="store", help="Genome sequence in fasta format")
    parser.add_option( "-b", "--bam_file", dest="bam_file", action="store", help="BAM file for storing sequence read alignment")
    parser.add_option( "-D", "--result_dir", dest="result_dir", help="Result directory; storing output files.")
    parser.add_option( "-A", "--gff_file", dest="gff_file", action="store", type="str", help="Assembled isoforms output file; formate is GFF", default="trsk_genes.gff")
    parser.add_option( "-R", "--bed_file", dest="bed_file", action="store", type="str", help="Assembled isoforms region file; formate is BED", default="trsk_regions.bed")
    parser.add_option( "-I", "--max_intron_len", dest="max_intron_len", action="store", type="int", help="max intron length", default=20000)
    parser.add_option( "-M", "--mm_filter", dest="mm_filter", action="store", type="int", help="mm_filter", default=2)
    parser.add_option( "-E", "--el_filter", dest="el_filter", action="store", type="int", help="el_filter", default=8)
    parser.add_option( "-X", "--max_exon_len", dest="max_exon_len", action="store", type="int", help="max exon length", default=8000)
    parser.add_option( "--min_exon_len", dest="min_exon_len", action="store", type="int", help="min exon length", default=10)
    parser.add_option( "--exon_mean", dest="exon_mean", action="store", type="int", help="exon_mean", default=5)
    parser.add_option( "--exon_term_thresh", dest="exon_term_thresh", action="store", type="int", help="exon_mean", default=3)
    parser.add_option( "--exon_drop", dest="exon_drop", action="store", type="int", help="exon drop", default=5)
    parser.add_option( "--exon_cut", dest="exon_cut", action="store", type="int", help="exon cut", default=3)
    parser.add_option( "--intron_cut", dest="intron_cut", action="store", type="int", help="exon cut", default=3)
    parser.add_option( "--intron_conf", dest="intron_conf", action="store", type="int", help="exon cut", default=1)
    parser.add_option( "--intron_dist", dest="intron_dist", action="store", type="int", help="exon cut", default=0)
    parser.add_option( "--intron_seed_conf", dest="intron_seed_conf", action="store", type="int", help="exon cut", default=3)
    #parser.add_option( "--reject_retained_introns", dest="reject_retained_introns", action="store", type="str", help="exon cut", default=0)
    parser.add_option( "--term_filter", dest="term_filter", action="store", type="float", help="term filter", default=2.0)
    #parser.add_option( "--find_orf", dest="find_orf", action="store", type="str", help="find orf", default=1)
    parser.add_option( "-O", "--min_orf_len", dest="min_orf_len", action="store", type="int", help="min orf length", default=300)
    parser.add_option( "--min_orf_sep", dest="min_orf_sep", action="store", type="float", help="min orf sep", default=0.7)
    parser.add_option( "--term_offset", dest="term_offset", action="store", type="int", help="term offset", default=100)
    parser.add_option( "--region_rel_length", dest="region_rel_length", action="store", type="float", help="region_rel_length", default=0.25)
    parser.add_option( "--min_intergenic_len", dest="min_intergenic_len", action="store", type="int", help="min_intergenic_len", default=50)
    parser.add_option( "--max_intergenic_len", dest="max_intergenic_len", action="store", type="int", help="max_intergenic_len", default=20000)
    #parser.add_option( "--intergenic_win", dest="intergenic_win", action="store", type="int", help="intergenic_win", default=100)
    parser.add_option( "-s", "--strand_specific", dest="strand_specific", action="store", type="str", help="cDNA library preparation protocol - strand specific (T) or non strand specific (F) default: F", default="F")

    ( options, args ) = parser.parse_args()
    
    #TODO check condition for input variables, specifically we need genome, bam file 
    if not ( options.genome_sequence and options.bam_file and \
        options.result_dir ):
        parser.print_help()
        sys.exit(-1)

    #TODO genome preparation step 
    gio_path_temp = os.path.join(options.result_dir, "temp_gio")
    make_gio(options.genome_sequence, gio_path_temp)

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
        " % ( options.max_intron_len, options.mm_filter, options.el_filter, options.max_exon_len, \
            options.min_exon_len, options.exon_mean, options.exon_term_thresh, options.exon_drop, options.exon_cut, \
            options.intron_cut, options.intron_conf, options.intron_dist, options.intron_seed_conf, \
            options.term_filter, options.min_orf_len, options.min_orf_sep, options.term_offset, options.region_rel_length, \
            options.min_intergenic_len, options.max_intergenic_len ) 

    if options.strand_specific in ['T', 't', 'True', 'true']:
        advance_options = "%s -ss" % advance_options
    elif options.strand_specific in ['F', 'f', 'False', 'false']:
        advance_options = "%s -nss" % advance_options

    #print advance_options
    gtf_db = run_trsk(gio_path_temp, options.bam_file, options.result_dir, advance_options)
    print gtf_db 

    #TODO clean the results 
    print "starting a check through the predicted gene models"
    fine_transcripts = validate_pred_gene_models(gtf_db, options.genome_sequence, options.gff_file)
    print "filtered transcripts are stored at %s" % fine_transcripts
    
    print "cleaning..."
    shutil.rmtree(gio_path_temp)
    os.unlink(gtf_db)

    # TODO clean the region file based on the new gff file 
    print "done"


def validate_pred_gene_models(gff_name, fas_file, out_fname):
    """
    check the sequence consistency/quality of predicted fragment

    @args gff_name: result file gff format from TranscriptSkimmer
    @type gff_name: str
    @args fas_file: genome sequence in fasta format
    @type fas_file: str 
    @args out_fname: filtered gene models in gff format 
    @type out_fname: str 
    """

    ## getting the genome annotation from GFF file 
    gff_content = GFFParser.Parse(gff_name)
    
    ## getting the spliced transcripts from the predicted gene list 
    transcripts_region = collections.defaultdict(list)
    for gene_recd in gff_content:
        spliced_transcript = collections.defaultdict(list)

        for idx, sub_rec in enumerate(gene_recd['transcripts']):
            exon_cnt = len(gene_recd['exons'][idx])

            ## skipping the single-exon transcripts 
            if exon_cnt > 1: 
                for idk, ex in enumerate(gene_recd['exons'][idx]):
                    if idk == 0:
                        ex[0] = None 
                    if exon_cnt-1 == idk:
                        ex[1] = None

                    spliced_transcript[(gene_recd['name'], sub_rec[0], gene_recd['strand'])].append(ex)

        transcripts_region[gene_recd['chr']].append(spliced_transcript)

    print "check for the splice site consensus for predicted transcripts"
    ## check for splice site consensus sequence of predicted transcripts 
    get_gene_models = collections.defaultdict()
    for fas_rec in SeqIO.parse(fas_file, "fasta"):
        if fas_rec.id in transcripts_region:
            for details in transcripts_region[fas_rec.id]:
                for genes, regions in details.items():

                    acc_cons_cnt = 0 
                    don_cons_cnt = 0 

                    for region in regions:
                        if genes[-1] == '+':
                            ## acceptor splice site 
                            if not numpy.isnan(region[0]):
                                acc_seq = fas_rec.seq[int(region[0])-3:int(region[0])-1]
                                if str(acc_seq).upper() == "AG":
                                    acc_cons_cnt += 1 

                            if not numpy.isnan(region[1]):
                                don_seq = fas_rec.seq[int(region[1]):int(region[1])+2]
                                if str(don_seq).upper() == "GT":
                                    don_cons_cnt +=1 

                        elif genes[-1] == '-':
                            ## donor splice site 
                            if not numpy.isnan(region[0]):
                                don_seq = fas_rec.seq[int(region[0])-3:int(region[0])-1]
                                don_seq = don_seq.reverse_complement()
                                if str(don_seq).upper() == "GT":
                                    don_cons_cnt +=1 
                            
                            if not numpy.isnan(region[1]):
                                acc_seq = fas_rec.seq[int(region[1]):int(region[1])+2]
                                acc_seq = acc_seq.reverse_complement()
                                if str(acc_seq).upper() == "AG":
                                    acc_cons_cnt += 1 
                    ## check for half of the consensus sites 
                    if acc_cons_cnt > (len(regions)/2) and don_cons_cnt > (len(regions)/2):
                        get_gene_models[(fas_rec.id, genes[0], genes[1], genes[2])] = 1   
    
    gff_cont = GFFParser.Parse(gff_name)

    ## filter out the best gene models based on the consensus 
    print "writing the fine tuned transctipts to the the file: %s " % out_fname 
    out_fh = open(out_fname, "w")
    for recd in gff_cont:
        trans_indices = [] 

        for idx, sub_rec in enumerate(recd['transcripts']):
            if (recd['chr'], recd['name'], sub_rec[0], recd['strand']) in get_gene_models:
                trans_indices.append(idx)

        if trans_indices:
            chr_name = recd['chr']
            strand = recd['strand']
            start = recd['start']
            stop = recd['stop']
            source = recd['source']
            ID = recd['name']
            Name = recd['gene_info']['Name']
            Name = ID if Name != None else Name  
            out_fh.write('%s\t%s\tgene\t%d\t%d\t.\t%s\t.\tID=%s;Name=%s\n' % (chr_name, source, start, stop, strand, ID, Name))
                
            for idz, tid in enumerate(recd['transcripts']):
                if idz in trans_indices:

                    t_start = recd['exons'][idz][0][0]
                    t_stop = recd['exons'][idz][-1][-1]
                    t_type = recd['transcript_type'][idz] 

                    out_fh.write('%s\t%s\t%s\t%d\t%d\t.\t%s\t.\tID=%s;Parent=%s\n' % (chr_name, source, t_type, t_start, t_stop, strand, tid[0], ID))
                    
                    for ex_cod in recd['utr5_exons'][idz]:
                        out_fh.write('%s\t%s\tfive_prime_UTR\t%d\t%d\t.\t%s\t.\tParent=%s\n' % (chr_name, source, ex_cod[0], ex_cod[1], strand, tid[0])) 
                    for ex_cod in recd['cds_exons'][idz]:
                        out_fh.write('%s\t%s\tCDS\t%d\t%d\t.\t%s\t%d\tParent=%s\n' % (chr_name, source, ex_cod[0], ex_cod[1], strand, ex_cod[2], tid[0])) 
                    for ex_cod in recd['utr3_exons'][idz]:
                        out_fh.write('%s\t%s\tthree_prime_UTR\t%d\t%d\t.\t%s\t.\tParent=%s\n' % (chr_name, source, ex_cod[0], ex_cod[1], strand, tid[0]))
                    for ex_cod in recd['exons'][idz]:
                        out_fh.write('%s\t%s\texon\t%d\t%d\t.\t%s\t.\tParent=%s\n' % (chr_name, source, ex_cod[0], ex_cod[1], strand, tid[0])) 
    out_fh.close()
    return out_fname


def run_trsk(gio_file, bam_file, res_path, options, tmp_gff_file="_tmp_trsk_genes.gff", tmp_reg_file="_tmp_trsk_regions.bed"):
    """
    run TransriptSkimmer for provided genome

    """

    print bam_file

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

    os.chdir(res_path) 

    cli_trsk = "infer_genes -gio %s -bam %s -gff %s -reg %s %s" % (gio_file, bam_file, tmp_gff_file, tmp_reg_file, options)  
    sys.stdout.write('\trun TranscriptSkimmer as: %s \n' % cli_trsk)

    process = subprocess.Popen(cli_trsk, shell=True) 
    process.wait()

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
