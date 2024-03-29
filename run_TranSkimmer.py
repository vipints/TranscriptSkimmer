#!/usr/bin/env python
"""
wrapper for transcriptskimmer  

requirements:
    biopython :- http://biopython.org 
    gfftools  :- http://github.com/vipints/gfftools  
"""

import os 
import sys 
import numpy 
import pysam 
import shutil
import subprocess
from Bio import SeqIO 
from utils import GFFParser, helper 
from collections import defaultdict

from optparse import OptionParser, OptionGroup

def main():
    """
    main inputs to TranscriptSkimmer
    """

    #TODO options to the trsk program 
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

    ## genome preparation step 
    gio_path_temp = os.path.join(options.result_dir, "temp_gio")
    make_gio(options.fasta_file, gio_path_temp)

    #TODO run the program 
    #options="-maxel %d -ss -reglen 0.66 -maxic %d -minic 20 -maxin %d -mm 2 -exm 3 -indt 150 -exd 20 -tf 0.5 -inscf 3 -excut 3 -toff 100 -el 15" % (max_exon_length, max_intergenic_region, max_intron_length)
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

    ## running the core transcript_skimmer engine 
    gtf_db = run_trsk(gio_path_temp, options.bam_file, options.result_dir, advance_options)

    ## check the consistency of splice sites and write the result to outfile 
    outfilename = os.path.join(options.result_dir, options.gff_file) 
    splice_site_check(options.result_dir, gtf_db, options.min_orf_len, options.fasta_file, outfilename)
    
    ## cleanup to the genome files 
    shutil.rmtree(gio_path_temp)


def splice_site_check(gff_name, min_orf_length, fafile, outFile):
    """
    checking the consistency of splice site on predicted transcript models
    """
    
    gff_content = GFFParser.Parse(gff_name) ## getting the transcript predictions
    spliced_cand = 0
    sing_exon_gen = 0
    transcripts_region = defaultdict(list)
    for gene_recd in gff_content: ## screening the spliced transcripts
        spliced_transcript = defaultdict(list)
        for idx, sub_rec in enumerate(gene_recd['transcripts']):
            try:
                exon_cnt = len(gene_recd['exons'][idx])
            except:
                continue

            if exon_cnt > 1: ## skipping the single-exon transcripts 
                orf_length = 0 
                for idk, ex in enumerate(gene_recd['exons'][idx]):
                    orf_length += ex[1]-(ex[0]-1)

                    if idk == 0:
                        spliced_transcript[(gene_recd['name'], sub_rec[0], gene_recd['strand'])].append((None, ex[1]))
                    elif exon_cnt-1 == idk:
                        spliced_transcript[(gene_recd['name'], sub_rec[0], gene_recd['strand'])].append((ex[0], None))
                    else:
                        spliced_transcript[(gene_recd['name'], sub_rec[0], gene_recd['strand'])].append((ex[0], ex[1]))

                if orf_length <= min_orf_length: ## min orf length for the transcripts 
                    del spliced_transcript[(gene_recd['name'], sub_rec[0], gene_recd['strand'])] ## clearing that transcript details
                    continue
                    
                spliced_cand +=1
            else:
                sing_exon_gen +=1 
                single_exon_len = 0 
                for idk, ex in enumerate(gene_recd['exons'][idx]):
                    single_exon_len = ex[1]-(ex[0]-1)

                if single_exon_len > 900:
                    spliced_transcript[(gene_recd['name'], sub_rec[0], gene_recd['strand'])].append((ex[0], ex[1]))

        if spliced_transcript: 
            transcripts_region[gene_recd['chr']].append(spliced_transcript)

    sys.stdout.write("...considering %d spliced transcripts\n" % spliced_cand)
    sys.stdout.write("discarding transcripts...\n\t%d transcripts with single exon\n" % sing_exon_gen)
    ## got the best transcript to be considerd, now starting the consensus seq check
    genemodels = splice_site_consensus(fafile, transcripts_region)
    write_gene_models_file(gff_content, genemodels, outFile)


def write_gene_models_file(gff_cont, gene_models, outFileName):
    """
    writing the filtered gene models to the result file
    """
    sys.stdout.write("writing filtered gene models to %s ...\n" % outFileName)
    true_genes = 0 
    true_transcripts = 0 
    out_fh = open(outFileName, "w")
    for recd in gff_cont:
        trans_indices = [] 
        for idx, sub_rec in enumerate(recd['transcripts']):
            if (recd['chr'], recd['name'], sub_rec[0], recd['strand']) in gene_models:
                trans_indices.append(idx)
        if trans_indices:
            true_genes += 1 
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
                    true_transcripts += 1 
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
    sys.stdout.write("...done\n")
    sys.stdout.write("number of genes considered  %d\n" % true_genes)
    sys.stdout.write("number of transcripts considered  %d" % true_transcripts)


def splice_site_consensus(fas_file, splice_region):
    """
    splice site consensus check
    """
    sys.stdout.write("splice site sequence consensus check started...\n")
    get_gene_models = defaultdict()
    splice_site_con = 0 
    fas_fh = helper.open_file(fas_file)
    for fas_rec in SeqIO.parse(fas_fh, "fasta"):
        if fas_rec.id in splice_region:
            for details in splice_region[fas_rec.id]:
                for genes, regions in details.items():
                    acc_cons_cnt = 0 
                    don_cons_cnt = 0 

                    for region in regions:
                        if genes[-1] == '+':
                            if region[0]:## acceptor splice site 
                                acc_seq = fas_rec.seq[int(region[0])-3:int(region[0])-1]
                                if str(acc_seq).upper() == "AG":
                                    acc_cons_cnt += 1 

                            if region[1]:
                                don_seq = fas_rec.seq[int(region[1]):int(region[1])+2]
                                if str(don_seq).upper() == "GT":
                                    don_cons_cnt +=1 

                        elif genes[-1] == '-':
                            if region[0]: ## donor splice site 
                                don_seq = fas_rec.seq[int(region[0])-3:int(region[0])-1]
                                don_seq = don_seq.reverse_complement()
                                if str(don_seq).upper() == "GT":
                                    don_cons_cnt +=1 
                            
                            if region[1]:
                                acc_seq = fas_rec.seq[int(region[1]):int(region[1])+2]
                                acc_seq = acc_seq.reverse_complement()
                                if str(acc_seq).upper() == "AG":
                                    acc_cons_cnt += 1 
                    ## check for half of the consensus sites 
                    if acc_cons_cnt > (len(regions)/2) and don_cons_cnt > (len(regions)/2):
                        get_gene_models[(fas_rec.id, genes[0], genes[1], genes[2])] = 1   
                    else:
                        splice_site_con +=1 
    fas_fh.close()
    sys.stdout.write("...considering %d best transcripts\n" % len(get_gene_models))
    sys.stdout.write("discarding transcripts...\n\t%d splice-site consensus sequence missing\n" % splice_site_con)


def run_trsk(gio_file, bam_file, res_path, options, tmp_gff_file="tmp_trsk_genes.gff", tmp_reg_file="tmp_trsk_regions.bed"):
    """
    run TransriptSkimmer for provided genome
    """
    try:
        subprocess.call(["infer_genes"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except:
        exit("Please make sure that the `infer_genes` binary is in your $PATH")
    
    if not os.path.exists(bam_file + ".bai"): ## check for the index 
        file_prefix, ext = os.path.splitext(bam_file)
        sorted_bam = "%s_sortbyCoord" % file_prefix

        sys.stdout.wrtite("trying to sort by the coordinates with output prefix as: %s\n" % sorted_bam)
        if not os.path.isfile("%s.bam" % sorted_bam):
            pysam.sort(bam_file, sorted_bam)
            
        sorted_bam = "%s.bam" % sorted_bam
        sys.stdout.wrtite("indexing %s\n" % sorted_bam)
        if not os.path.exists(sorted_bam + ".bai"):
            pysam.index(sorted_bam) 

        bam_file = sorted_bam 

    sys.stdout.wrtite("bam file using is %s" % bam_file)

    gio_file = "%s/genome.config"  % gio_file 
    if not os.path.isfile("%s" % gio_file):
        sys.stdout.wrtite("error: failed to fetch genome index object file from %s\n" % gio_file ) 
        sys.exit(-1)
        
    cli_trsk = "infer_genes -gio %s -bam %s -gff %s -reg %s %s" % (gio_file, bam_file, tmp_gff_file, tmp_reg_file, options)  
    sys.stdout.write('\trun TranscriptSkimmer as: %s\n' % cli_trsk)
    
    try:
        os.chdir(res_path) 
        process = subprocess.Popen(cli_trsk, shell=True) 
        returncode = process.wait()
        
        if returncode !=0:
            raise Exception, "Exit status return code = %i" % returncode
        sys.stdout.write("TranscriptSkimmer run finished\n")

    except Exception, e:
        sys.stdout.write("Error running TranscriptSkimmer.\n%s" %  str( e ))
        sys.exit(-1)

    return os.path.join(res_path, tmp_gff_file)


def make_gio(in_file_name, gio_path):
    """
    make_gio builds a genome information object for an input fasta file.
    takes 2 arguments:
    @args fasta_file: is the input file in fasta format
    @type fasta_file: str
    @args gio_path: genome information object will be written to
    @type gio_path: dir
    """

    try:
         f_in = open(in_file_name, "r")
    except Exception as msg:
         print(msg)
         exit("cannot open input file %s" % in_file_name)

    write_dna = 1 
    flat_path = os.path.join(gio_path, "genome")

    try:
        if os.path.exists(flat_path):
            print("directory %s exists already." % flat_path)
        else:
            os.makedirs(flat_path)
    except Exception as msg:
        print(msg)
        exit("cannot create path %s" % flat_path)

    f_out = None
    f_out_dna = None
    contig_list = []
    num_bases = 0

    for line in f_in:
        if line.isspace():
            print("warning: ignoring empty line in file %s" % in_file_name)
            continue 
        if line[0].isspace():
            exit("wrong format: leading white space in file %s" % in_file_name)

        if line.startswith(">"):
            if f_out != None:
                f_out.close()
            if f_out_dna != None:
                f_out_dna.close()

            contig_list.append(line[1:-1].split()[0])
            out_name = os.path.join(flat_path, contig_list[-1] + ".flat")
            out_dna_name = os.path.join(flat_path, contig_list[-1] + ".dna")
            
            try:
                f_out = open(out_name, "w")
                if write_dna==1:
                    f_out_dna = open(out_dna_name, "w")
                    f_out_dna.write(line)
            except Exception as msg:
                print(msg) 
                exit("cannot open file %s" % out_name)
        else:
            try:
                f_out.write(line[0:-1].lower())
                if write_dna==1:
                    f_out_dna.write(line.lower())
            except Exception as msg: 
                if f_out != None:
                    print(msg)
                    exit("cannot write to file %s" % out_name)
                else:
                    exit("improper input format. No header in first line")

            num_bases += len(line)-1
    f_out.close()

    try:
        print("creating file %s" % os.path.join(gio_path, "genome.config"))
        f_conf = open(os.path.join(gio_path, "genome.config"), "w")
        f_conf.write("BASEDIR " +  os.path.abspath(gio_path) +"\n\n")
        f_conf.write("CONTIGS " +  str(len(contig_list)) +"\n")
        for c in contig_list:
            f_conf.write(c + "\tgenome/" + c + ".flat\tgenome/" + c + ".dna\n")
        f_conf.write("\nALPHABET acgt\n\n")
        f_conf.write("ESTFILES 0\n\n")
        f_conf.write("CDNAFILES 0\n\n")
        f_conf.write("ANNOTATIONFILES 0\n")
        f_conf.close()
    except Exception as msg:
        print(msg) 
        exit("cannot create file %s" % os.path.join(gio_path, "genome.config"))



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
