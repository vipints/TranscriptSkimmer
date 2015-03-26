CXX = g++ -O9 -g
GENOME = genome
NGS = ngs_tools
NGS_TOOLS = ${NGS}/get_reads_direct.cpp ${NGS}/read.cpp
#TODO ask for the samtool dependency 
SAMDIR = /share/software/samtools/samtools-0.1.18_fPIC/
SAMINCL = -I${SAMDIR} -L${SAMDIR} -lbam -lz -lcurses

all: infer_genes

infer_genes: infer_genes_cmd.cpp infer_genes.h infer_genes.cpp region.cpp region.h gene.cpp gene.h gene_tools.h gene_tools.cpp Config.cpp Config.h
	$(CXX) infer_genes_cmd.cpp infer_genes.cpp region.cpp gene.cpp gene_tools.cpp Config.cpp ${NGS_TOOLS} ${GENOME}/genome.cpp ${GENOME}/tools.cpp -I${GENOME} -I${NGS} ${SAMINCL} -o infer_genes
