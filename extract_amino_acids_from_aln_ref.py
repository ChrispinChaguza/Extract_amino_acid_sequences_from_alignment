#!/usr/bin/env python
import glob
import os
import sys
from Bio import SeqIO
from math import ceil
import argparse

def main():

    options=argparse.ArgumentParser(sys.argv[0],
    description='Script for extracting codon and amino acid sequences given amino acid position in a list of reference genomes in GenBank format (downloaded from NCBI or in NCBI-compliant format)',
        prefix_chars='-',
        add_help=True,
        epilog='Written by Chrispin Chaguza, Yale School of Public Health, Yale University, 2021')
    options.add_argument('-s','--sequences',action='store',required=True,nargs=1,
        metavar='input_sequences',dest='input_sequences',help='Input (multi-) fasta file containing aligned nucleotide sequences (alignment should be of same length as the reference genome)')
    options.add_argument('-r','--reference',action='store',required=True,nargs=1,
        metavar='input_reference',dest='input_reference',help='Input file containing locations to the reference genome (Genbank) to be used for annotation')
    options.add_argument('-c','--aa_position',action='store',required=True,nargs=1,
        metavar='amino_acid_position',dest='amino_acid_position',help='Tab-delimited input file containing the gene name (1st column) and amino acid position on each line (no header)')
    options.add_argument('-o','--output',action='store',required=True,nargs=1,
        metavar='output_file',dest='output_file',help='Output file containing the extracted nucleotide/amino acid sequences')
    options.add_argument('-v','--verbose',action='store_true',required=False,
        help='Additionally show output on the terminal')

    options=options.parse_args()

    input_sequences=options.input_sequences[0:][0]
    reference_genome=options.input_reference[0:][0]
    amino_acid_position=options.amino_acid_position[0:][0]
    output_file=options.output_file[0:][0]
    verbose=options.verbose

    alignFile = [i for i in SeqIO.parse(input_sequences,"fasta")]
    snpPositions = [str(i).strip().split("\t") for i in open(amino_acid_position,"r")]
    annotFile = SeqIO.read(reference_genome,"genbank")

    out_fhandle=open(output_file,"w")

    if verbose: 
        print("SAMPLE_ID\tGENOME_POS\tCODON_NUMBER\tCODON_START\tCODON_END\tREF_CODON\tALT_CODON\t\
              REF_AMINO\tALT_AMINO\tGENE\tNOTE")

    out_fhandle.write("SAMPLE_ID\tGENOME_POS\tCODON_NUMBER\tCODON_START\tCODON_END\tREF_CODON\tALT_CODON\t\
                       REF_AMINO\tALT_AMINO\tGENE\tNOTE\n")

    for geneName,snpPos1 in snpPositions:
        
        for seqNum,seqRec in enumerate(alignFile):
            
            for annotFeat in annotFile.features:
                start=0
                end=0
                snpPos1=int(snpPos1)

                if annotFeat.type == "CDS":
                    if annotFeat.location.start.position<=annotFeat.location.end.position:
                        start = annotFeat.location.start.position
                        end = annotFeat.location.end.position
                    else:
                        end = annotFeat.location.start.position
                        start = annotFeat.location.end.position

                    gene = ""
                    if 'gene' in annotFeat.qualifiers:
                        gene=annotFeat.qualifiers['gene'][0]
                    
                    note = ""
                    if 'note' in annotFeat.qualifiers:
                        note=annotFeat.qualifiers['note'][0]

                    snpPos = snpPos1
                    if geneName==gene:
                        if start<end:
                            snpPos=start+snpPos1*3
                        else:
                            snpPos=end+snpPos1*3

                        if snpPos>=start and snpPos<=end:
                            if geneName==gene:
                                if start<end:
                                    snpPos=start+snpPos1*3
                                else:
                                    snpPos=end+snpPos1*3

                                codonNum = snpPos1
                                codonStart = snpPos-3
                                codonEnd = snpPos
                                codonRef = str(annotFile.seq[codonStart:codonEnd])
                                codonAlt = str(seqRec.seq[codonStart:codonEnd])
                        
                                aminoRef = str(annotFile.seq[codonStart:codonEnd].translate())
                                aminoAlt = str(seqRec.seq[codonStart:codonEnd].translate())

                                SNP = str(seqRec[codonStart:codonEnd].translate())

                                if verbose: 
                                    print(str(seqRec.id)+"\t"+str(snpPos)+"\t"+str(codonNum)+"\t"+str(codonStart)+"\t"+str(codonEnd)+"\t"+\
                                          str(codonRef)+"\t"+str(codonAlt)+"\t"+str(aminoRef)+"\t"+str(aminoAlt)+"\t"+\
                                          str(gene)+"\t"+str(note))

                                out_fhandle.write(str(seqRec.id)+"\t"+str(snpPos)+"\t"+str(codonNum)+"\t"+str(codonStart)+"\t"+\
                                                  str(codonEnd)+"\t"+str(codonRef)+"\t"+str(codonAlt)+"\t"+str(aminoRef)+"\t"+\
                                                  str(aminoAlt)+"\t"+str(gene)+"\t"+str(note)+"\n")

                else:
                    pass

if __name__=="__main__":
    main()
