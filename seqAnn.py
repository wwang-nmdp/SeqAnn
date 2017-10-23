import argparse

import config
import blast
from config import GeneType
from align import align
from splitGene import splitGene


def main():
    parser = argparse.ArgumentParser(description='The seqAnn application')
    parser.add_argument('-i', action="store", dest="i", help="The input file is fasta with a single sequence")
    parser.add_argument('-g', action="store", dest="g", help="The gene type to process. It includes most of HLA and KIR types.For example: HLA-A, KIR")
    options = parser.parse_args()
    config.geneType = GeneType(options.g)
    print("Annotating\n")
    blast.blast(options.i, config.get_blast_file())
    print("The input sequence best matched with "+blast.get_match_gene()+"\n")
    alignHelp = align(options.i)
    alignHelp.searchGene(blast.get_match_gene())
    alignHelp.writeInputFile()
    alignHelp.runAlign()

    splitGeneHelp = splitGene(alignHelp.gene)
    splitGeneHelp.runSplit()


if __name__ == '__main__':
    main()


