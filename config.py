from enum import Enum

blast_output = "out/blast.txt"
alignInput = "alignMeta.txt"
alignOutput = "alignOutput.txt"
geneOutput = "geneOutput.txt"
geneType = ""

class GeneType(Enum):
    HLA_A="HLA-A"
    HLA_B="HLA-B"
    HLA_C="HLA-C"
    DRB1 = "HLA-DRB1"
    KIR = "KIR"


def get_blast_file():
    if geneType == GeneType.HLA_A:
        return "Subjects/A_gen.fasta"
    elif geneType == GeneType.HLA_B:
        return "Subjects/B_gen.fasta"
    elif geneType == GeneType.HLA_C:
        return "Subjects/C_gen.fasta"
    elif geneType == GeneType.DRB1:
        return "Subjects/DRB1_gen.fasta"
    elif geneType == GeneType.KIR:
        return "Subjects/KIR.fasta"
    else:
        return "Subjects/KIR.fasta"

def get_gene_data_file():
    if geneType == GeneType.HLA_A or geneType == GeneType.HLA_B or geneType == GeneType.HLA_C or geneType == GeneType.DRB1:
        return 'hla.dat'
    else :
        return 'kir.dat'

