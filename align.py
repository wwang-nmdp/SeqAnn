

import subprocess
import config
import os

class align:

    def __init__(self, input):
        self.input = input
        self.inF = open(config.get_gene_data_file(), 'r')
        self.gene = Gene("")

    def writeInputFile(self):
        inF = open(self.input, "r")
        fasta = open(config.alignInput, "w")
        for line in inF:
            fasta.write(line)
        fasta.write("\n")
        fasta.write(self.gene.toFasta())
        fasta.close()

    def goto_end(self):
        try:
            while True:
                line = self.inF.readline()
                if line.startswith("//"):
                    break
        except(StopIteration):
            return

    def findSeq(self):
        seq = ""
        while True:
            line = self.inF.readline()
            if line.startswith("//"):
                break
            else:
                seqFragments = line.split()
                seqFragments = seqFragments[:-1]
                seq += ''.join(seqFragments)

        return seq

    def readGene(self, geneID):
        gene = Gene(geneID)
        # find index
        while True:
            line = self.inF.readline()
            if line.startswith("FT   UTR"):
                gene.saveUTR(line)
            elif line.startswith("FT   exon"):
                gene.saveExon(line)
            elif line.startswith("FT   intron"):
                gene.saveIntron(line)
            elif line.startswith("SQ"):
                gene.saveSeq(self.findSeq())
                break
            else:
                continue
        self.gene = gene
        self.inF.close()


    def searchGene(self, gene):
        while True:
            line = self.inF.readline()
            if line.find(gene) == -1:
                # go to end of this gene block
                self.goto_end()
            else:
                self.readGene(gene)
                break

    def runAlign(self):
        os.remove(config.alignOutput)
        subprocess.run(["bin/clustalo", "-i", config.alignInput , "--outfmt=clu", "-o",
                        config.alignOutput, "--wrap=50000"])


class Section:
    def __init__(self, length, name):
        self.length = length
        self.name = name


class Gene:
    data = []
    utr = 0
    exon = 1
    intron = 1
    seq = ""

    def __init__(self, name):
        self.name = name

    def saveUTR(self, line):
        data = line.split()
        dataIndex = data[2].split("..")
        if self.utr == 0:
            name = 'UTR3'
        else:
            name = 'UTR5'
        self.utr += 1
        section = Section(int(dataIndex[1])-int(dataIndex[0])+1, name)
        self.data.append(section)

    def saveExon(self, line):
        data = line.split()
        dataIndex = data[2].split("..")
        section = Section(int(dataIndex[1])-int(dataIndex[0])+1, 'Exon' + str(self.exon))
        self.data.append(section)
        self.exon += 1

    def saveIntron(self, line):
        data = line.split()
        dataIndex = data[2].split("..")
        section = Section(int(dataIndex[1]) - int(dataIndex[0]) + 1, 'Intron' + str(self.intron))
        self.data.append(section)
        self.intron += 1

    def saveSeq(self, seq):
        self.seq = seq

    def toFasta(self):
        fasta = ">sequence \n"
        fasta += self.seq
        return fasta

