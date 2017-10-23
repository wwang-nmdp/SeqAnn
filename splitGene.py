import config
from align import Section
class splitGene:
    def __init__(self, gene):
        self.gene = gene
        self.map = dict()
        self.buildMap()

    def runSplit(self):
        print("Parsing exon/intron boundaries......\n")
        inF = open(config.alignOutput, 'r')
        lines = inF.readlines()
        inputSeq = lines[3].split()[1]
        reference = lines[4].split()[1]
        index = self.countIndex(reference)

        # print output
        outF = open(config.geneOutput, 'w')

        exons = []

        utr3 = True
        # print data
        for section in index:
            outF.write(str(config.geneType) + ",")
            if(section.name == "UTR3"):
                outF.write("Five_prime-UTR,1,")
            elif(section.name == "UTR5"):
                outF.write("Three_prime-UTR,1,")
            else:
                outF.write(section.name[:-1] + ",")
                outF.write(section.name[-1] + ",")
            i = section.length
            geneSegment = inputSeq[0:i]
            geneSegment = geneSegment.replace('-', '')
            outF.write( geneSegment + ','+"\n")
            inputSeq = inputSeq[i:]

            if section.name.startswith('Exon'):
                exons.append(geneSegment)

        outF.close()
        print("Validating the parsed coding sequences......")
        cds = ''.join(exons)
        cds = cds.upper()
        if self.validateCds(cds) or self.validateCds(cds[1:]) or self.validateCds(cds[2:]):
            print("Passed validation!  \n")
        else:
            print("Failed validation: stop codon(s) detected inter the CDS, please double check the input sequence.\n")


    def countIndex(self, reference):
        index = []
        for section in self.gene.data:
            count = 0
            i = section.length
            while i != 0:
                if reference[i] != '-':
                    i -= 1
                    count += 1
                else:
                    count += 1
            index.append(Section(count, section.name))

        return index

    def buildMap(self):
        self.map.update({"TCC": "S"})
        self.map.update({"TCA": "S"})
        self.map.update({"TCG": "S"})
        self.map.update({"TCT": "S"})

        self.map.update({"TTC": "F"})
        self.map.update({"TTT": "F"})
        self.map.update({"TTA": "L"})
        self.map.update({"TTG": "L"})

        self.map.update({"TAC": "Y"})
        self.map.update({"TAT": "Y"})
        self.map.update({"TAA": "*"})
        self.map.update({"TAG": "*"})

        self.map.update({"TGC": "C"})
        self.map.update({"TGT": "C"})
        self.map.update({"TGA": "*"})
        self.map.update({"TGG": "W"})

        self.map.update({"CTA": "L"})
        self.map.update({"CTC": "L"})
        self.map.update({"CTG": "L"})
        self.map.update({"CTT": "L"})

        self.map.update({"CCA": "P"})
        self.map.update({"CCC": "P"})
        self.map.update({"CCT": "P"})
        self.map.update({"CCG": "P"})

        self.map.update({"CAC": "H"})
        self.map.update({"CAT": "H"})
        self.map.update({"CAA": "Q"})
        self.map.update({"CAG": "Q"})

        self.map.update({"CGA": "R"})
        self.map.update({"CGC": "R"})
        self.map.update({"CGG": "R"})
        self.map.update({"CGT": "R"})

        self.map.update({"ATA": "I"})
        self.map.update({"ATC": "I"})
        self.map.update({"ATT": "I"})
        self.map.update({"ATG": "M"})

        self.map.update({"ACA": "T"})
        self.map.update({"ACC": "T"})
        self.map.update({"ACG": "T"})
        self.map.update({"ACT": "T"})

        self.map.update({"AAC": "N"})
        self.map.update({"AAT": "N"})
        self.map.update({"AAA": "K"})
        self.map.update({"AAG": "K"})

        self.map.update({"AGC": "S"})
        self.map.update({"AGT": "S"})
        self.map.update({"AGA": "R"})
        self.map.update({"AGG": "R"})

        self.map.update({"GTA": "V"})
        self.map.update({"GTC": "V"})
        self.map.update({"GTG": "V"})
        self.map.update({"GTT": "V"})

        self.map.update({"GCA": "A"})
        self.map.update({"GCC": "A"})
        self.map.update({"GCG": "A"})
        self.map.update({"GCT": "A"})

        self.map.update({"GAC": "D"})
        self.map.update({"GAT": "D"})
        self.map.update({"GAA": "E"})
        self.map.update({"GAG": "E"})

        self.map.update({"GGA": "G"})
        self.map.update({"GGC": "G"})
        self.map.update({"GGG": "G"})
        self.map.update({"GGT": "G"})

    def validateCds(self, cds):
        chunks = [cds[i:i + 3] for i in range(0, len(cds), 3)]
        for i in range(0, len(chunks)):
            try:
                # if last chunk doesn't contain three nts, don't translate it
                if i == len(chunks)-1 and len(chunks[i]) != 3:
                    return True
                aa = self.map[chunks[i]]
            except KeyError:
                print(chunks[i])
                return False
            if aa == "*" and i != len(chunks)-1:
                return False

        return True
