import subprocess
import config


def blast(input_file, gene):
    subprocess.run(["bin/blastn", "-query", input_file, "-subject", gene, "-out",
                    config.blast_output, "-outfmt", "3", "-num_descriptions", "1", "-num_alignments", "0"])


def get_match_gene():
    with open(config.blast_output) as f:
        content = f.readlines()
    size = len(content)
    loop = 0
    myList = []
    while loop != size:
        if content[loop].startswith("Sequences producing"):
            loop = loop +2
            line = content[loop]
            line = line.strip()
            geneName = line.split()[0]
            myList.append(geneName.split(":")[1])
        else:
            loop = loop +1
    return myList
