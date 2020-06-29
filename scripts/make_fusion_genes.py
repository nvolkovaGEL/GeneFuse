"""
Generate a user-defined fusion file 

Input file should be a text file containing 
1. gene name and 
2. transcript name (optional)

For example:
gene1   transcript1
gene2
gene3   transcript3
"""
__author__ = "Kai"
__date__ = "20/05/2020"

import argparse


def read_genelist(inputf):
    with open(inputf, "r") as fh:
        for line in fh:
            yield line.rstrip("\n").split()


def make_fusion_gene(gene, fw, refflat):
    # no transcript specified --> use the longest transcript
    if len(gene) == 1:
        transcripts = {}
        with open(refflat, "r") as fh:
            for line in fh:
                cur_gene, transcript, chrom, strand, start, end, _, _, _, exonstart, exonend = line.rstrip("\n").split("\t")
                if cur_gene != gene[0]:
                    continue
                transcripts[transcript] = (chrom, strand, start, end, exonstart, exonend)
        if transcripts == {}:
            raise ValueError('This gene symbol cannot be found in refFlat.txt:' + gene[0])
        transcript = get_longest_transcript(transcripts.keys(), refflat)
        chrom, strand, start, end, exonstart, exonend  = transcripts[transcript]

    # use user-specified transcript
    elif len(gene) == 2:
        found = 0
        with open(refflat, "r") as fh:
            for line in fh:
                cur_gene, transcript, chrom, strand, start, end, _, _, _, exonstart, exonend = line.rstrip("\n").split("\t")
                if cur_gene == gene[0] and transcript == gene[1]:
                    found = 1
                    break
            if found == 0:
                 raise ValueError('Wrong gene symbol or transcript provided:' + gene[0] + ', ' + gene[1])

    # write to a file
    header = '>' + str(gene[0]) + '_' + str(transcript) + ',' + str(chrom) + ':' + str(start) + '-' +  str(end) + '\n'
    fw.write(header)
    exons = list(zip(exonstart.split(","), exonend.split(",")))[:-1]
    if strand == '-':
        exons = exons[::-1]
    for index, each_exon in enumerate(exons, start=1):
        fw.write(str(index) + ',' + str(each_exon[0]) + ',' + str(each_exon[1]) + '\n')
    fw.write("\n")

def get_longest_transcript(transcripts, refflat):
    longest_length = 0
    longest_transcript = ""
    with open(refflat, "r") as fh:
        for line in fh:
            line = line.strip().split()
            if line[1] in transcripts:
                length = int(line[5]) - int(line[4])
                if length > longest_length:
                    longest_length = length
                    longest_transcript = line[1]
    return longest_transcript


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="Input filename")
    parser.add_argument("-r", "--refflat", required=True, help="Path to the refFlat.txt file, need to be downloaded from UCSC in advance")
    parser.add_argument("-o", "--output", required=True, help="The output filename")
    args = parser.parse_args()
    
    with open(args.output, "w") as fw:
        for gene in read_genelist(args.input):
            make_fusion_gene(gene, fw, args.refflat)
