#!/usr/bin/env python

# Parse multialleic sites in an input vcf to seperate records for each of the respective alternate alleles.
# INFO column is simply copied.
# Genotypes are adjusted such that each non-alt allele is reset to 0 unless it is the alt allele in question
# Unlike previous versions this script makes no assumptions about ploidy
# Validate to work fine on test data

import gzip
import re
import sys


def parse_site(line):
    """
    Function for parsing vcf records [site info only]
    """

    i = re.split("\t", line)
    alleles = re.split(",", i[4])
    if len(alleles) == 1:
        return [line]
    return ["\t".join(i[:4] + [alleles[i2]] + i[5:]) for i2 in range(len(alleles))]


def get_alleles(line):
    """
    extract the allels (5th column of vcf file)
    returns a list of allels
    """

    # cut off is set to not split complete (very long) line
    cutoff = 250
    fields = line[:250].split("\t")
    # if not enough field are found, try to split a larger string
    while len(fields) <= 6:
        cutoff = cutoff + 250
        fields = line[:cutoff].split("\t")

    return fields[4].split(",")


def parse(line):
    """
    Function for parsing vcf records
    """
    line = re.sub(r"\./\.:[^\t\n]*", r"./.", line)
    alleles = get_alleles(line)
    if len(alleles) == 1:
        return [line]
    i = line.split("\t")

    gti = i[8].split(":").index("GT")  # Position of GT in format field
    i[-1] = i[-1].strip()

    # genereed the  fixed fields
    outlist = [
        ["\t".join(i[:4] + [alleles[i2]] + i[5:9])] for i2 in range(len(alleles))
    ]

    for genotype_field in i[9:]:
        genotype_field = genotype_field.split(":")
        genotypes = re.split(r"/|\|", genotype_field[gti])
        for i2 in range(len(alleles)):
            if "." in genotypes:  # Ignore sites with missing genotype data
                None
            else:
                a = str(i2 + 1)
                genotype_field[gti] = "/".join(
                    ["1" if allele == a else "0" for allele in genotypes]
                )
            outlist[i2].append(":".join(genotype_field))

    return [r + "\n" for r in ["\t".join(out) for out in outlist]]


def main():
    handle = sys.argv[1]
    if handle == "-":
        handle = sys.stdin
    elif handle.endswith(".gz"):
        handle = gzip.open(handle, "rt", encoding="utf8")
    else:
        handle = open(handle, "rt")

    # Load file
    line = handle.readline()

    # write header
    while line[0] == "#":
        sys.stdout.write(line)
        line = handle.readline()

    # If site only vcf
    # print(line)
    if len(re.split("\t", line)) == 8:
        while line != "":
            for i in parse_site(line):
                sys.stdout.write(i)
                None
            line = handle.readline()

    else:
        while line != "":
            for i in parse(line):
                sys.stdout.write(i)
                None
            line = handle.readline()


if __name__ == "__main__":
    main()
