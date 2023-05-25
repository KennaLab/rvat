#!/usr/bin/env python

# Parse multialleic sites in an input vcf to seperate records for each of the respective alternate alleles.
# INFO column is simply copied.
# Genotypes are adjusted such that each non-alt allele is reset to 0 unless it is the alt allele in question
# Unlike previous versions this script makes no assumptions about ploidy
# Validate to work fine on test data

import sys, re, gzip
handle=sys.argv[1]
if handle=="-":
	handle=sys.stdin
elif handle[-3:]==".gz":
	handle=gzip.open(handle,'r')
else:
	handle=open(handle,'r')


def parse_site(line):
	# Function for parsing vcf records [site info only]
	i=re.split("\t",line)
	alleles=re.split(",",i[4])
	if len(alleles)==1:
		return([line])
	else:
		lines=[]
		for i2 in range(len(alleles)):
			a=str(i2+1)
			out="\t".join(i[:4]+[alleles[i2]]+i[5:])
			lines.append(out)
		return(lines)


def parse(line):
	# Function for parsing vcf records
	line=re.sub("\./\.:[^\t\n]*","./.",line)
	i=re.split("\t",line)
	alleles=re.split(",",i[4])
	if len(alleles)==1:
		return([line])
	else:
		lines=[]
		gti=re.split(":",i[8]).index("GT")	# Position of GT in format field
		i[-1]=i[-1].strip()
		for i2 in range(len(alleles)):
			a=str(i2+1)
			out="\t".join(i[:4]+[alleles[i2]]+i[5:9])
			for i3 in i[9:]:				
				i3=re.split(":",i3)
				if "." in re.split("/|\|",i3[gti]):	# Ignore sites with missing genotype data
					None
				else:
					i3[gti]="/".join(["1" if allele==a else "0" for allele in re.split("/|\|",i3[gti])])
				out+="\t"+":".join(i3)
			lines.append(out+"\n")
		return(lines)



# Load file
line=handle.readline()

# write header
while line[0]=="#":
	sys.stdout.write(line)
	line=handle.readline()

# If site only vcf
if len(re.split("\t",line))==8:
	while True:
		if line=="":
			break
		for i in parse_site(line):
			sys.stdout.write(i)
			None
		line=handle.readline()

# If genotypes present
else:
	while True:
		if line=="":
			break
		for i in parse(line):
			sys.stdout.write(i)
			None
		line=handle.readline()

