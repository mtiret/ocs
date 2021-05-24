# Optimal Contribution Selection

## How to run

Rscript simulate.R _ngen_ _alpha_ _beta_ _gamma_ _qtnfile_ _mate_ _shuffle_

where _ngen_ is the number of generations, _alpha_ is the weight of coancestry vs genetic gain, _beta_ and _gamma_ are transformation parameters of the realized genomic relationship matrix, _qtnfile_ is the name of the file with allelic effects, _mate_ a boolean indicating whether apply mate allocation, and _shuffle_ a boolean indicating whether considering the shuffled file. 

The folder where the program runs should also contain the genotype (in.gen), the genetic map (in.map.small), and the population simulator (meiosis, source provided). 

## Example input files

in.gen: genotypes of each SNPs (0/1, diploid). Each line is a locus, and each column is one haploid strand. 

in.map.small: a genetic map, with the starting position equal to 0 (in Morgan).

in.qtn: allelic effects of each locus (it is the effect of 1-encoded alleles). 

## Population simulator

**meiosis** is a program that simulates the next generation of a population, knowing the pedigree (i.e., the mating plan), the genotypes of the parents, and the genetic map. 

Compilation command: 
g++ main.cpp -L. -lconsole_option -o meiosis -static
