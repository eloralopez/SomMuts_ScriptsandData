#!/usr/bin/env python


#####
#DESCRIPTION: take in VCF, filters out the sites that have the same genotype in all branches, outputs genoytpes for the variable sites

import sys


input = open(sys.argv[1],'r') #example: /Users/eloralopez/Documents/SomaticMutations/OfuAug/AH09_Aug2011_cohort_20180106_merged_filtered_nomissing_nosym_noblanks.vcf
#another example: /Users/eloralopez/Documents/SomaticMutations/gnomex/OfuAug2013_75_cohort_merged_filtered.vcf
output = sys.argv[2] #example: /Users/eloralopez/Documents/SomaticMutations/OfuAug/AH09datatable.txt
#another example: /Users/eloralopez/Documents/SomaticMutations/OfuAug/AH88datatable.txt

with open(output, "w") as out_file:
    #header = ["chrom.pos", "ref", "alt", "AH88_1", "AH88_10", "AH88_11", "AH88_12", "AH88_13", "AH88_14", "AH88_15", "AH88_16", "AH88_17", "AH88_2", "AH88_3","AH88_4", "AH88_5", "AH88_6", "AH88_7", "AH88_8", "AH88_9"]
    header = ["chrom.pos", "ref", "alt", "AH09_10A", "AH09_10B", "AH09_11A", "AH09_11B", "AH09_1A","AH09_1B","AH09_2B","AH09_2C","AH09_3A","AH09_3C","AH09_4A","AH09_4C","AH09_5A","AH09_5B","AH09_6A","AH09_6C","AH09_7A","AH09_7B","AH09_8B","AH09_8C","AH09_9A","AH09_9C"]
    header_output = "\t".join(header)
    out_file.write(header_output + '\n')
        
    for line in input:

        if line.startswith('#'): #ignores all the header lines

            continue
        if not line.startswith('contig'): #ignores all those that are symbiont (NOT CORAL) contigs
            continue

        else:

            line = line.strip()

            items = line.split('\t')

            contig = items[0]

            position = items[1]
        
            concatenated = contig + "." + position

            ref_allele = items[3]


            alt_allele = items[4]

            genotypes = items[9:]

            geno_list=[]
            
            genoplusdepthlist=[]

            for genotype in genotypes:

                genos = genotype.split(':')

                alleles = genos[0].split('/')

                alleledepths = genos[1].split(',') #depth per allele at locus
                
                #FROM GATK: AD is the unfiltered allele depth, i.e. the number of reads that support each of the reported alleles. All reads at the position (including reads that did not pass the variant caller’s filters) are included in this number, except reads that were considered uninformative. Reads are considered uninformative when they do not provide enough statistical evidence to support one allele over another.
                
                refdepth = alleledepths[0] 
                altdepth = alleledepths[1]
                
                totaldepth_atlocus = genos[2] #totaldepth at locus

                #FROM GATK: DP is the filtered depth, at the sample level. This gives you the number of filtered reads that support each of the reported alleles. You can check the variant caller’s documentation to see which filters are applied by default. Only reads that passed the variant caller’s filters are included in this number. However, unlike the AD calculation, uninformative reads are included in DP.
                
                if int(alleles[0]) == 0:

                    a1 = ref_allele

    #
                else:

                    a1 = alt_allele

                if int(alleles[1]) == 0:

                    a2 = ref_allele

                else:

                    a2 = alt_allele
                genostring = a1 + '/' + a2
                genoPLUSdepth = genostring + "," + totaldepth_atlocus + "," + refdepth + ","  + altdepth
                # print(genoPLUSdepth)
                genoplusdepthlist.append(genoPLUSdepth)
                genoPLUSdepth = '\t'.join(genoplusdepthlist)
                # print(genoPLUSdepth)
                geno_list.append(genostring)
                genostring = '\t'.join(geno_list)
                # print(genostring)
                # genoplusdepthlist.append(genoPLUSdepth)
                
    #
            #genoPLUSdepth = '\t'.join(genoplusdepthlist)

            
            if all(x==geno_list[0] for x in geno_list):
                continue #ignores all sites that are not variable within-colony 
            else:
            
                outlist = [concatenated, ref_allele, alt_allele, genoPLUSdepth]
                
                outstring = '\t'.join(outlist)

                out_file.write(outstring+'\n')

input.close()

out_file.close()