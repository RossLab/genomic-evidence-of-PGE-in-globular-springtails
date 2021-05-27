from sys import stdout
from sys import argv

stdout.write("""
##agp-version	2.0
# ORGANISM: Allacma fusca
# TAX_ID: 39272
# ASSEMBLY NAME: Afus1_asm01
# ASSEMBLY DATE: 09-November-2019
# GENOME CENTER: Edinburgh Genomics
# DESCRIPTION: Whole genome of Allacma fusca, decontaminated, reassembled and scaffolded using transcriptomic evidence using SCUBAT
""")

with open(argv[1], 'r') as ff:
    for line in ff:
        if line.startswith(">"):
            header = line.rstrip('\n')
            scf, contigs = header.split('[')
            scf = scf[1:-1] # do I want > in there?
            contigs = contigs.rstrip(']').split(',')
            part_number = 1
            scf_from = 1
            for i, contig in enumerate(contigs):
                contig = contig.split("'")[1]
                if i == 0:
                    strand_1 = contig.split('_')[6]
                    strand_relative_to_1 = "+"
                else:
                    if strand_1 == contig.split('_')[6]:
                        strand_relative_to_1 = "+"
                    else:
                        strand_relative_to_1 = "-"
                ctg = "_".join(contig.split('_')[0:2])
                contig_len = int(contig.split('_')[3])
                # print('contig')
                stdout.write("{}\t{}\t{}\t{}\tW\t{}\t1\t{}\t{}\n".format(scf, scf_from, scf_from + contig_len - 1, part_number, ctg, contig_len, strand_relative_to_1))
                scf_from = scf_from + contig_len
                part_number += 1
                if i + 1 < len(contigs):
                    # print('gap')
                    stdout.write("{}\t{}\t{}\t{}\tU\t100\tscaffold\tyes\talign_trnscpt\n".format(scf, scf_from, scf_from + 99, part_number))
                    part_number += 1
                    scf_from += 100
        #     seq = ''
        #
        # else:
        #     seq += line.rstrip('\n')
        #     break

