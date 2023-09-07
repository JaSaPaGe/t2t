import click as ck
import math
import os
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

@ck.command()
@ck.option(
    '--alignment-file', '-af', default='minimap/verkko_trio_aln_hap2.srt.paf',
    help='Sorted output of minimap2 alignment to CHM13 without centromeres')
@ck.option('--identity',  '-i', default=90)
@ck.option('--seq-len', '-s', default=10000)
def main(alignment_file, identity, seq_len):
    """Takes sorted output of alignment file. Use sort -k6,6 -k8,8n """
    chroms = [f'chr{i}' for i in range(1, 23)]
    chroms.append('chrX')
    chroms.append('chrY')
    chroms.append('chrM')

    with open(alignment_file) as f:
        for line in f:
            it = line.strip().split('\t')
            ctg, ctg_len, ctg_start, ctg_end = it[0], int(it[1]), int(it[2]), int(it[3])
            ch, ch_len, ch_start, ch_end = it[5], int(it[6]), int(it[7]), int(it[8])
            match_num = int(it[9])
            base_num = int(it[10])
            ident = round(match_num / base_num * 100, 2)
            seq_n = ch_end - ch_start + 1
            # Ignore hits with low identity and short length
            if ident < identity or seq_n < seq_len:
                continue
            print(ch, ch_len, ch_start, ch_end, ctg_len, ctg_start, ctg_end, identity, ctg, seq_len, it[4])


if __name__ == '__main__':
    main()
