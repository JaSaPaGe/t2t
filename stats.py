import click as ck
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def load_lens(filename):
    res = {}
    with open(filename) as f:
        for line in f:
            it = line.strip().split()
            chr_id = it[0]
            length = int(it[1])
            res[chr_id] = length
    return res

def count_contigs(filename):
    res = {}
    with open(filename) as f:
        sequences = SeqIO.parse(f, 'fasta')
        for record in sequences:
            if not record.id.startswith('chr'):
                continue
            i = 0
            seq = record.seq
            j = seq.find('N', i)
            ctgs = 1
            while i < len(seq) and j != -1:
                ctgs += 1
                while seq[j] == 'N':
                    j += 1
                i = j
                j = seq.find('N', i)
            res[record.id] = ctgs
    return res

def lengths(filename):
    with open(filename) as f:
        sequences = SeqIO.parse(f, 'fasta')
        for record in sequences:
            print(record.id, len(record.seq))

def load_qv(filename):
    res = {}
    with open(filename) as f:
        for line in f:
            it = line.strip().split('\t')
            res[it[0]] = float(it[3])
    return res

@ck.command()
def main():
    lengths('data/assembly_ksa002/ksa002.asm.hap1.centro.fa')
    return
    ksa_len = load_lens('ksa001v0.3.0.len')
    chm_len = load_lens('chm13.len')
    qv = load_qv('merqury_v0.3.0.ksa001v0.3.0.qv')
    chroms = [f'chr{i}' for i in range(1, 23)]
    chroms.append('chrX')
    chroms.append('chrM')
    total_ksa = 0
    total_chm = 0
    total_ctgs = 0
    for chrom in chroms:
        # print(f'{chrom} & 1 & {chm_len[chrom]:,} & 1 & {ksa_len[chrom]:,} & {qv[chrom]:.2f} \\\\')
        print(f'|{chrom} | 1 | {chm_len[chrom]:,} | 1 | {ksa_len[chrom]:,}|')
        total_ksa += ksa_len[chrom]
        total_chm += chm_len[chrom]
    total_chm += chm_len['chrY']
    # print(f'chrY & 1 & {chm_len["chrY"]:,} & 1 & 0 & - \\\\')
    # print(f'Total & 25 & {total_chm:,} & 24 & {total_ksa:,} & 47.84 \\\\')
    print(f'|chrY | 1 & {chm_len["chrY"]:,} | 0 | 0 |')
    print(f'|Total | 25 | {total_chm:,} | 24 | {total_ksa:,}|')

if __name__ == '__main__':
    main()
