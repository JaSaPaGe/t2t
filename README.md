# Assembly of genomes with Hifiasm and Masurca
## Example of hifiasm assembly with HiFi and HIC data
```
hifiasm -o genome.asm -t64 --h1 hic-reads1.fastq.gz  --h2 hic-reads2.fastq.gz hifi.fastq.gz
```
This command will produce the assemblies for haplotypes
two files:
* genome.asm.hic.hap1.p_ctg.gfa
* genome.asm.hic.hap2.p_ctg.gfa

And merged genome will be saved in:
* genome.asm.hic.p_ctg.gfa

The *.gfa files can be converted to FASTA with the following command:
```awk '/^S/{print ">"$2;print $3}' genome.asm.hic.p_ctg.gfa > genome.asm.fa```

## Chromosome scaffolding with [MaSuRCA](https://github.com/alekseyzimin/masurca)

```chromosome_scaffolder.sh -t 64 -r chm13v2.0.fa -q genome.asm.fa -s hifi.fastq.gz -ch 50```

The scaffolded genome will be saved in *.reconciled.fa file