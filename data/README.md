# Data used as part of the pipeline for the analysis of mitochondrial reads

## Simulated datasets

We use the [wgsim](https://github.com/lh3/wgsim) utility to simulate reads with desired parameters.

### nuclear_data

Data simulated with the followig wgsim commands:

Single end reads, length 100bp and 100,000 reads created each time.
No indels, no mutations.
Error rate varying: 0.0 0.1 0.5 1.0 2.0 5.0 1.0

## Reference mitochondria files

chrMT and shifted_chrMT are taken from the Revised Cambridge Reference Sequence of the human mitochondrial DNA ([rCRS](https://www.mitomap.org/MITOMAP/HumanMitoSeq)) - NCBI [entry](https://www.ncbi.nlm.nih.gov/nuccore/251831106). The breakpoint in the shifted version was artificially offset by **8000** bases to take into account the circularity of mitochondria.

## Numts files

Reference numts files as coordinates, stored in the [bed](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) format. If available, two versions are used: one with the coordinates of numts on the nucleus (\*.nuclear.bed), and one with the corresponding location of the numts on the mitochondria (\*.mitochondrial.bed).

### Calabrese et. al

HG19 numts reference coordinates as implemented in the UCSC Genome Browser.

> Calabrese, F. M., Simone, D., & Attimonelli, M. (2012). Primates and mouse NumtS in the UCSC Genome Browser. BMC bioinformatics, 13(4), S15.

### Dayama et. al

Numts reference coordinates as used in the [dinumt](https://github.com/mills-lab/dinumt) program.

> Dayama, G., Emery, S. B., Kidd, J. M., & Mills, R. E. (2014). The genomic landscape of polymorphic human nuclear mitochondrial insertions. Nucleic acids research, 42(20), 12640-12649.
