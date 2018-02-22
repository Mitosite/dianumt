# Data used in the dianumt pipeline for analysis of mitochondrial reads, as implented on the HAMSTR webpage


## Reference mitochondria files

chrMT and shifted_chrMT are taken from the Revised Cambridge Reference Sequence of the human mitochondrial DNA ([rCRS](https://www.mitomap.org/MITOMAP/HumanMitoSeq)) - NCBI [entry](https://www.ncbi.nlm.nih.gov/nuccore/251831106). The breakpoint in the shifted version was artificially offset by **8285** bases to take into account the circularity of mitochondria.

## Numts files

Reference numts files as coordinates, stored in the [bed](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) format.

### Calabrese et. al

HG19 numts reference coordinates as implemented in the UCSC Genome Browser.

> Calabrese, F. M., Simone, D., & Attimonelli, M. (2012). Primates and mouse NumtS in the UCSC Genome Browser. BMC bioinformatics, 13(4), S15.

### Dayama et. al

Numts reference coordinates as used in the [dinumt](https://github.com/mills-lab/dinumt) program.

> Dayama, G., Emery, S. B., Kidd, J. M., & Mills, R. E. (2014). The genomic landscape of polymorphic human nuclear mitochondrial insertions. Nucleic acids research, 42(20), 12640-12649.
