# Mitochondrial reads analysis pipeline

This pipeline was developped to investigate the effect numts have on mitochondrial DNA NGS reads. It is implemented on the HAMSTR webpage, available through Imperial College London's network. Three main aspects were investigated to understand how they could lead to confounding when mapping mitochondrial reads.

## Aligners

Different aligners produce different results, so we decided to compare a few commonly used programs. We implemented 6 aligners in this pipeline: [BWA](http://bio-bwa.sourceforge.net/), [SOAP3-dp](https://github.com/aquaskyline/SOAP3-dp), [NovoAlign](http://www.novocraft.com/products/novoalign/), [NextGenMap](https://github.com/Cibiv/NextGenMap), [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), and [SMALT](http://www.sanger.ac.uk/science/tools/smalt-0).

## Remapping reads mapping to numts to the mitochondria

In order to investigate how discarding reads mapping to numts impacts variants, we built the pipeline as a 2 steps mapping process: reads mapping to numts regions in the nucleus are not taken into account in the first step, but remapped to the mitochondria in the second.

## Mitochondria circularity

Since mitochondria are circular molecules, mapping mitochondrial reads requires specific dispositions. In this pipeline, we map mitochondrial reads to two mitochondrial genomes: one is an unaltered reference mitochondrion (taken from Revised Cambridge Reference Sequence of the human mitochondrial DNA ([rCRS](https://www.mitomap.org/MITOMAP/HumanMitoSeq)), the other is the same reference but shifted by 8000 bases. Variants are then called on the two files using [mitoCaller](https://lgsun.irp.nia.nih.gov/hsgu/software/mitoAnalyzer/mitoAnalyzer.htm) and the resulting outputs combined.
