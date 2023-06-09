# Database-links
This word document contains useful collected databases

# Plant Transcription Factor Databases and Tools

This repository contains a curated list of plant transcription factor databases and tools. Below are the categorized links:

## Plant (Arabidopsis) Transcription Factors

- [RARTF](http://rarge.gsc.riken.jp/rartf/)

- [AGRIS](http://arabidopsis.med.ohio-state.edu/AtTFDB/)

- [DATF](http://datf.cbi.pku.edu.cn/)

- [PlnTFDB](http://plntfdb.bio.uni-potsdam.de/v2.0/index.php?sp_id=ATH) - A part of a plant transcription factor database. Data of other plants are also stored.

- [Molbiol Tools](https://molbiol-tools.ca/Transcriptional_factors.htm)

## Conserved Domain Search

- [InterProScan](http://www.ebi.ac.uk/Tools/InterProScan/) - Perform against various Conserved Domain databases and provides sophisticated graphical output for known motifs.

- [MEME](http://meme.sdsc.edu/meme/intro.html) - Finding known or unknown Conserved Domains among a set of proteins for discovering unknown motifs.

- [SALAD database](http://salad.dna.affrc.go.jp/salad/en/) - For known and unknown motifs.

## Homology Search

- [TAIR BLAST](http://www.arabidopsis.org/Blast/index.jsp) - For Arabidopsis.

- [NCBI BLAST](http://blast.ncbi.nlm.nih.gov/Blast.cgi) - For multispecies search.

## Prediction of Subcellular Localization

- [SUBAII](http://www.plantenergy.uwa.edu.au/suba2/) - Provides hydropathy plots of all Arabidopsis proteins. Experimental data are also stored.

- SubLoc, TargetP, WoLF PSORT

- [LocDB](https://www.rostlab.org/services/locDB/index.php) - Protein Localization Database for Human and Arabidopsis.

## Protein-Protein Interaction

- [Arabidopsis predicted interactome](http://www.arabidopsis.org/portals/proteome/proteinInteract.jsp)

- [EBI IntAct](http://www.ebi.ac.uk/intact/site/index.jsf) - Stores continuously updated Protein-Protein Interaction (PPI) information of all organisms based on literature curation.

- [AtPID](http://atpid.biosino.org/index.php) - A search facility with graphical output against a predicted and literature-curated Arabidopsis PPI data set.

## Small RNAs

- [ASRP](http://asrp.cgrb.oregonstate.edu/db/) - Includes data of miRNA, siRNA, and ta-siRNA.

## Repository of Microarray Data

- [NCBI GEO](http://www.ncbi.nlm.nih.gov/geo/)

- [EBI ArrayExpress](http://www.ebi.ac.uk/microarray-as/ae/)

- [NASCArrays](http://affymetrix.arabidopsis.info/narrays/experimentbrowse.pl)

- [ATTED-II](http://atted.jp/)

- [Genevestigator](https://www.genevestigator.com/gv/index.jsp)

- [BAR eFP browser](http://bbc.botany.utoronto.ca/efp/cgi-bin/efpWeb.cgi)

## Finding Novel cis-elements

- [TAIR motif analysis](http://www.arabidopsis.org/tools/bulk/motiffinder/index.jsp)

- [AGRIS ATCISDB](http://arabidopsis.med.ohio-state.edu/AtcisDB/)

## GO Categorization

- [TAIR GO](http://www.arabidopsis.org/tools/bulk/go/index.jsp) - Annotation search.

- [agriGO](http://bioinfo.cau.edu.cn/agriGO/index.php) - GO Analysis Toolkit and Database (for Agricultural Community).


## Tools/Methods:
1. Quality control
   - FastQC: [Link](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
   - RSeQC: [Link](http://rseqc.sourceforge.net/)
   - QoRTs: [Link](http://hartleys.github.io/QoRTs/)

2. Read Mappers
   - Bowtie/Tophat/Tophat2: [Link](https://ccb.jhu.edu/software/tophat/index.shtml)
   - STAR: [Link](https://code.google.com/p/rna-star/)
   - HISAT: [Link](http://www.ccb.jhu.edu/software/hisat/index.shtml)
   - BWA: No link provided
   - Kallisto: [Link](https://pachterlab.github.io/kallisto/about.html)
   - Salmon: [Link](http://combine-lab.github.io/salmon/)

3. Read counting tools
   - HTseq: [Link](http://www-huber.embl.de/HTSeq/doc/overview.html)
   - FeatureCounts: [Link](http://bioinf.wehi.edu.au/featureCounts/)
   - SpliceNet: [Link](http://jjwanglab.org/SpliceNet/)

4. Normalization
   - FPKM/RPKM: No link provided
   - TPM: No link provided
   - TMM: No link provided
   - RAIDA: No link provided
   - DEseq2: No link provided

5. Correction for batch effects
   - Limma-remove BatchEffect: No link provided
   - Svaseq: [Link](https://github.com/jtleek/svaseq)
   - Combat: [Link](http://www.bu.edu/jlab/wp-assets/ComBat/Abstract.html)

6. Co-expression module detection
   - WGCNA: [Link](https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/)
   - DiffCoEx: No link provided
   - DICER: No link provided
   - CoXpress: [Link](http://coxpress.sourceforge.net/)
   - DINGO: No link provided
   - GSCNA: No link provided
   - GSVD: No link provided
   - HO-GSVD: [Link](https://github.com/aanchan/hogsvd-python/blob/master/README.md)
   - Biclustering: No link provided

7. Functional enrichment
   - DAVID: [Link](https://david.ncifcrf.gov/)
   - PANTHER: [Link](http://pantherdb.org/)
   - g:Profiler: [Link](http://biit.cs.ut.ee/gprofiler/)
   - ClusterProfiler: [Link](https://github.com/GuangchuangYu/clusterProfiler/blob/master/vignettes/clusterProfiler.Rmd)
   - Enrichr: [Link](http://amp.pharm.mssm.edu/Enrichr/)
   - ToppGene: [Link](https://toppgene.cchmc.org/)

8. Regulatory network inference
   - ARACNE: [Link](http://califano.c2b2.columbia.edu/aracne)
   - Genie3: [Link](https://bioconductor.org/packages/release/bioc/html/GENIE3.html)
   - CoRegNet: No link provided
   - cMonkey: No link provided

9. Visualization
   - Cystoscape: [Link](http://www.cytoscape.org/)
   - BioLayout: [Link](http://www.biolayout.org/)

10. Co-expression databases
    - COXPRESdb: [Link](http://coxpresdb.jp/)
    - GeneFriends: [Link](http://www.genefriends.org/)
    - GeneMANIA: [Link](http://www.genemania.org/)
    - GENEVESTIGATOR: [Link](https://genevestigator.com/gv/)
    - GIANT: [Link](http://giant.princeton.edu)




