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


## Functional Analysis of Transcription Factors in Arabidopsis Methods and tools for RNA-seq-based co-expression network analysis
1. Quality control
   - FastQC: [Link](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
   - RSeQC: [Link](http://rseqc.sourceforge.net/)
   - QoRTs: [Link](http://hartleys.github.io/QoRTs/)

2. Read Mappers
   - Bowtie/Tophat/Tophat2: [Link](https://ccb.jhu.edu/software/tophat/index.shtml)
   - STAR: [Link](https://code.google.com/p/rna-star/)
   - HISAT: [Link](http://www.ccb.jhu.edu/software/hisat/index.shtml)
   - BWA: 
   - Kallisto: [Link](https://pachterlab.github.io/kallisto/about.html)
   - Salmon: [Link](http://combine-lab.github.io/salmon/)

3. Read counting tools
   - HTseq: [Link](http://www-huber.embl.de/HTSeq/doc/overview.html)
   - FeatureCounts: [Link](http://bioinf.wehi.edu.au/featureCounts/)
   - SpliceNet: [Link](http://jjwanglab.org/SpliceNet/)

4. Normalization
   - FPKM/RPKM: 
   - TPM: 
   - TMM: 
   - RAIDA: 
   - DEseq2: 

5. Correction for batch effects
   - Limma-remove BatchEffect: 
   - Svaseq: [Link](https://github.com/jtleek/svaseq)
   - Combat: [Link](http://www.bu.edu/jlab/wp-assets/ComBat/Abstract.html)

6. Co-expression module detection
   - WGCNA: [Link](https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/)
   - DiffCoEx: 
   - DICER: 
   - CoXpress: [Link](http://coxpress.sourceforge.net/)
   - DINGO: 
   - GSCNA: 
   - GSVD: 
   - HO-GSVD: [Link](https://github.com/aanchan/hogsvd-python/blob/master/README.md)
   - Biclustering: 

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
   - CoRegNet: 
   - cMonkey: 

9. Visualization
   - Cystoscape: [Link](http://www.cytoscape.org/)
   - BioLayout: [Link](http://www.biolayout.org/)

10. Co-expression databases
    - COXPRESdb: [Link](http://coxpresdb.jp/)
    - GeneFriends: [Link](http://www.genefriends.org/)
    - GeneMANIA: [Link](http://www.genemania.org/)
    - GENEVESTIGATOR: [Link](https://genevestigator.com/gv/)
    - GIANT: [Link](http://giant.princeton.edu)

## Proteomics

1. **Protein Dynamics**

   - [MakeMultimer](http://watcut.uwaterloo.ca/tools/makemultimer/)
     - Description: Calculates the coordinates of the missing subunits of PDB structures that represent multimeric molecules.

   - [Protein Interaction Server](http://pic.mbu.iisc.ernet.in)
     - Description: Recognizes various kinds of interactions within a protein or between proteins in a complex.

   - [The ConSurf Server](http://consurf.tau.ac.il/2016/)
     - Description: Identification of function regions of a protein.

   - [meta-PPISP](http://pipe.scs.fsu.edu/meta-ppisp.html)
     - Description: Binding site prediction - Protein.

   - [Protein Plus](https://proteins.plus/)
     - Description: Binding site prediction - Small Ligands.

   - [Hotspot Wizard](https://loschmidt.chemi.muni.cz/hotspotwizard/)
     - Description: Design of mutations and smart libraries in protein engineering.

   - [CAVER](http://www.caver.cz/index.php)
     - Description: Analysis and visualization of tunnels and channels in protein structures.

2. **Nucleic Acid Resources**

   - [Nucleic Acid Database](http://ndbserver.rutgers.edu)
     - Description: Experimentally-determined nucleic acids and complex assemblies.

   - [EM Databank](http://www.emdatabank.org)
     - Description: Unified data resource for 3-dimensional electron microscopy.

3. **Functional Annotation**

   - [DAVID](https://david.ncifcrf.gov/home.jsp)
     - Description: Provides a comprehensive set of functional annotation tools for investigators to understand the biological meaning behind a large list of genes.

4. **Biological Macromolecular Structures**

   - [PDB - Europe](https://www.ebi.ac.uk/pdbe/)
     - Description: European resource for the collection, organization, and dissemination of data on biological macromolecular structures.

## other misc links

**Gene Set Analysis:**
1. MSigDB - [Link](http://software.broadinstitute.org/gsea/msigdb/index.jsp)

**Motif Search Tools:**
2. Melina II - [Link](http://melina2.hgc.jp/public/index.html)
3. Prosite - [Link](https://prosite.expasy.org/prosite.html)

**Protein Sequence Analysis and Classification:**
4. Interpro - [Link](http://www.ebi.ac.uk/interpro/)

**Transcription Factor Binding Site Analysis:**
5. Pscan - [Link](http://159.149.160.88/pscan/)

**RNA Analysis:**
6. psRNATarget - [Link](http://plantgrn.noble.org/psRNATarget/)
7. ViennaRNA Web Services - [Link](http://rna.tbi.univie.ac.at)

**Plant-Specific Databases:**
8. P3DB - [Link](http://p3db.org/index.php)
9. PMRD: Plant microRNA Database - [Link](http://bioinformatics.cau.edu.cn/PMRD/)
10. Arabidopsis miRNA candidates - [Link](http://sundarlab.ucdavis.edu/mirna/search_candidates.html)
11. Trava - [Link](http://travadb.org)
12. Arabidopsis eFP Browser - [Link](http://bar.utoronto.ca/efp/cgi-bin/efpWeb.cgi)
13. Arabidopsis Hormone Database - [Link](http://ahd.cbi.pku.edu.cn)
14. AtPIN: Arabidopsis thaliana protein interaction network - [Link](https://atpin.bioinfoguy.net/cgi-bin/atpin.pl)

**miRNA Analysis:**
15. miARma-Seq - [Link](https://www.nature.com/articles/srep25749)
16. miRGen - [Link](http://carolina.imis.athena-innovation.gr/diana_tools/web/index.php?r=mirgenv3%2Findex)
17. TarBase - [Link](http://carolina.imis.athena-innovation.gr/diana_tools/web/index.php?r=tarbasev8%2Findex/)
18. miRNEST 2.0 - [Link](http://rhesus.amu.edu.pl/mirnest/copy/home.php)
19. pssRNAMiner - [Link](http://bioinfo3.noble.org/pssRNAMiner/)

**Metabolomics Analysis:**
20. MetaboAnalyst - [Link](http://www.metaboanalyst.ca)

**RNA-Seq Analysis:**
21. RNA-seqlopedia - [Link](https://rnaseq.uoregon.edu/#figure3.4)
22. RNA Seq Tutorial - [Link](https://github.com/griffithlab/rnaseq_tutorial)
23. CoExpress - [Link](http://sablab.net/coexpress.html)
24. GREIN: GEO RNA-seq Experiments Interactive Navigator - [Link](https://shiny.ilincs.org/grein/?gse=)

**Protein Analysis:**
25. UbPred - [Link](http://www.ubpred.org)
26. iPTMNet - [Link](https://research.bioinformatics.udel.edu/iptmnet/)
27. PaxDB - [Link](https://pax-db.org)
28. MicroRPM - [Link](http://microrpm.itps.ncku.edu.tw)
29. Interactome 3D - [Link](https://interactome3d.irbbarcelona.org/index.php)

**Other Tools:**
30. Java Heatchart - [Link](http://www.javaheatmap.com)
31. CyTrargetLinker - [Link](https://projects.bigcat.unimaas.nl/cytargetlinker/tutorial-1/)
32. Cyverse - [Link](https://wiki.cyverse.org/wiki/dashboard.action)
33. Ensembl Plants - [Link](https://plants.ensembl.org/info/website/ftp/index.html)
34. OmicsDB - [Link](https://www.omicsdi.org/home)
35. GSA: Genome Sequence Archive - [Link](http://bigd.big.ac.cn/gsa/)
36. PMirKB: Plant microRNA Knowledge Base - [Link](http://bis.zju.edu.cn/pmirkb/index.php)
37. RNA Tools - [Link 1](https://web.njit.edu/~wangj/rna/sequence.htm), [Link 2](http://rna.informatik.uni-freiburg.de/INFORNA/Input.jsp)
38. Expression Angler - [Link](http://bar.utoronto.ca/ntools/cgi-bin/ntools_expression_angler.cgi)
39. KAAS - KEGG Automatic Annotation Server - [Link](https://www.genome.jp/tools/kaas/)




