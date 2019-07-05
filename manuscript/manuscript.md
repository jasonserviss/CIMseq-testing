<font size="5"><strong>Single cell interaction mapping captures spatially defined architecture in the murine intestine</strong></font>

Jason T. Serviss<sup>1•</sup>, Nathanael Andrews<sup>1•</sup>, …, Marco Gerling<sup>2</sup>, Martin Enge<sup>1</sup>

<sup>1</sup> Department of Oncology-Pathology, Karolinska Institutet, Stockholm, Solna, Sweden.  
<sup>2</sup> Department of Biosciences and Nutrition, Center of Innovative Medicine, Karolinska Institutet, Huddinge, Sweden.  
• Equal contribution  

Correspondence to: Martin Enge (martin.enge@ki.se)

### Abstract

Advances in single-cell biology has enabled investigation of isolated single cells at an unprecedented scale and resolution. However, cells in multicellular organisms are largely defined by their spatial organisation within organ structures, which means that general methods for studying direct cell interaction on a large scale are needed. Here we propose a novel method, Cell Interaction by Multiplet Sequencing (CIM-Seq), that uses computational deconvolution of RNA-seq data from partially dissociated tissue to create cell interaction maps. We applied CIM-seq to murine small intestine, demonstrating that it recapitulates known cell interactions such as stem cell-paneth cell contacts. Furthermore, we discover ... Thus, CIM-Seq is a general method for cell interaction studies that can be used on cell types defined to an arbitrary resolution allowing identification of interacting sub-cell types or cell states. (I think we really need to stress here or somewhere the benifits of our methods compared to others) INCLUDE LINK TO SOFTWARE

### Introduction

The lower gastrointestinal (GI) tract serves as an important element in the digestive system facilitating water and nutrient uptake, waste collection in addition to being a microbial boundary. The lower GI tract is further divided into the small and large intestine each having additional distinct sub-divisions, with the small intestine containing the duodenum, jejunum, ileum and the large intestine containing the cecum, colon, and rectum. These sub-divisions differ in their functionality, histology, gene expression profile and, in some cases, developmental origin (source). The micro-architecture of the small and large intestine is composed of columnar epithelial cells which form crypts, invaginations into the underlying lamina propria, at the base of which resident stem cells give rise to the multiple cell types which comprise the tissue. Mention e.g. stem - paneth interaction here? In the small intestine these cells are arranged in the form of villi - finger shaped structures that greatly increase the available surface area and thus the tissues nutrient absorption capability (figure 1x). Villi are comprised of multiple cell types including absorptive columnar cells, hormone-producing chromaffin cells, and mucus-secreting goblet cells. In contrast, the large intestine lacks villi although it does contain crypts, similar to the small intestine, where stem cells produce major cell types akin to those found in the small intestine. Recent single cell RNA sequencing (scRNA-seq) studies have provided an unparalleled glimpse into the plethora of cell types and cell states in the lower GI tract (source, source) indicating previously unobserved diversity in this tissue. Despite the accumulation of knowledge concerning the archtecture of the lower GI tract, the molecular mechanisms underlying the observed heterogeneity within these structures is largely unknown. Importantly, disease susceptibility, clinical features, and patient outcome vary dramatically between the different structures. Colonic crypts are the source of colorectal carcinomas although tumors differ between proximal and distal colon, as do patient outcomes (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2650032/, https://www.nature.com/articles/nrc2392, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3037515/). Regional ileitis and celiac disease as well as other human intestinal disorders are largely confined to particular segments. (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3210370/). Initiating cells and processes that trigger disease onset in other gastrointestinal diseases such as inflammatory bowel disease, ulcerative colitis, and Crohn’s disease remain unknown. As such, an increased understanding of the cellular heterogeneity, tissue architecture, and underlying molecular signatures of the lower GI tract is imperative to further our understanding of the diseases incurred in these tissues.

Single-cell mRNA-seq (scRNA-seq) methods can measure subtle changes in cell state, revealing specialised subpopulations within a cell type, while providing, in parallel, the underlying gene expression giving rise to such changes. However, functional characterisation of such subpopulations remains challenging. A promising approach to untangling cellular function is spatial transcriptomics (Crosetto, Bienko, and van Oudenaarden 2014, J. Lundeberg paper?), where gene expression measurements are complemented with spatial attributes. Cell-cell interaction is known to be important in establishment and maintenance of diverse organ structures (Sato and Clevers 2013) and, as such, delineation of these interactions allow inference of cellular function. Unfortunately, current spatial transcriptomic approaches suffer either from a lack of sensitivity, low resolution, are limited to a restricted panel of genes, or require complicated and laborious protocols, often requiring specialised equipment. Add additional comparisons to other spatial methods:

R.Satija,J.A.Farrell,D.Gennert,A.F.Schier,andA.Regev,“Spatialreconstructionofsingle-cellgeneexpressiondata,” Nature biotechnology, vol. 33, no. 5, p. 495, 2015.

K.Achim,J.-B.Pettit,L.R.Saraiva,D.Gavriouchkina,T.Larsson,D.Arendt,andJ.C.Marioni,“High-throughputspatial mapping of single-cell rna-seq data to tissue of origin,” Nature biotechnology, vol. 33, no. 5, p. 503, 2015.

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5321580/ (need some spatial "index"?)

https://www.biorxiv.org/content/biorxiv/early/2018/10/30/456350.full.pdf (need to make lots of assumptions)

Here we describe Cellular Interactions via Multiplet Sequencing (CIM-seq) a novel single cell resolution high-throughput method for unsupervised interrogation of physical cell-cell interactions in a tissue. Using fluorescence activated cell sorting (FACS) and *in silico* deconvolution, CIM-seq utilises single cells and multiplets, a commonly discarded bi-product of tissue disassociation, to provide a detailed map of the cellular connectome. We utilise CIM-seq to characterise cell interactions in the mouse small intestine identifying a previously known enrichment in small intestinal stem cell and paneth cell interactions. In addition, we examine cell-cell interactions along the span of the mouse colon establishing... Finally, we...

### Results

##### CIMseq: a method facilitating the interogation of tissue architecture

Single-cell RNA-seq methods begin by obtaining a suspension of dissociated single cells, usually by enzymatic digestion of the extracellular matrix. Cells which are not fully dissociated (multiplets) are typically removed in a subsequent sorting step. In CIM-Seq, we instead use these multiplets to determine which cells were directly adjacent to each other in the intact tissue. CIM-seq consists of three main parts (Fig. 1a). In the first step, we perform partial dissociation of the target tissue, followed by cell sorting of singlets and multiplets into multiwell plates and conventional Smart-Seq2 library preperation. Second, we use the sequenced and quantified singlet data to perform automated feature selection followed by graph-based clustering to construct a blueprint of cell types and states in the tissue. In the third step, we employ computational deconvolution to project the multiplets onto the previously defined cell types. Specifically, the deconvolution uses particle swarm optimization (PSO) to optimize a cost function that calculates the probability of the observed multiplet expression values given a set of synthetic multiplet values derived from the singlet samples. 

We first tested the propensity of previously disassociated single cells to re-associate in suspension which would result in connections which are not representative of real connections in the tissue. To do this, we first sorted single cells from disassociated mouse small intestine into the same tube. At time points 0, 30, and 120 minutes, we re-analyzed the cells using the previously defined singlets gating scheme finding that 0.04, 0.42, and 0.24% of the cells had re-associated as multiplets indicating a negligable rate of re-association after initial disassociation (Supplementary Fig. 1a). 

To measure the performance of CIM-seq in a controlled setting, we first used three distinct cell lines and sorted these as both singlets and multiplets of a known composition. The final dataset included 79 singlets and 83 multiplets (31, 27, and 25 multiplets consiting of 2, 3, and 4 cells, respectivley). Evaluation of the estimated cell number based on the percentage of ERCC spike-ins showed a clear average increase in cells contained in multiplets compared to singlets (Fig. 1b). The estimated cell number was also shown to correspond well to the number of sorted cells per well (Fig. 1b). Together these results indicate that the ERCC-based estimation of cell numbers accuratley reflect the number of cells contained in the ... Next, we examined cell specific marker genes and found them to be rarely coexpressed in singlets but frequently coexpressed in a subset of multiplets (Supplementary Fig. 1b). Unsupervised classification of the singlets reflected the known identity of the cells (Supplementary Fig. 1c). Gene exp heatmap here? During the ongoing deconvolution, the cost function estimates the probability bla bla bla (Fig. 1c). As the cost function is minimized during the ongoing optimization the swarm bla bla to find the lowest cost (Fig. 1d). CIM-seq mediated deconvolution revealed a high level of correspondence between the expected and detected number of connections (Supplementary Fig. 1d) as well as a mean misclassification rate () of less than or equal to 0.12 in all tested multiplet compositions (Fig. 1e and Supplementary Fig. 1e). These results demonstrate the ability of CIM-seq to sucessfully bla bla bla.


##### CIMseq accuratley identifies the cell-cell interactions in the small intestine stem cell niche

With the extensive amount of previous research investigating the architecture of the mammal small intestine, we next choose to evaluate the performance of CIMseq in the mouse small intestine epithelium. 

The cell types of the small intestinal epithelium all derive from a resident stem cell population at the bottom of the small intestinal crypts (ref). As differentiation takes place, cells migrate up the vilus and towards the intestinal lumen. Characterization of the small intestine stem cell niche has previously shown the interaction between stem cells and paneth cells, and paracrine Wnt signaling, to be crucial in facilitating stem cell maintainance. To interogate the cell-cell interactions of the small intestine epithelium XXX singlets and XXX multiplets were isolated from XXX C57BL/J mice and sequenced using Smart-Seq2. Clustering and classification of the resulting data showed cell types characteristic to the small intestine and in agreement with previous studies in mouse small intestine (Umap fig, gene exp?). The results indicate a clear enrichment between Lgr5+ stem cells and Lyz1+ paneth cells. In addition, multiple connections between the multiple differentiation stages, from stem cell to all differentiated cell types, are also detected. In summary, ...


* Mention the manual counting of cell numbers contained in multiplets. 
* Include detailed (sub-classification) analysis in colon
* Potentially differential expression results will want to fit in here as well

*Introduce CIM-seq results*

* Use complete analysis to highlight low number of false positives
* Zoom in on small intestine and show detected known connections
* Zoom in on colon and include the "no colon paneth" result
* How to fit in the "residuals" analysis? 
* End with indicating general differences in colon and small intestine architecture

*Introduce segmented colon analysis in bulk and with multiplets*

#### Discussion

Recently, several methods that allow analysis of highly multiplexed gene expression and spatial information have been developed. Multiplexed mRNA staining methods(Crosetto, Bienko, and van Oudenaarden 2014) generally provide very high spatial accuracy, allowing for cellular or even sub-cellular localization of transcripts, while sacrificing the number of genes that can be interrogated. Importantly, the constraint that the set of genes to measure is determined beforehand limits their use to testing known hypotheses rather than generating novel ones, although higher multiplexing allows a large set of hypotheses to be tested in parallel. Index based methods (Ståhl et al. 2016) are not limited to a pre-defined set of genes, but are instead limited in spatial accuracy by the size of the barcoded features precluding true single-cell accuracy. Also, the relatively low number of barcoded adapters limits the sensitivity of these methods to detect low-expressed genes.

CIM-seq solves many of the problems of previous methods. All transcripts are measured, similarly to the index based methods, but with the advantage that multiplets are investigated with no loss of sensitivity compared to conventional scRNA-seq. This is important, since it means that we can mine the multiplet data for interesting genes using differential gene expression methods (eg. see Fig 4c). Even subtle gene expression changes due to interaction, that do not result in a new cell type, can potentially be detected by mining for genes with a large residual in the deconvolution step (Fig 3e). Also, since CIM-seq relies on actual physical interaction of cells in intact tissue, it has single-cell spatial accuracy. The major drawback is that it only allows us to obtain information on direct interactions and cannot detect higher-order structures, which may be important for organ development by, for example, creating a concentration gradient of signalling molecules in the tissue.

When analyzing the fetal pancreas data, we find known interactions such as ductal/acinar cells interacting at the interphase between ducts and acini. As has been previously suggested, the mesenchymal cells include a variety of sub-cell types. By analyzing these mesenchymal subtypes as individual clusters, we find that only a particular subtype has strong connections to endocrine progenitors. It is known from cell transplantation experiments that the mesenchymal context is crucial for the differentiation of ductal progenitors into endocrine cells. However, the identity of the mesenchymal cells that provide this context has not been elucidated, and no molecular mechanism has been shown. 

Thus, we have designed and tested a novel method, CIM-Seq, for creating cell interaction maps from primary tissue. The method accurately predicts known cellular interactions, and suggests a specific subset of mesenchymal cells provide the cellular context required for endocrine cell differentiation. Mesenchymal cells have been previously used when re-establishing endocrine function in mice with type I diabetes (Figliuzzi et al. 2014). Identifying the molecular mechanism of how this is achieved would be an important step towards designing reconstitution protocol for diabetes-associated beta cell dysfunction.

CIM-seq is a general method capable of analyzing any solid tissue, and we expect it to be useful in a wide variety of scientific questions. This includes organ development and developmental diseases, decrease of fitness in aging which has been partly attributed to loss of organ integrity, and tumor-stromal interactions in malignancies. Thus, we anticipate that CIM-Seq will be broadly used to generate novel hypotheses about how specific cell interactions influence cell function and identity.

#### Materials and Methods

* Add info for cell lines, i.e. purchased from, verified, mycoplasma tested.

**_Mice_**  
All animal procedures were performed in accordance with Karolinska Institutet guidelines. Adult, age matched C57BL/J (origin) mice were used for extraction of colonic and small intestinal crypts.

**_Crypt isolation for CIM-Seq_**  
The intestinal tract, ranging from the stomach to the rectum, was removed from C57BL/J wild-type mice and kept on ice. Small intestine and colon were subsequently removed and remaining parts were discarded. Lumina of colon and small intestine were washed three times with Dulbecco’s Phosphate Buffered Saline (PBS) and adipose tissue connected to the exterior of the small intestine or colon was removed. Colon and small intestine were opened longitudinally, and any remaining mucous on covering the epithelial layer was gently rubbed off. Tissues were washed once in PBS before being cut into 0.5 - 1 mm long fragments. Colon fragments were immersed in 10mM EDTA-PBS and incubated for 105 minutes on ice, with shaking every 15 minutes. Small intestinal fragments were similarly immersed in 10mM EDTA-PBS on ice and shaken. After 15 minutes small intestinal fragments were allowed to settle at the bottom and supernatant was discarded and replaced with new 10mM EDTA-PBS and shaken vigorously. This procedure was repeated 1-2 times, with supernatant fractions investigated through light microscopy before discarding, in order to enrich for crypts. Small intestinal fragments were then incubated on ice for 45-60 minutes for a total of 105 minutes. Following EDTA-PBS treatment fractions were triturated 10-15 times. Colon fractions were strained through a 100 um filter while small intestinal fractions were strained through a 70 um filter. Fractions were centrifuged at 300g for 5 minutes and dissociated using TrypLE Express (Invitrogen) at 37°C for 15-20 minutes. Enzymatic dissociation was supervised using light microscopy during regular intervals to obtain an appropriate amount of single cells and multiplets. Dead cells were removed using Dead Cell Removal Kit (Miltenyi Biotec) according to the manufacturer's protocol before staining for FACS sorting was performed.

**_Cell sorting_**  
Cells were stained with FITC anti-CD326 and 7AAD (Sony Biotechnology) in order to monitor number of epithelial cells and select live cells. Additionally, Paneth cells were enriched using Brilliant Violent anti-CD24. Single cells were distinguished from multiplets using FSC-H and FSC-W gating and the gating scheme was visually confirmed via phase contrast microscopy (Supplementary Fig. XXX). Cells were sorted into 384- or 96-well plates containing hypotonic lysis buffer using a SH800S FACS sorter (Sony). Plates were sealed immediately after sorting with microseal F foil (Biorad) and centrifuged at 4000g at 4°C for 5 min before being frozen on dry ice and stored at -80°C.

**_scRNA-seq_**  
Single-cell and multiplet RNA-Seq libraries were generated as described (Picelli et al. 2014). Briefly, single-cells collected in 384-well plates were lysed, followed by reverse transcription with template-switch using an LNA-modified template switch oligo to generate cDNA. After 21 cycles of pre-amplification, DNA was purified and analyzed on an automated Fragment Analyzer (Advanced Analytical). Each cell’s cDNA fragment profile was individually inspected and only wells with successful amplification products (concentration higher than 0.06 ng/ul) and with no detectable RNA degradation were selected for final library preparation. Tagmentation assays and barcoded sequencing libraries were prepared using Nextera XT kit (Illumina) according to the manufacturer's instructions. Barcoded libraries were pooled and subjected to 75 bp paired-end sequencing on the Illumina NextSeq 500 instrument.

Sequencing reads were trimmed, adapter sequences removed and the reads aligned to the hg19 reference assembly using STAR (Dobin et al. 2013) with default parameters. Duplicate reads were removed using picard (McKenna et al. 2010). Transcript counts were obtained using HTSeq (Anders, Pyl, and Huber 2014) and hg19 UCSC exon/transcript annotations. Transcript counts were normalized into log transformed counts per million. Single cell profiles with the following features were deemed to be of poor quality and removed: 1) cells with less than a specified total number of valid counts in exonic regions (1.8x104 and 4x104 for the mouse gut and sorted multiplets dataset, respectively), 2) cells with low actin TPM, and 3) cells with high fractions of ERCC reads. To determine a cutoff for actin TPM, we used the normal distribution with empirical mean and standard deviation from log2 transformed actin TPM, which is normally distributed in successful experiments. The cutoff was set to the 0.01 quantile (eg. the lower 0.01 % of the bell curve). A similar strategy was used for fractions of ERCC reads, samples with a log2 transformed fraction of ERCC reads above the 0.99 % quantile were rejected. In total XXX singlets (out of XXX analyzed wells, XXX% success rate) and XXX multiplets (out of XXX, XXX% success rate) were retained for further analysis.

**_CIM-seq_**  
CIM-seq is implemented in the R statistical programming language (source) and is freely available at (where) under a GPL-3 license. The CIM-seq method consists of two distinct stages; a preparatory stage, including feature selection, dimensionality reduction, and classification, and a deconvolution stage.

_preparatory stage_   
CIM-seq permits the elements of the preparatory stage to be carried out using appropriate methods selected by the end user. For the analysis of the data included here we choose to utilize the Seurat package version 2.3.4 (source) for all three preparatory stage elements; i.e. feature selection, dimensionality reduction, and classification. Dimensionality reduction for visualization was performed using the UMAP algorithm. Unsupervised classification of cell types and cell states was performed using a graph-based method. Briefly, data was scaled and centered using the *ScaleData* function after which feature selection was performed using *FindVariableGenes* function, selecting only those genes with a mean expression > 1 and dispersion > 1 and otherwise using default parameters. Principal component analysis with the resulting features was performed followed by Jack straw testing to determine significant (p < 0.001) principal components. Significant principal components were subsequently used as input to the graph-based classification algorithm via the *FindClusters* function. Specifically, the louvain community detection algorithm was utilized to detect cell types with 100 random starts. Classification was often performed in an iterative fashion with various settings of the resolution parameter with final classifications judged to be sufficient based on the partitioning of known cell types and marker genes. Feature selection for downstream deconvolution was carried out via the *FindAllMarkers* function using the AUC test and setting the minimum difference in the fraction of detection parameter to 0.4 and log fold change threshold to log(2).

Include number of features included per analysis. Include Louvain resolution for each analysis. Include number of significant principal components per analysis.

_Deconvolution stage_  
In the deconvolution stage we utilize the blueprint of the possible cell types/states in the tissue from the preparatory stage to determine the cellular composition of each of the multiplets. The deconvolution takes advantage of particle swarm optimization (PSO) where [0, 1] constrained candidate solutions (swarm particles), consisting of a vector of fractions with one value per cell type, are optimized with respect to a cost over a number of iterations. PSO makes few assumptions about the problem being optimized, and does not require a differentiable optimization problem, while still being able to search a large space for candidate solutions (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4436220/). 

The CIMseq cost function is based on the probability, p(m | s,c), of observing the gene expression profile of a multiplet (m), given a candidate solution (s), and the singlet profiles from each class (c). The probability is determined empirically by creating a number of *in silico* multiplet profiles derived from gene expression values of one randomly selected singlet from each class (cv), multiplied by the candidate solution vector (s), and summed over all classes. Each of these gene expression point values are treated as individual poisson processes with λ = round(s * cv) with the joint point p-value (probability mass function, pmf) for a gene j in m (m<sub>j</sub>), over n randomly generated multiplets given by: 

<style>
.sumOver {
    position: relative;
    font-size: 2.5em;
    display: inline-block;
    line-height: .7em;
    vertical-align: middle;
}

.sumOver:before {
    font-size: 12px;
    display: block;
    position absolute;
    left: 0;
    top: 0;
    content: "n";
    width: 22px;
    text-align: center;
}

.sumOver:after {
    font-size: 12px;
    display: block;
    position absolute;
    left: 0;
    bottom: 0;
    content: "i = 1";
    width: 27px;
    text-align: center;
}
</style>

Equation (1)
<p align="center" style="font-size: 130%;">
	<span class="sumOver">&Sigma;</span>
	(<sup class="rest">&lambda;<sup>m<sub>j</sub></sup> e<sup>-&lambda;</sup></sup> &frasl; <sub><sup>m<sub>j</sub></sup></sub>!) / n
</p>​

The cost is defined as the sum of the -log10 of the probabilities, p(m | s,c), given by Equation (1). The final result of the deconvolution procedure gives one solution, a vector of fractions, per multiplet and a corresponding cost. 

Converting the solution vector of fractions into connections is achieved by first normalizing the elements of the solution vector corresponding to each cell type by the ERCC-estimated median relative RNA contribution for each cell type. The solution vector is then multiplied by the ERCC-estimated cell number for the multiplet, to take into account the number of approximated cells in the multiplet, and then rounded. Values in the resulting vector that are greater than 0 indicate the presence of the corresponding cell types in the multiplet and, simultaneously, a connection between those cell types.

CIMseq uses a modified SPSO2007 reference implementation of the PSO algorithm with two additions designed to improve performance when used in conjunction with CIMseq. These modifications include 1) increased user control over early stopping criteria allowing for early termination of the optimization and 2) acceptance of user supplied swarm particle starting positions. For 1) we allowed early termination in the case that the cost did not improve by 1 (costs tend to be on the scale of thousands) in 5 iterations for all deconvolutions. For 2), we supplied precalculated swarm positions for all possible cell class combinations of one and two and added random normally distributed noise with a mean of 0 and standard deviation of 1 / the number of cell classes. In each of the deconvolutions we allow a maximum of 100 iterations with a minimum of 400 swarm particles and 2000 synthetic multiplets.

Hypothesis testing
Determining the probability that a connection is significantly enriched is accomplished by first estimating the expected number of connections between all cell types. We assume that if connections are randomly distributed between cell types that the number of connections between two individual cell types would follow the relative abundance of those cell types. Therefore, we first estimate the relative abundance of each cell type from the deconvoluted multiplets and then subsequently calculate the expected number of edges based on that relative abundance and the total number of connections detected in the data. The probability of the obserrving the detected number of connections in light of the expected number of connections is then calculated using the .... bla bla expected to follow a ... distribution

Assumptions. Null.

**_Data availability and Reproducibility statement_**

#### Acknowledgements

#### References

Anders, Simon, Paul Theodor Pyl, and Wolfgang Huber. 2014. “HTSeq - A Python Framework to Work with High-Throughput Sequencing Data.” doi:10.1101/002824.  
Crosetto, Nicola, Magda Bienko, and Alexander van Oudenaarden. 2014. “Spatially Resolved Transcriptomics and beyond.” Nature Reviews. Genetics 16 (1): 57–66.  
Dobin, Alexander, Carrie A. Davis, Felix Schlesinger, Jorg Drenkow, Chris Zaleski, Sonali Jha, Philippe Batut, Mark Chaisson, and Thomas R. Gingeras. 2013. “STAR: Ultrafast Universal RNA-Seq Aligner.” Bioinformatics  29 (1): 15–21.  
Falix, Farah A., Daniël C. Aronson, Wouter H. Lamers, and Ingrid C. Gaemers. 2012. “Possible Roles of DLK1 in the Notch Pathway during Development and Disease.” Biochimica et Biophysica Acta 1822 (6): 988–95.  
Figliuzzi, Marina, Barbara Bonandrini, Sara Silvani, and Andrea Remuzzi. 2014. “Mesenchymal Stem Cells Help Pancreatic Islet Transplantation to Control Type 1 Diabetes.” World Journal of Stem Cells 6 (2): 163–72.  
McKenna, Aaron, Matthew Hanna, Eric Banks, Andrey Sivachenko, Kristian Cibulskis, Andrew Kernytsky, Kiran Garimella, et al. 2010. “The Genome Analysis Toolkit: A MapReduce Framework for Analyzing next-Generation DNA Sequencing Data.” Genome Research 20 (9): 1297–1303.  
MCLUST: Software for Model-Based Clustering, Density Estimation and Discriminant Analysis. 2002.  
Moriya, Takuya, Atsuko Kasajima, Kazuyuki Ishida, Yoshiyuki Kariya, Jun-Ichi Akahira, Mareyuki Endoh, Mika Watanabe, and Hironobu Sasano. 2006. “New Trends of Immunohistochemistry for Making Differential Diagnosis of Breast Lesions.” Medical Molecular Morphology 39 (1): 8–13.  
Nombela-Arrieta, César, Jerome Ritz, and Leslie E. Silberstein. 2011. “The Elusive Nature and Function of Mesenchymal Stem Cells.” Nature Reviews. Molecular Cell Biology 12 (2): 126–31.  
Picelli, Simone, Omid R. Faridani, Asa K. Björklund, Gösta Winberg, Sven Sagasser, and Rickard Sandberg. 2014. “Full-Length RNA-Seq from Single Cells Using Smart-seq2.” Nature Protocols 9 (1): 171–81.  
Sato, Toshiro, and Hans Clevers. 2013. “Growing Self-Organizing Mini-Guts from a Single Intestinal Stem Cell: Mechanism and Applications.” Science 340 (6137): 1190–94.  
Ståhl, Patrik L., Fredrik Salmén, Sanja Vickovic, Anna Lundmark, José Fernández Navarro, Jens Magnusson, Stefania Giacomello, et al. 2016. “Visualization and Analysis of Gene Expression in Tissue Sections by Spatial Transcriptomics.” Science 353 (6294): 78–82.  
Stanger, Ben Z., Bangyan Stiles, Gregory Y. Lauwers, Nabeel Bardeesy, Michael Mendoza, Ying Wang, Amy Greenwood, et al. 2005. “Pten Constrains Centroacinar Cell Expansion and Malignant Transformation in the Pancreas.” Cancer Cell 8 (3): 185–95.  


Supplementary Materials

 
Method Testing: Synthetic dataset
Here we utilize random data generated from the negative binominal distribution to simulate scRNAseq single cells. The dataset includes 1000 singlet expression profiles with 2000 genes each and 10 different cell types.
 

Supplementary Figure 1. t-SNE results in 2 dimensional space from the synthetic dataset with each color representing a synthetic cell type.
 
Singlet cells were subjected to t-SNE to provide a cell-state blueprint to the deconvolution process. The results (Supplementary Fig. S1) indicate all cell types group together with no overlap between the groups thus providing a “well behaved” dataset as a starting point for testing the CIM-seq method. Single cells were then combined using the mean of each gene expression value for all possible combinations of 2, 3, or 4 cells. 15 of these multiplets were randomly chosen to use in the testing process.
 

Supplementary Figure 2. Connections resulting from the deconvolution of the synthetic datasets multiplets. The connections are overlaid onto the previous results from t-SNE. The plot includes each individual cell (small colored points) and the mean for each classification group (large colored points). The number of detected connections between each cell type classification is represented by the weight of the connection edges.

The detected connections from the deconvolution process were then visualized in a network plot (Supplementary Fig. S2) and the detected connections were compared to the real connections (Supplementary Fig. S3). The results indicate that the CIM-seq method was successfully able to deconvolute all multiplets and detect the connections present in the multiplets.


Supplementary Figure 3. Testing results from the synthetic dataset. Real connections, as determined from the singlets used to construct each multiplet, are shown on the x axis. Connections detected by CIM-seq are shown on the y axis. The number of connections is represented by the color. 

Finally, we tested the CIM-seq method by adding increasing amounts of randomly chosen values to the gene expression values of the multiplets consisting of two, three, and four cells. 10 multiplets were randomly chosen for each fraction of noise tested. The deconvolution was subsequently run for each multiplet with the results indicating that CIM-seq successfully detects all connections even as noise levels approach 80% (Supplementary Fig S4).   

 
Supplementary Figure 4. CIM-seq performance on the synthetic dataset with increasing amounts of random noise. The x axis shows the fraction of random noise and the y axis shows the percentage of identified connections.The lower and upper hinges correspond to the 25th and 75th percentiles.The upper whisker extends from the hinge to the largest value no further than 1.5 * interquartile range from the hinge. The lower whisker extends from the hinge to the smallest value at most 1.5 * interquartile range of the hinge. Data beyond the end of the whiskers are outlying points and are plotted individually. Please note that some of the box plots appear as a line due to the fact that the result was the same for all 10 replicates.  

 
 
Method Testing: Pancreas dataset
Next we desired to test the CIM-seq method using conditions that more accurately mirror a real scenario. For this we utilized the single cells from the fetal pancreas dataset. Single cells were subjected to dimensionality reduction and classification of the resulting groups of cells (Fig. 2a). Multiplets were then generated using the mean gene expression values for each classified group in combinations of 2, 3, or 4 cells after which 30 of these multiplets were chosen randomly for testing.
The resulting connection network was plotted and examined for the expected connections (Supplementary Fig. S5). Testing results indicated that CIM-seq successfully identified all connections present in the multiplets (Supplementary Fig. S6).
 

Supplementary Figure 5. Connections resulting from the deconvolution of the pancreas datasets in silico generated multiplets. The connections are overlaid onto the previous results from t-SNE. The plot includes each individual cell (small colored points) and the mean for each classification group (large colored points). The number of detected connections between each cell type classification is represented by the weight of the connection edges.
 

Supplementary Figure 6. Testing results from the pancreas datasets in silico generated multiplets. Real connections, as determined from the singlets used to construct each multiplet, are shown on the x axis. Connections detected by CIM-seq are shown on the y axis. The number of connections is represented by the color. 















Class specific gene expression in fetal pancreas.


Supplementary Figure 7. A heatmap showing expression of genes detected to be specifically expressed in all classified cell types. Colors represent the sum of quantile based Kolmogorov–Smirnov test statistics for comparisons between each individual cell type and all other cell types.
 
Supplementary Figure 8. A heatmap showing expression of genes detected to be specifically expressed in THY1 expressing mesenchymal cell types. Colors represent the sum of quantile based Kolmogorov–Smirnov test statistics for comparisons between each individual cell type and all other cell types.


Fetal pancreas markers

Supplementary Figure 9. Expression of insulin and the Thy1-cell surface antigen in multiplets and single cells from the fetal pancreas dataset. Due to the fact that, in the developing pancreas, INS is specifically expressed in pancreatic endocrine cells and THY1 in mesenchymal cells, co-expression of these markers is only observed in the multiplets. 



