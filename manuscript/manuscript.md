### Single cell interaction mapping captures spatially defined architecture in the murine intestine

Jason T. Serviss1•, Nathanael Andrews1•, … random people … , Marco Gerling, Martin Enge1

1 Department of Oncology-Pathology Karolinska Institutet, Stockholm, Sweden.  

• Equal contribution  
Correspondence to: Martin Enge (martin.enge@ki.se)

#### Abstract

Advances in single-cell biology has enabled investigation of isolated single cells at an unprecedented scale and resolution. However, cells in multicellular organisms are largely defined by their spatial organisation within organ structures, which means that general methods for studying direct cell interaction on a large scale are needed. Here we propose a novel method, Cell Interaction by Multiplet Sequencing (CIM-Seq), that uses computational deconvolution of RNA-seq data from partially dissociated tissue to create cell interaction maps. We applied CIM-seq to murine small intestine, demonstrating that it recapitulates known cell interactions such as stem cell-paneth cell contacts. Furthermore, we discover ... Thus, CIM-Seq is a general method for cell interaction studies that can be used on cell types defined to an arbitrary resolution allowing identification of interacting sub-cell types or cell states.

#### Introduction

The lower gastrointestinal (GI) tract serves as an important element in the digestive system facilitating water and nutrient uptake, waste collection in addition to being a microbial boundary. The lower GI tract is further divided into the small and large intestine each having additional distinct sub-divisions, with the small intestine containing the duodenum, jejunum, ileum and the large intestine containing the cecum, colon, and rectum. These sub-divisions differ in their functionality, histology, gene expression profile and, in some cases, developmental origin (source). The micro-architecture of the small and large intestine is composed of columnar epithelial cells which form crypts, invaginations into the underlying lamina propria, at the base of which resident stem cells give rise to the multiple cell types which comprise the tissue. Mention e.g. stem - paneth interaction here? In the small intestine these cells are arranged in the form of villi - finger shaped structures that greatly increase the available surface area and thus the tissues nutrient absorption capability (figure 1x). Villi are comprised of multiple cell types including absorptive columnar cells, hormone-producing chromaffin cells, and mucus-secreting goblet cells. In contrast, the large intestine lacks villi although it does contain crypts, similar to the small intestine, where stem cells produce major cell types akin to those found in the small intestine. Recent single cell RNA sequencing (scRNA-seq) studies have provided an unparalleled glimpse into the plethora of cell types and cell states in the lower GI tract (source, source) indicating previously unobserved diversity in this tissue. Despite the accumulation of knowledge concerning the archtecture of the lower GI tract, the molecular mechanisms underlying the observed heterogeneity within these structures is largely unknown. Importantly, disease susceptibility, clinical features, and patient outcome vary dramatically between the different structures. Colonic crypts are the source of colorectal carcinomas although tumors differ between proximal and distal colon, as do patient outcomes (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2650032/, https://www.nature.com/articles/nrc2392, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3037515/). Regional ileitis and celiac disease as well as other human intestinal disorders are largely confined to particular segments. (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3210370/). Initiating cells and processes that trigger disease onset in other gastrointestinal diseases such as inflammatory bowel disease, ulcerative colitis, and Crohn’s disease remain unknown. As such, an increased understanding of the cellular heterogeneity, tissue architecture, and underlying molecular mechanisms of the lower GI tract is imperative to further our understanding of the diseases incurred in these tissues.

Single-cell mRNA-seq (scRNA-seq) methods can measure subtle changes in cell state, revealing specialised subpopulations within a cell type, while providing, in parallel, the underlying gene expression giving rise to such changes. However, functional characterisation of such subpopulations remains challenging. A promising approach to untangling cellular function is spatial transcriptomics (Crosetto, Bienko, and van Oudenaarden 2014, J. Lundeberg paper?), where gene expression measurements are complemented with spatial attributes. Cell-cell interaction is known to be important in establishment and maintenance of diverse organ structures (Sato and Clevers 2013) and, as such, these approaches allow inference of cellular function via their cellular interactions. Unfortunately, current spatial transcriptomic approaches suffer either from a lack of sensitivity, low resolution, are limited to a restricted panel of genes, or require complicated and laborious protocols, often requiring specialised equipment.

Here we describe Cellular Interactions via Multiplet Sequencing (CIM-seq) a novel single cell resolution high-throughput method for interrogating physical cell-cell interactions. Using fluorescence activated cell sorting (FACS) and *in silico* deconvolution, CIM-seq utilises single cells and multiplets, a commonly discarded bi-product of tissue disassociation, to provide a detailed map of the cellular connectome. We utilise CIM-seq to characterise cell interactions in the mouse small intestine identifying a previously known enrichment in small intestinal stem cell and paneth cell interactions. In addition, we examine cell-cell interactions along the span of the mouse colon establishing...

----

In addition to the villi, small intestine crypts are responsible for the regeneration of the cell types of the villi taking place every 3-5 days (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4883796/). Stem cells residing in the crypt give rise to the 

At least one paragraph about the architecture of small intestine and colon. Start with general funciton (nutrient and water uptake, microbial boundary), mention fast cell turnover (every 3 days). Mention known macro (pathological definitions; small intestine (duodenum, jejunum, ileum: note this differs in mouse); colon (cecum, colon, and rectum)) and micro architecture (crypts present in both SI and colon, villi - villi not present in colon) cell types (absorptive columnar cells; mucus-secreting goblet cells; and hormone- producing enteroendocrine/chromaffin cells) and their interactions (give stem-paneth as e.g., use this as an opportunity to highlight the importance of cell-cell interactions). 
The proximal and distal colon also differ developmentally, with the midgut forming the proximal region, and the hindgut forming the distal region. These regions also display distinct histology and function, with the proximal colon serving more to absorb water, and distal to store feces. Colorectal cancers also differ between these regions, as do patient outcomes (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2650032/). Importantly, colonic crypts are also the source of colorectal carcinomas (https://www.nature.com/articles/nrc2392). Also potentially mention lymphoid structures such as Peyer’s patches (small intestine) and isolated lymphoid follicles (colon) but probably not because of lack of blood in seq data. Start with small intestine since it is better characterized and lead into colon with transition mentioning that colon is not as well studied. Mention increased understanding of gastrointestinal diseases such as inflammatory bowel disease, ulcerative colitis, Crohn’s disease, intestinal cancer; include that initiating cells and processes that trigger disease remain unknown. "Certain human intestinal disorders, such as regional ileitis and celiac disease, are largely confined to particular segments." (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3210370/). "In another illustration of sub-mucosal heterogeneity, Hirschsprung disease, a congenital segmental aganglionosis that usually results from defective expansion or migration of neural crest-derived enteric neurons, affects the terminal colon far more often than proximal regions [36]." (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3210370/)

----

#### Results

In order to deconvolute the cell connectome, we use FACS to sort multiplets of incompletely dissociated cells along with single cells from the same sample. The single cells provide a set of blueprints of possible sub-cell types, and subsequent machine learning approaches allow us to infer the composition of the multiplets. Using this method we can measure the composition of a very large number of multiplets onto a lattice consisting of sub-cell types of arbitrary complexity.

The first step in single-cell RNA-seq methods is obtaining a suspension of dissociated single cells, usually by enzymatic digestion of the extracellular matrix holding the tissue together. Cells which are not fully dissociated (cell multiplets) are typically removed in the subsequent sorting step. In CIM-Seq, we instead use these multiplets to determine which cells were directly adjacent to each other in the intact tissue. CIM-Seq consists of three main parts (Fig 1a). In the first step, we perform partial dissociation of the target tissue, followed by cell sorting of singlets and multiplets into multiwell plates and conventional Smart-Seq2 library preperation. Second, we use automated feature selection followed by dimensionality reduction and graph-based clustering to construct a set of cell type templates from the singlet data. These cell clusters represent a blueprint of available cell types or cell states that we want to project the multiplets onto. In the third step, we employ computational deconvolution of multiplets into the cell types defined previously. Cell types that are frequently found together in multiplets can be inferred to interact in the intact tissue. 

As a dataset for testing our method, we obtained 131 singlet and 69 multiplet RNA-seq libraries by sorting cells from a fetal pancreatic rudiment and analyzing their gene expression using SmartSeq2. 

Creation of cell-type blueprint and performance testing on simulated data
We used dimensionality reduction followed by model-based clustering to obtain an unbiased set of thirteen expression profiles corresponding to different cell types or cell states. By mapping the expression values of known marker genes to each cell, we determined the cell-type identities of the clusters (Fig 2, Fig. S7). Interestingly, mesenchymal cells exhibited a large degree of heterogeneity, forming multiple discrete clusters representing distinct cell states (or sub-cell types) and which were readily identified by our unbiased clustering method (Fig 2a). Such diversity in mesenchymal cells is known to exist, but has so far not been thoroughly studied due to technical difficulties (Nombela-Arrieta, Ritz, and Silberstein 2011). 

To measure the performance of CIM-seq in a controlled setting, we first used simulated data where the cell type composition of each multiplet is known. As shown in Fig. S1-4, CIM-seq successfully identified the correct cell type composition of in silico multiplets even when substituting 80% of expression values with random noise. Similarly, multiplets generated in silico from random combinations of cell types derived from our dataset of fetal pancreatic cells, were correctly identified (Fig S5-6).

CIM-seq correctly identifies known interacting cell types
Next, we applied CIM-seq to the 69 cell multiplets sorted from the same pancreatic rudiment as the singlets used to construct the blueprint of cell types. To determine the approximate number of cells that make up each multiplet, we added a set of synthetic mRNA controls at a known concentration to each well. As shown in Fig. 1b, the fold difference in the fraction of synthetic mRNA indicated that most multiplets were comprised of 3-5 cells. 

CIM-seq indicated that acinar and ductal cells frequently interact, in line with current knowledge of pancreas development. Also, endocrine progenitor cells were strongly connected to endothelial cells, again corroborating prior knowledge that islets are highly vascularized. As expected due to the small number of multiplets investigated, not all cell types were represented in the 69 multiplets (indicated by star in legend in Fig 3a,b Fig 4a). The missing cell types were all rare also in the singlet data, indicating that bigger numbers of multiplets would be needed to reliably cover all cell types. 

Model residual analysis suggests specific features of centroacinar cells
The ductal cells most proximal to the pancreatic acinus are called centroacinar cells. These cells have generated considerable interest in cancer research, since they have been suggested as a frequent cell of origin for pancreatic ductal adenocarcinoma (PDAC) (Stanger et al. 2005). Centroacinar cells are relatively rare and very similar to ductal cells, and so did not produce a discrete cluster in our analysis. However, ductal-like interaction partners of acinar cells are expected to be enriched in centroacinar cells. Since our optimization strategy is missing a blueprint for these cells, centroacinar-specific genes will fit worse to the predicted model than the genes expressed in both duct and centroacinar cells.

We exploited this side-effect of the method and determined the expression level normalized residual between model fit and observed gene expression values, for multiplets that contained both acinar (cluster H1) and ductal (cluster I1) blueprints. As shown in figure 3e, the top genes include DLK1 which is known to be expressed in distal growing pancreatic ducts(Falix et al. 2012), and KRT8 expressed in lobular breast carcinoma(Moriya et al. 2006). Further studies will be required to determine whether genes with high residual are due to specific changes in duct growth during organ development or if they represent features of centroacinar cells.

A mesenchymal cell subset strongly interacts with endocrine progenitors
A specific subtype of mesenchymal cells is known to be required for differentiation of endocrine progenitor cells into hormone producing islet cells, but the identity of this subtype is not known. One of the most highly significant interactions identified by CIM-seq was between mesenchymal subtype L1 and endocrine progenitors, indicating that this subtype constitutes the cellular context for cells that undergo differentiation into pancreatic endocrine cells. Thus, our data suggests a specific identity of the subtype of mesenchymal cells that is required for the creation of endocrine islet cells.

We performed differential gene expression analysis using a version of the Kolmogorov-Smirnov statistic to determine which genes are uniquely expressed in each mesenchymal cluster (Fig 4c). Subtype L1, which interacts with endocrine progenitors, expresses diverse signalling molecules such as PEBP1, known to modulate the MAPK pathway, and focal adhesion molecules such as TENC1. Future experimental evidence is required to determine in what way L1-specific genes drive endocrine differentiation, and whether the identified signaling molecules can be used in therapy, for example by pushing duct cells to differentiate into insulin producing beta cells in diabetes type I patients.




Figure 1 The CIM-Seq method. a) Overview of the method. b) Fraction of reads from ERCC spike-in controls in singlets (left) and multiplets (right). A constant amount of non-human control RNA is spiked into each well when the plate is prepared. Relative abundance of synthetic spike-in RNA indicates number of cells per multiplet. Left axis: Number of cells, based on fold change over mean spike in fraction in singlets. Right axis: raw values of ERCC control as a fraction of total RNA. 
 


Figure 2 2d tSNE dimensionality reduction of singlets based on pairwise pearson correlation between the transcriptional profile of the 2000 genes most highly expressed in a cell. Each cell is represented by a point. a) Points are colored by cluster identity. Clusters were determined using mclust, based on pariwise distances in tSNE (2d projected) spaceb) color is by expression of marker genes indicated at top of respective panels.


Figure 3 CIM-seq uncovers known cell-cell interactions in the developing pancreas. a) Connections involving endocrine progenitor cells (Cluster B1).  Lines indicate connections between clusters (large colored points) of cells (small colored points). Line thickness indicates the number of connections. Stars in legend indicate cell clusters that were not found in the multiplet data. b) Statistical significance of connections between endocrine progenitors (B1) and every cluster (A1-M1). Statistical significance is derived from a poisson test, with lambda set to median number of connections to B1.  c) Connections involving acinar cells. Legend as in (a). d) Statistical significance of connections between acinar cell cluster (H1) to clusters (A1-M1). Legend as in (b). e) Model residual analysis of multiplets incorporating acinar cell cluster (H1) and ductal cell clusters (I1). Average normalized residual is plotted for all genes, outliers are labeled.




Figure 4 Mesenchymal subtype L1 is strongly connected to endocrine progenitors. a) Connections involving mesenchymal subtype H1.  Lines indicate connections between clusters (large colored points) of cells (small colored points). Line thickness indicates the number of connections. Stars in legend indicate cell clusters that were not found in the multiplet data. b) Statistical significance of connections to mesenchymal subtype H1 and cluster A1-M1. Statistical significance is derived from a poisson test, with lambda set to median number of connections to H1. c) Unique genes in each mesenchymal subcluster were extracted using K-S statistics.


#### Materials and Methods

**_Murine gut procurment and FACS sorting_**
Human Pancreas and Islet Procurement and FACS sorting. All studies involving human pancreas or islets were conducted in accordance with Stanford University Institutional Review Board guidelines. The de-identified human pancreatic rudiment, gestational week 14, was obtained from the University of Washington School of Medicine.

Isolated human islets were dissociated into single cells by enzymatic digestion using Liberase DL (Sigma). Prior to antibody staining, cells were incubated with blocking solution containing FACS buffer (2% v/v fetal bovine serum in PBS and goat IgG [Jackson Labs], 11.2 μg per million cells). LIVE/DEAD Fixable Aqua Dead Cell Dye (Life Technologies) was used as a viability marker. Cells were then stained with appropriate antibodies at 1:100 (v/v) final concentration. Anti human EpCAM-APC (Biolegend, 324208) was use for staining doublets. Cells were sorted on a special order 5-laser FACS Aria II (BD Biosciences) using a 100 m nozzle. Cell singlets were determined by using the FSC-H and FSC-A channels with singlets and multiplets sorted into separate plates.

Sorted single cells were collected directly into 96-well plates (Bio-Rad cat #: HSP9601) containing 4 µL of lysis buffer with dNTPs and ERCC spike-in controls.

**_scRNA-seq_**
Single-cell and multiplet RNA-Seq libraries were generated as described (Picelli et al. 2014). Briefly, single-cells collected in 96-well plates were lysed, followed by reverse transcription with template-switch using an LNA-modified template switch oligo to generate cDNA. After 21 cycles of pre-amplification, DNA was purified and analyzed on an automated Fragment Analyzer (Advanced Analytical). Each cell’s cDNA fragment profile was individually inspected and only wells with successful amplification products (concentration higher than 0.06 ng/ul) and with no detectable RNA degradation were selected for final library preparation. Tagmentation assays and barcoded sequencing libraries were prepared using Nextera XT kit (Illumina) according to the manufacturer's instructions. Barcoded libraries were pooled and subjected to 75 bp paired-end sequencing on the Illumina NextSeq 500 instrument.

Sequencing reads were trimmed, adapter sequences removed and the reads aligned to the hg19 reference assembly using STAR (Dobin et al. 2013) with default parameters. Duplicate reads were removed using picard (McKenna et al. 2010). Transcript counts were obtained using HTSeq (Anders, Pyl, and Huber 2014) and hg19 UCSC exon/transcript annotations. Transcript counts were normalized into log transformed counts per million (log2(counts * 1 000 000 / total_counts  + 1). Single cell profiles with the following features were deemed to be of poor quality and removed : 1) cells with less than 100.000 total number of valid counts on exonic regions. 2) cells with very low actin CPM. To determine a cutoff for actin CPM, we used the normal distribution with empirical mean and standard deviation from log2 transformed actin CPM, which is normally distributed in successful experiments. The cutoff was set to the 0.01 quantile (eg. the lower 0.01 % of the bell curve). In total 131 singlets (out of 180 analyzed wells, 72.8% success rate) and 69 multiplets (out of 95, 72.6% success rate) were retained for further analysis.

**_CIM-seq_**

CIM-seq is implemented in the R statistical programming language (source) and is freely available at (where) under a GPL-3 license. The CIM-seq method consists of two distinct stages; a preperatory stage, including feature selection, dimensionality reduction, and classification, and a deconvolution stage.

Preperatory stage  
CIM-seq permits the elements of the preperatory stage to be carried out using any appropriate methods selected by the end user. For the analysis of the data included here we choose to utilize the Seurat package version (source) for all 3 preperatory stage elements; i.e. feature selection, dimensionality reduction, and classification. Dimensionality reduction for visualization was performed using the t-SNE (source) (or UMAP) algorithm. Unsupervised classification of cell types and cell states was performed using a graph-based method. Breifly, feature selection was performed using *FindVariableGenes* function, selecting only those genes with a mean expression > 1 and dispersion > 1 and using default parameters otherwise. Principal component analysis (PCA) with the resulting features was performed followed by Jack straw testing to determine significant (p < XXX) principal components. Significant principal components were subsequently used as input to the graph-based classification algorithm via the *FindClusters* function. Specifically, the louvain community detection algorithm was utilized to detect cell types with 100 random starts. Classification was often performed in an iterative fashion with various settings of the resolution parameter with final classifications judged to be sufficient based on the partitioning of known cell types and marker genes. Feature selection for downstream deconvolution was carried out via the *FindAllMarkers* function using the AUC test and setting the minimum difference in the fraction of detection parameter to 0.4 and log fold change threshold to log(2).

Include number of features included per analysis here. Include Louvain resolution for each analysis here.

Deconvolution stage  
In the deconvolution stage we utilize the blueprint of the possible cell types/states in the tissue from the preperatory stage to determine the cellular composition of each of the multiplets. The deconvolution takes advantage of particle swarm optimization (PSO) where a candidate solution (swarm particle), consisting of a vector of fractions with one value per cell type, is optimized with respect to a cost over a number of iterations. PSO is optimal for solving this problem due to the fact that it makes few assumptions about the problem being optimized and can search a large space for candidate solutions (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4436220/). The cost function essentially works by synthesizing synthetic multiplets using the singlet gene expression data and subsequently generating a gene-wise probability of the observed multiplet expression value given the synthetic multiplet values. To accomplish this, a defined number of singlet sets are composed with one randomly selected singlet from each of the predefined classes present in each set. For each iteration each singlet set is multiplied with the corresponding fraction in the candidate solution. The sum for each gene is calculated in each of the sets thus generating one synthetic multiplet per set. The Poisson density for each gene in the observed multiplet is then calculated based on the synthetic multiplet values for the corresponding gene giving a vector of densities corresponding to the probability of the observed expression value given the synthetic multiplet values. Subsequently summing the log10 mean densities for each gene and taking the inverse of the sum generates the cost. The final result of the deconvolution procedure gives one vector of fractions per multiplet and a corresponding cost.

We utilized the SPSO2007 reference implementation of the PSO algorithm provided by the psoptim function in the pso (source) R package. In each of the deconvolutions we allow a maximum iteration of 100 with 500 swarm particles and 400 synthetic multiplets are provided to the cost function. Candidate solutions are [0, 1] constrained.

**_Data availability and Reproducibility statement_**




#### Discussion

Recently, several methods that allow analysis of highly multiplexed gene expression and spatial information have been developed. Multiplexed mRNA staining methods(Crosetto, Bienko, and van Oudenaarden 2014) generally provide very high spatial accuracy, allowing for cellular or even sub-cellular localization of transcripts, while sacrificing the number of genes that can be interrogated. Importantly, the constraint that the set of genes to measure is determined beforehand limits their use to testing known hypotheses rather than generating novel ones, although higher multiplexing allows a large set of hypotheses to be tested in parallel. Index based methods (Ståhl et al. 2016) are not limited to a pre-defined set of genes, but are instead limited in spatial accuracy by the size of the barcoded features precluding true single-cell accuracy. Also, the relatively low number of barcoded adapters limits the sensitivity of these methods to detect low-expressed genes.

CIM-seq solves many of the problems of previous methods. All transcripts are measured, similarly to the index based methods, but with the advantage that multiplets are investigated with no loss of sensitivity compared to conventional scRNA-seq. This is important, since it means that we can mine the multiplet data for interesting genes using differential gene expression methods (eg. see Fig 4c). Even subtle gene expression changes due to interaction, that do not result in a new cell type, can potentially be detected by mining for genes with a large residual in the deconvolution step (Fig 3e). Also, since CIM-seq relies on actual physical interaction of cells in intact tissue, it has single-cell spatial accuracy. The major drawback is that it only allows us to obtain information on direct interactions and cannot detect higher-order structures, which may be important for organ development by, for example, creating a concentration gradient of signalling molecules in the tissue.

When analyzing the fetal pancreas data, we find known interactions such as ductal/acinar cells interacting at the interphase between ducts and acini. As has been previously suggested, the mesenchymal cells include a variety of sub-cell types. By analyzing these mesenchymal subtypes as individual clusters, we find that only a particular subtype has strong connections to endocrine progenitors. It is known from cell transplantation experiments that the mesenchymal context is crucial for the differentiation of ductal progenitors into endocrine cells. However, the identity of the mesenchymal cells that provide this context has not been elucidated, and no molecular mechanism has been shown. 

Thus, we have designed and tested a novel method, CIM-Seq, for creating cell interaction maps from primary tissue. The method accurately predicts known cellular interactions, and suggests a specific subset of mesenchymal cells provide the cellular context required for endocrine cell differentiation. Mesenchymal cells have been previously used when re-establishing endocrine function in mice with type I diabetes (Figliuzzi et al. 2014). Identifying the molecular mechanism of how this is achieved would be an important step towards designing reconstitution protocol for diabetes-associated beta cell dysfunction.

CIM-seq is a general method capable of analyzing any solid tissue, and we expect it to be useful in a wide variety of scientific questions. This includes organ development and developmental diseases, decrease of fitness in aging which has been partly attributed to loss of organ integrity, and tumor-stromal interactions in malignancies. Thus, we anticipate that CIM-Seq will be broadly used to generate novel hypotheses about how specific cell interactions influence cell function and identity.

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



