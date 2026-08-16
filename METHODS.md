# Methods

## Genome dataset and quality control

Publicly available *Microcystis* genome assemblies were retrieved from the NCBI GenBank and RefSeq databases. The initial collection was dereplicated with dRep to reduce redundancy among highly similar assemblies. Genome quality was assessed with CheckM2, and assemblies with estimated completeness below 97% or contamination above 5% were excluded. The final dataset comprised 140 genomes. Retained assemblies were reannotated consistently with Prokka v1.14.6 to generate standardized annotations for orthology, pangenome and gene-content analyses. Geographic coordinates and sampling information were compiled from NCBI BioSample records and, where necessary, the corresponding original publications.

## Phylogenomic reconstruction and tree comparison

Phylogenomic relationships were reconstructed independently from single-copy orthologues and whole-genome SNPs. Single-copy orthologues were identified across the retained genomes with OrthoFinder. Their aligned sequences were used to infer a maximum-likelihood phylogeny with IQ-TREE v3.0. The best-fitting substitution model was selected with ModelFinder, and branch support was evaluated with 2,000 ultrafast bootstrap replicates.

For the SNP-based analysis, genome sequences were compared with the *Microcystis aeruginosa* NIES-843 reference genome (GCA_000010625.1), and high-quality SNPs were identified with Snippy v4.6.0. Indels and non-biallelic variants were excluded. A maximum-likelihood phylogeny was inferred from the whole-genome SNP alignment with IQ-TREE v3.0, again using ModelFinder and 2,000 ultrafast bootstrap replicates. Homologous recombination was identified with Gubbins v3.4.3, and a recombination-filtered SNP tree was reconstructed to assess whether the major topology was sensitive to recombinant regions.

Tree comparisons were restricted to shared genome tips. Robinson-Foulds distances were calculated with phangorn and expressed relative to the maximum possible RF distance. Cophenetic phylogenetic distance matrices were obtained with `cophenetic.phylo` in ape and compared using Pearson Mantel tests with 9,999 permutations. Mutual Clustering Information, Clustering Information Distance and Shared Phylogenetic Information were additionally calculated with TreeDist as supplementary measures of tree agreement.

## Population structure, ANI and geographic distribution

Population structure was evaluated from the recombination-filtered SNP data using two clustering approaches. Hierarchical Bayesian clustering was performed with fastBAPS using the optimized symmetric prior. An independent clustering analysis was performed with `snapclust` in adegenet. Candidate values from K = 1 to 10 were evaluated with `snapclust.choose.k`, and the preferred partition was selected using the Akaike information criterion.

Pairwise average nucleotide identity was calculated with FastANI and summarized as an all-by-all similarity matrix. Conventional ANI thresholds, including 95%, were examined to determine how threshold-based partitions compared with the higher-order population structure. Because fastBAPS and snapclust both supported five major clusters, these were retained as the conservative lineage framework C1-C5 for subsequent analyses rather than interpreted as formal species assignments.

Geographic structure was evaluated for genomes with valid lineage assignments and sampling coordinates. Sampling locations were grouped into six continent-level regions: Asia, Europe, North America, South America, Africa and Oceania. Associations between lineage membership and continent were tested using Fisher's exact test with Monte Carlo simulation (9,999 replicates). Differences among lineages in latitude, absolute latitude and longitude were evaluated using Kruskal-Wallis tests.

## Genome-wide diversity, differentiation and neutrality statistics

Population-genetic statistics were estimated with PopGenome v2.7.5. Nucleotide diversity (π), the fixation index (FST) and absolute nucleotide divergence (Dxy) were calculated in overlapping 50-kb windows with a 12.5-kb step across the genome. Genome-wide lineage and lineage-pair values were summarized from the corresponding window-level estimates. Differences in nucleotide diversity among lineages were evaluated using Kruskal-Wallis tests.

Tajima's D and Fu and Li's F were calculated for each lineage in overlapping 10-kb windows with a 2.5-kb step using PopGenome. Departure of window-level neutrality statistics from zero was evaluated with one-sample t-tests. Because these statistics can reflect multiple demographic and selective processes, they were used to characterize heterogeneity in population-genetic background rather than to assign a unique demographic history to each lineage.

## Multidimensional differentiation continuum

Pairwise differentiation among C1-C5 was summarized across five complementary dimensions: mean FST, mean Dxy, ANI distance, gene-content distance and pairwise gene flow. ANI distance was oriented so that larger values represented greater genome dissimilarity. Gene flow was reversed so that larger values represented stronger resistance to exchange. All five variables therefore increased in the direction of greater differentiation.

Each metric was standardized as a z score across the ten lineage pairs. An equal-weight composite score was calculated as the arithmetic mean of the five standardized variables. The same five variables were independently summarized by principal component analysis after centering and scaling. The direction of PC1 was reversed when necessary so that larger PC1 values corresponded to greater differentiation. The standardized PC1 score and equal-weight score were then averaged to produce a consensus continuum score, which was rescaled from 0 to 1 for visualization and ranking. The score was used as a relative summary of multidimensional lineage differentiation and was not interpreted as a probability of speciation or as a formal taxonomic scale.

To examine whether the gene-flow component simply recapitulated sequence-divergence information, gene-flow resistance was compared separately with FST, Dxy and ANI distance using Spearman matrix correlations. Significance was assessed with a label-permutation procedure in which the row and column labels of one symmetric matrix were permuted together. Each comparison used 9,999 permutations, and two-sided empirical P values were calculated from the absolute observed and permuted correlations.

## Homologous recombination and gene flow

Homologous recombination was inferred with Gubbins v3.4.3 from the whole-genome alignment. Recombinant intervals and the taxa associated with each event were extracted from the Gubbins output and linked to C1-C5 lineage assignments. Overlapping intervals were reduced with IRanges before affected sequence lengths were calculated.

For each focal genome, a recombinant region contributed to within-lineage exchange when the corresponding event included at least one additional genome from the same lineage. It contributed to outside-lineage exchange when the event included at least one genome from another lineage. For events involving both categories, the region could contribute to both summaries. Unique bases assigned to each category were divided by genome length and expressed as percentages. Lineage-level summaries were obtained from the genome-level values.

Pairwise shared-recombination matrices were constructed by calculating, for each genome, the proportion of sequence associated with recombinant events shared with each target lineage. Reciprocal lineage comparisons were averaged to obtain a symmetric lineage-by-lineage matrix. Resistance to gene flow was defined as decreasing shared recombination, with larger values representing less inferred exchange. This measure was used as a relative connectivity metric.

Genome-level estimates of ρ/θ, the relative rate of recombination to mutation, and r/m, the relative contribution of recombination and mutation to introduced substitutions, were summarized by lineage. Differences among C1-C5 were evaluated using Kruskal-Wallis tests followed by Dunn pairwise comparisons. Dunn-test P values were adjusted using the Benjamini-Hochberg procedure. The recombination-filtered alignment was also used to reconstruct a SNP phylogeny and assess the sensitivity of major phylogenomic relationships to recombinant regions.

## Pangenome, gene associations and horizontal gene transfer

The *Microcystis* pangenome was reconstructed with Panaroo v1.5.2 using the standardized Prokka GFF3 annotations. Gene families were classified by prevalence as core genes (present in 99-100% of genomes), soft-core genes (95-99%), shell genes (15-95%) or cloud genes (<15%).

Non-random associations among flexible gene families were evaluated with Coinfinder v1.1 using the Panaroo gene presence/absence matrix and the single-copy core-gene phylogeny. Significant associations and dissociations were identified using a Bonferroni-adjusted threshold of 0.05. These analyses were used to characterize organization within the flexible genome and were not used to define the differentiation continuum.

Putative horizontally transferred genes were identified with HGTector v2.0b using DIAMOND similarity searches against a locally compiled NCBI RefSeq protein database containing bacterial, archaeal, viral, fungal and protistan sequences. The database was compiled in September 2025. Because inferred HGT counts can scale with genome size, the relationship between raw HGT counts and genome size was first evaluated with Spearman rank correlation. HGT frequency was then standardized as the number of inferred HGT genes per megabase of genome sequence. Differences in normalized HGT frequency among C1-C5 were evaluated using the Kruskal-Wallis test.

## Environmental data and lineage-level environmental differences

Contemporary environmental conditions were described using bioclimatic, solar-radiation and UV-B variables. BIO1-BIO19 and 12 monthly solar-radiation layers were obtained from WorldClim v2.1 at 2.5-arc-min spatial resolution. Six UV-B variables were obtained from the glUV dataset. Environmental values were extracted at each genome sampling coordinate using terra.

Differences among C1-C5 in individual environmental variables were initially evaluated with Kruskal-Wallis tests. Where pairwise comparisons were required, Dunn tests with Bonferroni correction were used. Phylogenetic signal in environmental variables was evaluated using Pagel's λ with `fitContinuous` in geiger and Blomberg's K with `phylosignal` in picante; significance was assessed using likelihood-ratio or randomization tests, respectively.

## Geographic and environmental associations with genomic differentiation

For integrated landscape-genomic analyses, variables were divided into three environmental modules: BIO1-BIO19, monthly solar radiation and UV-B. Variables within each module were standardized and analysed separately by principal component analysis. The first two axes from each module were retained, giving six environmental axes in total. Pairwise environmental distance was calculated as Euclidean distance in this six-dimensional PCA space.

Three complementary genomic-distance matrices were constructed after matching genomes across metadata, ANI, recombination-filtered SNP and Panaroo datasets. ANI distance was calculated as 1 - ANI/100. SNP p-distance was calculated as the proportion of differing nucleotides among valid aligned sites, excluding gaps and ambiguous N positions. Gene-content distance was calculated as binary Jaccard distance from the Panaroo presence/absence matrix. Geographic distance was calculated as great-circle distance in kilometres from sampling coordinates using `distGeo` in geosphere.

Relationships among genomic, geographic and environmental distances were evaluated with Mantel and partial Mantel tests in vegan using Pearson correlations and 9,999 permutations. Partial Mantel analyses tested environmental distance while controlling for geographic distance and, conversely, geographic distance while controlling for environmental distance. Analyses were performed independently for ANI, SNP and gene-content distances.

Multivariate analyses were used to evaluate the environmental and spatial components of genomic variation. Spatial predictors comprised standardized longitude and latitude, their squared terms and their interaction. Marginal effects of the six environmental PCA axes and five spatial predictors were assessed with `adonis2` using 9,999 permutations and `by = "margin"`. Distance-based redundancy analysis was performed with `capscale`, fitting the six environmental axes while conditioning on the spatial terms. Overall model and constrained-axis significance were assessed with 9,999 permutations. Variation partitioning with `varpart` separated adjusted pure environmental, pure spatial, shared environment-space and unexplained fractions.

To identify individual climatic variables contributing to the multivariate BIO signal, candidate variables were selected from the five largest absolute loadings on BIO_PC1 and BIO_PC2. Unique candidate variables were then tested separately with marginal `adonis2` models containing the same spatial polynomial terms. These analyses were performed independently for ANI, SNP and gene-content distances with 9,999 permutations.

## Genotype-environment association analyses

Genotype-environment associations were evaluated using LFMM2 and partial redundancy analysis. Highly correlated bioclimatic predictors were removed using an absolute Pearson correlation threshold of |r| > 0.8, leaving BIO1, BIO2, BIO8, BIO12, BIO15 and BIO19. SNPs were restricted to biallelic sites and filtered to remove loci with more than 10% missing data, minor allele frequency below 0.05 or zero variance.

LFMM2 analyses were performed with LEA. Latent population structure was evaluated with `snmf` over K = 1-25, and the preferred K was selected using cross-entropy. Missing genotypes were imputed by the modal genotype before LFMM2 fitting. Associations between SNPs and each retained environmental predictor were tested with `lfmm2` and `lfmm2.test`. P values were adjusted using the Benjamini-Hochberg false-discovery-rate procedure, and loci with q < 0.05 for at least one predictor were retained as candidate environment-associated SNPs.

Partial redundancy analysis provided an independent multivariate assessment. The first two SNP principal components were included as conditioning variables to account for broad population structure. Significance of the overall constrained model, constrained axes and environmental terms was assessed by permutation with `anova.cca`. SNPs with loadings more than three standard deviations from the mean on significant constrained axes were recorded as partial-RDA candidates. Partial RDA was used as complementary evidence and was not used to restrict the LFMM2 candidate set entered into Gradient Forest.

## Gradient Forest modelling and future genomic offset

Gradient Forest models were fitted with the six retained bioclimatic predictors and the LFMM2 candidate environment-associated SNPs. Candidate identifiers were matched directly to the imputed SNP matrix, so all final LFMM2 candidates, rather than only loci also identified by partial RDA, were used as Gradient Forest response variables.

Models were fitted with the R package gradientForest using 500 trees per response variable and a correlation threshold of 0.5. The maximum tree level was set to max[1, floor(log2(0.368n/2))], where n was the number of matched genomes. Predictor importance and cumulative turnover functions were obtained from the fitted model. Contemporary environmental conditions were transformed through these functions to create a multivariate genomic-turnover space.

Future climate conditions were obtained from WorldClim v2.1 for the MPI-ESM1-2-HR global climate model under SSP2-4.5 for 2041-2060. Future values for BIO1, BIO2, BIO8, BIO12, BIO15 and BIO19 were extracted at the contemporary sampling coordinates with terra and transformed using the fitted Gradient Forest model.

For each genome, genomic offset was calculated as the Euclidean distance between its contemporary and future positions in the full Gradient Forest transformed space:

\[
\mathrm{Offset}_i =
\sqrt{\sum_{j=1}^{p}
\left(T_{ij}^{\mathrm{future}} - T_{ij}^{\mathrm{current}}\right)^2},
\]

where \(T_{ij}\) is the transformed value for environmental dimension \(j\). Genomic offset was therefore interpreted as relative projected displacement in contemporary genotype-environment relationships, not as a direct estimate of fitness loss or maladaptation. PCA of the contemporary and future transformed values was used only for visualization.

Spatial structure in genomic offset was evaluated after matching offsets to sampling coordinates. Spearman rank correlations were calculated between genomic offset and latitude and between genomic offset and absolute latitude. Geographic maps and continent-level summaries were used to describe broad spatial heterogeneity.

## Statistical analysis

Unless otherwise stated, statistical analyses were conducted in R. Non-parametric comparisons among lineages used Kruskal-Wallis tests, followed by Dunn tests where pairwise post hoc comparisons were required. Multiple comparisons were adjusted using the procedure specified for each analysis. Matrix-based and multivariate environmental analyses used permutation-based inference with 9,999 permutations unless otherwise stated. Spearman rank correlations were used for monotonic associations that did not require linearity or normality. All tests were two-sided.
