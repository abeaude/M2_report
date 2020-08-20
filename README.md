# Installation

## Install system dependencies
``` bash
sudo apt-get install libhdf5-dev libgsl0-dev libzzip-dev zlibc libc6 libcairo2-dev git curl libcurl4-openssl-dev libxml2-dev libudunits2-dev libgdal-dev libgeos-dev libproj-dev
pip install untangle
conda create -n scRNAseq_qc -c bioconda fastqc fastq-screen multiqc
conda activate scRNAseq_qc
```

Run `fastq_screen --get_genomes` in the folder where you want to store the different genomes used by `fastq_screen`. A folder called "FastQ_Screen_Genomes" is created. In this folder you will find a ".conf" file, this file is used by `fastq_screen` to run correctly. Open the file in a text editor, remove the "#" in front of the BOWTIE2 line. Now you need to set the correct path for the bowtie2 binary. It should be "~/anaconda3/envs/scRNAseq_qc/bin/bowtie2", replace any another path by this one in the conf file.

Now you can copy this file in the fastq_screen folder : "~/anaconda3/envs/scRNAseq_qc/share/fastq-screen-0.14.0-0"

## CellRanger installation 
Download the last version of bcl2fastq on this page (rpm) : https://support.illumina.com/downloads/bcl2fastq-conversion-software-v2-20.html

### bcl2fastq
``` bash
sudo apt-get install alien
sudo alien bcl2fastq2-v2.20.0.422-Linux-x86_64.rpm
sudo dpkg -i bcl2fastq2_0v2.20.0.422-2_amd64.deb
sudo apt-get install -f
```

### Cellranger
``` bash
cd /opt
sudo wget -O cellranger-3.1.0.tar.gz link
sudo tar -xzvf cellranger-3.1.0.tar.gz 
```

Add the following line to your .bashrc `export PATH="$PATH:/opt/cellranger-3.1.0"`. And source your bashrc to take changes into account. 

### List of correct index as of 06/2020
The file `cellranger_correct_index` contains all the correct index used by 10x genomics as of June 2020. They are use by `cellranger.py` to validate the index in the database of each project. 
To update the file just add new indexes at the end of the file with one index per line. 

## Setting up the project
Clone the repository : 
``` bash
git clone https://github.com/abeaude/scRNAseq.git
```

You need a working installation of anaconda. 
Init the project (choose 3)
``` R
renv::init()
```

Install tensorflow to use cellassign
``` R
install.packages('tensorflow')
tensorflow::install_tensorflow(method = "conda",extra_packages='tensorflow-probability', version = "2.2.0",envname = 'r-reticulate')
```

``` R
devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")
```

Install the required packages by using the `renv.lock` file:
``` R
renv::restore()
```

## Create the different conda environement
``` bash
conda install nb_conda_kernels
conda create --quiet --yes -n stream -c conda-forge -c bioconda python=3.6 ipython ipykernel r-irkernel stream=1.0 

conda create -n scanpy -c conda-forge -c bioconda python-igraph leidenalg seaborn scikit-learn statsmodels numba pytables scanpy=1.5.1 
conda activate scanpy
pip install fa2 loompy
conda deactivate

conda create -n cellphonedb python
conda activate cellphonedb
pip install cellphonedb  
conda deactivate
```

# Usage
___
## Cell ranger
The python script `cellranger.py` will help you running the cellranger pipeline on your data. Data are organised on per-project basis. A folder will store all the on going project (`PROJECT_DIR`), each project correspond to a folder containing all the experiment related to the project. At the root of the project folder, you will find an xslx file containing the list of the different experiment and some meta-information like the index used for sequencing. This file is important to run the programm as it will sepcify the correct index to use. 

```
usage: cellranger.py [-h] --run {fastq,count} --project PROJECT --experiment
                     EXPERIMENT --project-dir PROJECT_DIR [--qc-env QC_ENV]
                     [--transcriptome TRANSCRIPTOME]

optional arguments:
  -h, --help            show this help message and exit
  --run {fastq,count}   Specify if you want to run the fastq generation or the 
                        count step of the cell ranger pipeline
  --project PROJECT     Specify on which project you want to run the pipeline on
  --experiment EXPERIMENT
                        Specify on which experiment of the project you want to 
                        run the pipeline on
  --project-dir PROJECT_DIR
                        Folder where all the project are stored
  --qc-env QC_ENV       Conda environment to use
  --transcriptome TRANSCRIPTOME
                        Specify the transcriptome to use in the count step

```

After running the script the folder structure will be the following: 
```
PROJECT_DIR
├── PROJECT1
│   ├── EXPERIMENT1
│   │   │── BCL
│   │   │── fastqs
│   │   │── QC
│   │   │── counts
│   │   │── cellranger_run
│   ├── EXPERIMENT2
│   └── EXPERIMENT3
├── PROJECT2
│   ├── EXPERIMENT1
│   └── EXPERIMENT2
└── PROJECT3
```

## Trajectory

In order to run trajectories inference you a seurat object with different dimension reductions. The pipeline is divided in different step : 
1. Preparation of the data for the different method (you can call it as many times as you want)
2. Trajectory inference
3. Plotting and projection
4. Trajectory evaluation

To see a complete example on how to run trajectories on you can see : [Trajectory inference](Example_Scripts/traj_grafted_analysis.R)

### Usage : `data_preparation_trajectories`
> **[parameters]**
>* sobj: a seurat object
>* method: methods of trajectory inference to use (curretnly available : "slingshot","scorpius", "monocle3", "stream", "paga")
>* subset: A list to specify a subset on which run trajectory inference
>* cell_groups: How the cells should be grouped in the returned dynverse object (must be an available field in the metadata of the seurat object)
>* path: Path where to store the different object created
>* dimred: Which previously computed dimension reduction to use (Should be available in Seurat::Reductions(sobj)), if NULL, will use all the available dimensions reductions
>* start_cell: Specify the cell considered as a starting point of the trajectory
>* seed: Specify a seed to use, otherwise will use seed stored in sobj@misc$params$seed or finally 1337
>* filename: How to prefix all the output. Default : 'data_for_trajectories'

For each specified method, it will call the function `data_prep_method` to prepare the data in the requested format by the method and a yaml file containing different section. Each section represent a set of parameters to use in trajectory inference. The function will create some 'default' section but user can add afterwards as many section as needed. 

### How to specify a subset
It is not always necessary to run trajectory on all your cells, you can decide to restrict your analysis to a group of cells. This section will explain you how to use the the subset argument of the `data_preparation_trajectories` function. 

The subset argument must be a named list with names among 'cells', 'idents' or any metadata entry of your seurat object. If your list contains a 'cells' entry it will be used in priority ignoring any other entry. If 'cells' is not used but 'idents' is, it will be used in priority. Finally if neither 'cells' or 'idents' is used, it will use the different metadata entry and combine the result with a logical AND to select the cells of interest.  

When using 'cells' subsetting, your list must a contains a vector of valid cells in the 'cells' entry. When using 'idents' your list must a contains a vector of valid idents in the 'idents' entry. The same goes for any metadata entry you planned to use, otherwise the programm will throw an error. 

### How to choose the start cell
There are no magic method to select the starting point of the trajectory. You have to decide which cell to choose, based on your knowledge, the biology or arbitray. This section will only show you how to retrieve interactively the name of a cell. 

By running the following code, you will be able to hover on the different cells and see their name.
``` R
plot <- Seurat::DimPlot(sobj) # create a dim plot
Seurat::HoverLocator(plot = plot, information = Seurat::FetchData(sobj, vars = c("ident"))) # use the dimplot and make it interactive
```

A functian was added to allow a better workflow of cell selection : `select_cells`. You just need to call the function `select_cells(sobj, "seurat_clusters")` to get an interactive plot allowing you to choose on or more cells and get their name in return. This function can also be used to select multiple cells to be removed when subsetting, see [Trajectory inference](Example_Scripts/traj_grafted_analysis.R). 

Some methods needs a start cell to infer the the trajectory, other don't. But all trajecotries will be rooted using the start cell. This root is then used to compute pseudotime along the trajectory. 

### Usage : `run_ti`
> **[parameters]**
>* path: Path of an experiment
>* filename: filename used in `data_preparation_trajectories`

This function will run the trajecotry inference on the different object created with `data_preparation_trajectories`. For each method it will read the yaml file containing all the required parameters for the trajecotry inference. If a set of parameters was already used to produce a trajectory it will be skipped. This behaviour allow the user to add incrementally parameters sections while not recmoputing all trajecotries. 

### Usage : `plot_trajectories`
> **[parameters]**
>* trajectories : a list of trajectory as returned by `run_ti`
>* pb: a progress bar created with progress::progress_bar$new()
>* extra_grouping: A named list of dataframe. Each dataframe gives information about the grouping and must have 2 columns ("cell_id" and "group_id")

The function will return a list of plot with the same structure as the input. 

### Usage : `save_trajectories`
> **[parameters]**
>* trajectories: a list of trajectory as returned by `run_ti` or a list of plot as returned by `plot_trajectories`
>* path: Path of an experiment

This function is used to save both output of `run_ti` and `plot_trajectories`. See [Folder structure](#folder-structure) to see how trajectories and plot are saved. 

### Usage : `project_traj`
> **[parameters]**
>* trajectories: a list of trajectory as returned by `run_ti`
>* space: The dimensionality reduction of the cells. A matrix with the positions of cells (rows) in the dimensions (columns)
>* space_name: The name of the projection space (use in naming output plot)

This function will return a list of plot with the same structure as the input. To save the plot to the disk use `save_trajectories`

### Usage : `eval_trajectory`
> **[parameters]**
>* trajectories: a list of trajectory as returned by `run_ti`
>* metrics_id: Which metrics to evaluate. Check dyneval::metrics for a list of possible metrics. 

One use case for these metrics is to calculate the accuracy of a certain prediction compared to a reference trajectory. It will perform all the pairwise comaprison inside one subset. 

### Usage : `plot_trajectory_metrics`
> **[parameters]**
>* metric_res: a tibble or a dataframe as returned by `eval_trajectory`
>* metrics_id: Same metrics_id as used in `eval_trajectory`

Plot the results of the different trajecotry metrics. 

### Usage : `save_eval_plot`
> **[parameters]**
>* eval_plot: a list of plot as returned by `plot_trajectory_metrics`
>* path: Path of an experiment

This function is used to save the plot produced by `plot_trajectory_metrics`

### Usage : `load_all_traj`
> **[parameters]**
>* path: Path of an experiment
>* subset: Specify a particular subset
>* method: Specify a particular method of interest

If one need to reload a set of trajectory to do downstream analysis it can use this function. subset and method argument can be used independently or together.

### Folder structure
![alt text](images/trajectory_folder.png "Trajectory Folder Structure")

## Differentially expression
### Usage : `run_multiple_DE`
> **[parameters]**
>* parameters: If null will launch the shiny app
>* sobj: A seurat object
>* path: Path of an experiment
>* min.pct: only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations. Meant to speed up the function by not testing genes that are very infrequently expressed. Default is 0.1
>* .parallel: logical to indicate wheter or not run differential expression in parallel
>* verbose: Verbosity (TRUE/FALSE)

This function allow to run multiple differential expression at once. It will launch a shiny app allowing you to choose the comparisons to make. More information are avialble in the shiny app. 

### Usage : `run_DE`
> **[parameters]**
>* sobj: A seurat object
>* clusters: a character specifying a column metadata to use as cluster information or a dataframe with two columns : cells and cluster. The cluster column must be a factor. If NULL, will use the default idents of the seurat object (Idents(sobj))
>* DE_type: Can be one of "1vsAll", "1vs1", "SvsS", "conditions"
>* method: Can be one of "Limma_Voom", "Limma_Trend", "EdgeR_LRT", "EdgeR_QL"
>* testToUse: Can be one of "wilcox", "bimod", "t", "poisson", "negbinom", "LR", "MAST", "roc", "DESeq2". Those are test used by Seurat::FindMarkers
>* first_group: First group of clusters
>* second_group: Second group of clusters
>* path: Path of an experiment
>* min.pct: only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations. Meant to speed up the function by not testing genes that are very infrequently expressed. Default is 0.1
>* batch:
>* conditions_cluster : On which cluster run the differential expression with DE_type = "conditions". This will subset your seurat object to the according idents before any computation. 
>* .parallel: logical to indicate wheter or not run differential expression in parallel (Work only outside of rstudio)
>* verbose: Verbosity (TRUE/FALSE)

### Usage : `stat_DE`
> **[parameters]**
>* path: Path of an experiment
>* pval: Threshold of p-values to select significative genes. Default is 0.05
>* thresholdLogFC: Threshold of LogFC to select DE genes. Default is 0.58
>* thresholdpct: Threshold for percentage of untranslated genes for enrichment. It will throw a warning if above this threshold. Default is 5(%)
>* species: A supported species. To see supported species run `get_db()`.Default is "mouse"
>* verbose: Verbosity (TRUE/FALSE)

## Functional enrichment
### Usage : `run_functionnal_enrichment`
> **[parameters]**
>* path: Path of an experiment
>* pval: Threshold of p-values to select significative genes. Default is 0.05
>* thresholdLogFC: Threshold of LogFC to select DE genes. Default is 0.58
>* thresholdpct: Threshold for percentage of untranslated genes for enrichment. It will throw a warning if above this threshold. Default is 5(%)
>* species: A supported species. To see supported species run `get_db()`.Default is "mouse"
>* category_msigdb: One or more of "H", "C1", "C2", "C3", "C4", "C5", "C6", "C7"
>* method: One or more of "WikiPathway", "msigdb", "CellMarkers", "GOterms", "ReactomePA", "DiseaseOntology"
>* gmt_date: Release date to ue for gmt file (used for WikiPathways)
>* gmt_folder: Folder containing the gmt file
>* enrichment: Run enrichment or not. Default is TRUE
>* GSEA: Run GSEA or not. Default is TRUE
>* comparison: Run comparisons or not. Default is TRUE
>* verbose: Verbosity (TRUE/FALSE). Default is TRUE
>* force: Recompute everything even contrast already done. Default is FALSE

### Usage : `download_gmt`
> **[parameters]**
>* species: A supported species. To see supported species run `get_db()`
>* wikipathways_folder: Folder containing the different gmt files used by wiki pathways.
>* date: Specify the date of the wikipathway release to use, if NULL will download the latest version available. Default is NULL

### Usage : `get_db`
> **[parameters]**
>* species : One of "anopheles","arabidopsis","bovine","canine","celegans","chicken","chimp","coelicolor","ecolik12","ecsakai","fly","gondii","human","malaria","mouse","pig","rat","rhesus","xenopus","yeast","zebrafish"

If you run `get_db` without arguments you will obtain the list of supported species. This function is used to do matching between common species name and all the different database or name used whithin different packages. 

### Folder structure
![alt text](images/enrichment_folder.png "Enrichment Folder Structure")


## Interaction
### Usage : `run_CCInx`
> **[parameters]**
>* path: Path were to store the results
>* test: Test of differentially expression to use. If NULL use all available tests
>* clsuters: Cluster to use to infer interaction
>* species: Species to use, human or mouse

To run correctly, you need to run run previously the differential expression: `run_DE` between your different cluster of interest. 

### Usage : `run_cellphonedb`
> **[parameters]**
>* sobj: A seurat object
>* gene_type: Type of gene identifiers in the counts data [ensembl | gene_name | hgnc_symbol] 
>* species: Species to use, human or mouse
>* threshold: % of cells expressing the specific ligand/receptor
>* iterations: Number of iterations for the statistical analysis
>* pval: P-value threshold
>* path: Path were to store the results
>* filename: output filename prefix
>* verbose: Verbosity (TRUE/FALSE)

species = "human", threshold = 0.1, iterations = 1000, pval = 0.05
sobj, "gene_name", "mouse", path = path, filename = "test_phoneDB"

### Usage : `run_nichnetr`
> **[parameters]**
>* sobj: a seurat object
>* sender: Determine the potential sender cells. Name of cluster identity/identities of cells that presumably affect expression in the receiver cell type.
>* receiver: 	
Name of cluster identity/identities of cells that are presumably affected by intercellular communication with other cells
>* condition_oi: Condition of interest in which receiver cells were presumably affected by other cells. 
>* condition_reference: The second coPath were to store the resultsndition (e.g. reference or steady-state condition).
>* path: Path were to store the results
>* network_path: Path where are stored the network/datasets used by nichenetr
>* network_name: Name of the subfolder in network_path. Default: "mouse"
>* test: Test of differentially expression to use. If NULL use all available tests
>* pval: P-value threshold. Default : 0.05
>* logFC: log Fold-change threshold. Default: 0.58
>* expression.pct: To determine ligands and receptors expressed by sender and receiver cells, we consider genes expressed if they are expressed in at least a specific fraction of cells of a cluster. This number indicates this fraction. Must be the same as the one used in the differential expression part. Default: 0.10
>* top_n_ligands: Indicate how many ligands should be extracted as top-ligands after ligand activity analysis. Only for these ligands, target genes and receptors will be returned. Default: 20.
>* top_n_targets: To predict active, affected targets of the prioritized ligands, consider only DE genes if they also belong to the a priori top n ("top_n_targets") targets of a ligand. Default = 200.
>* cutoff_visualization: Because almost no ligand-target scores have a regulatory potential score of 0, we clarify the heatmap visualization by giving the links with the lowest scores a score of 0. The cutoff_visualization paramter indicates this fraction of links that are given a score of zero. Default = 0.33.

To run correctly, you need to run run previously the differential expression: `run_DE` between your different conditions on the receiver cells. 

## Integration

### Usage : `create_sobj_list`
> **[parameters]**
>* ... : Paths of the different seurat object to use. File extension must be .rda, .rds or .RData. 

### Usage : `liger_integration`
> **[parameters]**
>* sobj_list: A list of seurat object
>* k:
>* lambda:
>* split.by: 
>* project:
>* verbose:
>* seed: Specify a seed to use, otherwise will use seed stored in sobj_list[[1]]@misc$params$seed or finally 1337

sobj_list, k = 20, lambda = 5, split.by = "orig.ident", project = 'LIGER_integration', verbose = FALSE, seed = NULL

By default the liger intehration will use `k = 20` and `lambda = 5`. But liger provide function to estimate k and lambda for our data. In order to do so, you need to create a liger object with the function `create_liger_object(sobj_list)`, then you can run the liger pipeline to estimate those values.

```R
liger_object <- create_liger_object(sobj_list) %>%
  liger::normalize() %>%
  liger::selectGenes() %>%
  liger::scaleNotCenter()

liger::suggestK(liger_object, k.test = seq(5, 50, 5), rand.seed = 1337L, num.cores = 5)
liger::suggestLambda(liger_object, k = 20, lambda.test = seq(0.25,60,1), rand.seed = 1337L, num.cores = 5)
```

Those two functions will return a ggplot object, most appropriate k/lambda is likely around the "elbow" of the alignment plot (when alignment stops increasing). 

### Seurat integration

# Adding method for trajectory inference
___
The trajectory inference pipeline was constructed in way to ease the addition of new methods of inference. 