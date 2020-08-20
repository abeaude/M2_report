# Source useful function for the pipeline
source("utils.R")
# Source the folder corresponding to your analysis
sourceFolder("./Trajectory")

# Load your seurat object
load("/home/aurel/Nextcloud/Documents/Cours/M2_GENIOMHE/Stage/Project/NF1-MPNST/NF1-MPNST_20190212/ANALYSIS/DEFAULT/NF1-MPNST_20190212_SCT_DEFAULT_pca100.rda")

# Set the path of the project
# path <- "/home/aurel/Nextcloud/Documents/Cours/M2_GENIOMHE/Stage/Project/NF1-MPNST/NF1-MPNST_20190212/"
path <- '/home/aurel/test/'
# A subfolder TRAJECTORY will be created in this folder to store all your results

# Prepare the data
## Subset if neeed your object to only some cells, clusters
conditions <- list(idents = c(1,4,7,8,12))
root_cell <- select_cells(sobj,'seurat_clusters', "umap")

# Prepare data according to your subset
data_preparation_trajectories(sobj = sobj, method = c("slingshot","scorpius", "monocle3", "stream", "paga"), subset = conditions, cell_groups = "seurat_clusters", path = path, dimred = c("umap","pca"), start_cell = root_cell, seed = 1337L, filename = 'data_for_trajectories')

# Create another subset
cells <- Seurat::WhichCells(sobj,idents = c(1,4,7,8,12))

sobj_sub <- subset(sobj, idents = c(1,4,7,8,12))
cells_to_remove <- select_cells(sobj_sub,'seurat_clusters', "umap")
cells <- cells[-which(cells %in% cells_to_remove)]
conditions <- list('cells' = cells)

# Prepare data according to your subset
data_preparation_trajectories(sobj = sobj, method = c("slingshot","scorpius", "monocle3","paga"), subset = conditions, cell_groups = "seurat_clusters", path = path, dimred =  c("umap","pca"), start_cell = root_cell, seed = 1337L, filename = 'data_for_trajectories')

# Now run the trajectories inference
trajectories <- run_ti(path = path, filename = 'data_for_trajectories')
# You can now save each of your trajectory to the corresponding folder
save_trajectories(trajectories, path)
# Generate different plot for your trajectory
traj_plot <- plot_trajectories(trajectories)
save_trajectories(traj_plot,path)

# We can project all the trajecotries to a specific space like the initial UMAP
umap <- Seurat::Reductions(sobj,slot = "umap")@cell.embeddings
projected_traj <- project_traj(trajectories, space = umap, space_name = "UMAP")
save_trajectories(projected_traj,path)

# Let's have some metrics to comapre the trajecotry
metric_traj <- eval_trajectory(trajectories)
# plotting the results
plot_eval <- plot_trajectory_metrics(metric_traj)
save_eval_plot(plot_eval, path)

            # If you need to reload your result you can use this function. 
# You can also specify a subset and a methods
trajecotries <- load_all_traj(path)

# If you work on an integrated object you might want to visualize the different organisation of the sample
# on the trajectory

# make a list of extra gouping. Each dataframe must have 2 columns : "cell_id" and "group_id"
extra_grouping <- list(orig.ident = sobj@meta.data["orig.ident"] %>% 
       tibble::as_tibble() %>% 
       magrittr::set_colnames(c("cell_id", "group_id")))

# now do the plotting with this extra grouping
plot_trajectories(trajectories, extra_grouping = extra_grouping)