## CLustering backup
#!/usr/bin/env Rscript
library("tidyverse")
library("Seurat")
library("Matrix")
library("aricode")

rm(list=ls())

args = commandArgs(trailingOnly=TRUE)
case_pos <- args[1]
ctrl_pos <- args[2]
outdir <- args[3]

case_obj <- readRDS(file = case_pos)
ctrl_obj <- readRDS(file = ctrl_pos)

# Subset to common cells
common_cells <- intersect(colnames(case_obj),colnames(ctrl_obj))
ctrl_obj <- subset(ctrl_obj,cells=common_cells)
case_obj <- subset(case_obj,cells=common_cells)


# Normalize and find neighbors
case_obj <- NormalizeData(case_obj) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:30)
ctrl_obj <- NormalizeData(ctrl_obj) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:30)

resolution_list <- c(0.25,0.5,0.75,1.0)
cluster_names <- paste("cluster_res_",resolution_list,sep="")

for(i in 1:length(resolution_list)){
  case_obj <- FindClusters(case_obj,resolution = resolution_list[[i]],cluster.name = cluster_names[[i]])
  ctrl_obj <- FindClusters(ctrl_obj,resolution = resolution_list[[i]],cluster.name = cluster_names[[i]])
}


case_clusters <- lapply(cluster_names,function(x){
  FetchData(case_obj,vars = x)
}) %>% 
  bind_cols() %>% 
  rownames_to_column(var = "cell")

ctrl_clusters <- lapply(cluster_names,function(x){
  FetchData(ctrl_obj,vars = x)
}) %>% 
  bind_cols() %>% 
  rownames_to_column(var = "cell")

ARI(ctrl_clusters$cluster_res_0.2,case_clusters$cluster_res_0.2)

out_ARI <- lapply(cluster_names,function(x){
  ARI(case_clusters[[x]],ctrl_clusters[[x]])
})
names(out_ARI) <- resolution_list

json_obj <- list(
  module = "R3_clustering",
  timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ"),
  metrics = list(
    adjusted_rand_index = list(
      resolution_0.25 = out_ARI[["0.25"]],
      resolution_0.50 = out_ARI[["0.25"]],
      resolution_0.75 = out_ARI[["0.25"]],
      resolution_1.00 = out_ARI[["0.25"]]
    ),
    cluster_count = list(
      resolution_0.25_ctrl = length(unique(ctrl_clusters$cluster_res_0.25)),
      resolution_0.25_anon = length(unique(case_clusters$cluster_res_0.25)),
      resolution_0.5_ctrl = length(unique(ctrl_clusters$cluster_res_0.5)),
      resolution_0.5_anon = length(unique(case_clusters$cluster_res_0.5)),
      resolution_0.75_ctrl = length(unique(ctrl_clusters$cluster_res_0.75)),
      resolution_0.75_anon = length(unique(case_clusters$cluster_res_0.75)),
      resolution_1_ctrl = length(unique(ctrl_clusters$cluster_res_1)),
      resolution_1_anon = length(unique(case_clusters$cluster_res_1))
    )
  )
)

write_json(json_obj, outdir, pretty = TRUE, auto_unbox = TRUE)



