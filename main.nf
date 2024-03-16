#!/usr/bin/env nextflow
// Nextflow DSL2 pipeline for processing scRNA data using Seurat
//shebang line to define the environment

//enable DSL2 using NXF_VER=23.10.1
nextflow.enable.dsl=2

/**
 * Define the parameters for the script using the params object 
    * In Nextflow, `params` is a special object used to handle parameters that can be passed to the script at runtime. 
    * These parameters can be used to customize the behavior of the script without having to modify the script itself.
    * nextflow run myscript.nf --data_folder /path/to/data --output_folder /path/to/results (cn also run the script like this)
    * In this command, the --data_folder and --output_folder options are used to provide new values for the params.data_folder and params.output_folder parameters, 
    * which will override the default values set in the script.
 */

// Define parameters
params.data_folder = './data'
params.output_folder = './results'

/**
    * "Processing data from folder: ${params.data_folder}"
    * "Saving results to folder: ${params.output_folder}"
    * The log message is defined using a multi-line string, which is enclosed in triple quotes (""").
    * This allows the log message to span multiple lines, which makes it easier to read and understand.
*/

log.info """\
    S E U R A T - R  - P I P E L I N E
    ===================================
    data: ${params.data_folder}
    outdir: ${params.output_folder}
    """
    .stripIndent(true) // Strip leading whitespace from the log message

/**
 * Define a process named SeuratAnalysis to process the scRNA data using Seurat starts with process
    * In Nextflow, `process` is a special keyword used to define a process that will be executed as part of the workflow.
    * Each process is defined using a block of code that contains the process name, input and output declarations, and a script to be executed.
    * The input and output declarations are used to define the input and output channels for the process, which are used to connect the process to other processes in the workflow.
    * The script block contains the code that will be executed when the process is run.
    * The process name is used to identify the process in the workflow, and can be used to tag the process for better traceability.
 */
    process SeuratAnalysis { 
    tag "$h5_id"    // Tag the process with the sample ID (in this case h5_id) for better traceability otherwise sample_id can be used
    publishDir "${params.output_folder}", mode: 'copy'

    input:
    path h5_file // Input channel for the process, which accepts a list of files matching the pattern /path/to/data/*.h5 you can also name it as path input_file

    output:
    file("${h5_file}_seurat_results.rds")

    script:
    h5_id = h5_file.baseName
    """
    Rscript -e \"
                #!/usr/local/bin/ Rscript
                pacman::p_load(
                Seurat, ggplot2, patchwork, tidyverse, hdf5r, ggsci, 
                celldex, RColorBrewer, SingleCellExperiment, glmGamPoi, 
                reticulate, cowplot, viridis, pheatmap, scran, SingleR, 
                BiocParallel, DoubletFinder,presto); \\
                adj_matrix <- Read10X_h5('${h5_file}',use.names = T);\\
                seurat_obj <- CreateSeuratObject(counts = adj_matrix,project = 'TISCH2', min.cells = 3, min.features = 200);\\
                seurat_obj[['percent.mt']] <- PercentageFeatureSet(seurat_obj, pattern = '^MT');\\
                seurat_obj[['percent.rb']] <- PercentageFeatureSet(seurat_obj, pattern = '^RP[SL]');\\
                saveRDS(seurat_obj, file='${h5_file}_seurat_results.rds');\\
                adj.matrix <- NULL;\\
                str(seurat_obj);\\
                seurat_obj <- NormalizeData(seurat_obj, normalization.method = 'LogNormalize', scale.factor = 10000);\\
                seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = 'vst', nfeatures = 2000);\\
                seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj), vars.to.regress = NULL );\\
                seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj));\\
                seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10);\\
                seurat_obj <- FindClusters(seurat_obj, resolution = 0.5);\\
                seurat_obj <- RunUMAP(seurat_obj, dims = 1:10);\\
                seurat_obj <- FindAllMarkers(seurat_obj, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)"
    """
    }
workflow {
    Channel.fromPath("${params.data_folder}/*.h5") | SeuratAnalysis
}