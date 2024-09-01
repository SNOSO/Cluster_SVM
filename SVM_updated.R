library(e1071)
library(mclust)
library(caret)

load_and_process_data <- function(data_type = "DLPFC", gene_list_path, data_path){
  # Load gene list
  gl <- readRDS(gene_list_path)
  
  # Load the data depending on the data type
  if (data_type == "DLPFC") {
    DLPFC <- readRDS(data_path)
    DLPFC_676 <- DLPFC[["151676"]]
    
    # Subset data
    DLPFC_676@assays$SCT@scale.data <- DLPFC_676@assays$SCT@scale.data[
      which(rownames(DLPFC_676@assays$SCT@scale.data) %in% gl), ]
    
    matrix <- DLPFC_676@assays$SCT@scale.data
    coordinates <- DLPFC_676@images$slice1@coordinates[, c(1:3)]
    labels <- DLPFC_676$Labs
    
  } else if (data_type == "CODEX") {
    codex.obj <- readRDS(data_path)
    
    matrix <- codex.obj@assays$RNA@layers$scale.data
    coordinates <- cbind(colnames(codex.obj), 
                         codex.obj$Xcorr, 
                         codex.obj$Ycorr)
    labels <- codex.obj$celltype
    
  } else if (data_type == "Xenium") {
    xenium.obj <- readRDS(data_path)
    
    matrix <- xenium.obj@assays$SCT@scale.data
    coordinates <- cbind(as.numeric(colnames(xenium.obj)),
                         xenium.obj@images$crop@boundaries$centroids@coords)
    labels <- xenium.obj$niches
  } else {
    stop("Unknown data type")
  }
  
  # Ensure coordinates are numeric and integers
  coordinates <- data.frame(coordinates)
  coordinates[,2] <- round(as.numeric(coordinates[,2]))
  coordinates[,3] <- round(as.numeric(coordinates[,3]))
  
  return(list(matrix = matrix, coordinates = coordinates, labels = labels))
}

RunSVM <- function(matrix = matrix, coordinates = coordinates, labels = labels, nsize = 3, kernel_type = "linear"){
  #randomness seed
  set.seed(22)
  
  #train/test split
  ind <- sample(2, ncol(matrix), replace = TRUE, prob = c(0.7, 0.3))
  
  #isolate the input for training
  mat <- matrix[,ind == 1]
  labs <- labels[ind == 1]
  coords <- coordinates[ind == 1,]
  rn <- levels(factor(labs))
  
  avg_mat <- mat
  for (i in 1:length(rn)){
    wc <- coords
    mat_1 <- mat[,which(labs == rn[i])]
    wc <- wc[which(labs == rn[i]),]
    
    for (j in 1:ncol(mat_1)){
      roi <- wc[j,2]
      coi <- wc[j,3]
      allrows <- wc[,2]
      allcols <- wc[,3]
      neighs <- which((allrows %in% c((roi-nsize):(roi+nsize))) & 
                        (allcols %in% c((coi-nsize):(coi+nsize))))
      
      if (length(neighs) < 2){
        next
      }
      
      newj <- rowMeans(mat_1[,neighs])
      avg_mat[,colnames(mat_1)[j]] <- newj
    }
  }
  message("Finished training neighborhood averaging")
  
  df <- data.frame(cbind(labs, t(avg_mat)))
  for (i in 2:ncol(df)){
    df[,i] <- as.numeric(df[,i])
  }
  df$labs <- factor(df$labs)
  train <- df
  
  #run SVM with specified kernel type
  svm_model <- svm(labs ~ ., data = train, kernel = kernel_type, scale = TRUE)
  
  message("Finished model building")
  
  #isolate the input for testing
  mat <- matrix[,ind == 2]
  labs <- labels[ind == 2]
  coords <- coordinates[ind == 2,]
  
  wc <- coords
  avg_mat <- mat
  mat_1 <- mat
  
  for (j in 1:ncol(mat_1)){
    roi <- wc[j,2]
    coi <- wc[j,3]
    allrows <- wc[,2]
    allcols <- wc[,3]
    neighs <- which((allrows %in% c((roi-nsize):(roi+nsize))) & 
                      (allcols %in% c((coi-nsize):(coi+nsize))))
    
    if (length(neighs) < 2){
      next
    }
    
    newj <- rowMeans(mat_1[,neighs])
    avg_mat[,colnames(mat_1)[j]] <- newj
  }
  
  message("Finished testing neighborhood averaging")
  
  df <- data.frame(cbind(labs, t(avg_mat)))
  for (i in 2:ncol(df)){
    df[,i] <- as.numeric(df[,i])
  }
  df$labs <- factor(df$labs)
  test <- df
  
  #predict test with SVM
  p1 <- predict(svm_model, test)
  
  #calculate ARI
  ari <- adjustedRandIndex(test$labs, p1)
  
  #output results
  return(list(ARI = ari, ConfusionMatrix = confusionMatrix(p1, test$labs)))
}

data_list <- load_and_process_data(data_type = "DLPFC", 
                                   gene_list_path = "DLPFCGeneList.RDS", 
                                   data_path = "SpatialSeurats.RDS")

SVM <- RunSVM(matrix = data_list$matrix, 
              coordinates = data_list$coordinates, 
              labels = data_list$labels, 
              nsize = 3, 
              kernel_type = "radial")

