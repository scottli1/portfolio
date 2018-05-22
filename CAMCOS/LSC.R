# install.packages("pdist")
# install.packages("clue")
# install.packages("R.matlab")
# install.packages("RSpectra")
# install.packages("Matrix")

library(Matrix)
library(RSpectra)
library(pdist)
library(clue)
library(R.matlab)

LSC <- function(data,                       # input feature matrix 
                k,                          # number of clusters
                p = NULL,                   # number of landmarks
                r = NULL,                   # number of nearest landmarks 
                t = 0,                      # time-step parameter 
                method = "random",          # landmark selection method
                similarity = "cosine",      # similarity measure
                clustering = "traditional", # clustering method
                seed = 0,                   # initial seed
                iter.max = 100,             # max iterations for k-means
                nstart = 10,                # sets of initializations for k-means
                reclassify = TRUE,          # reclassify outliers or give zero label
                alpha1 = 0,                 # data outlier removal
                alpha2 = 0,                 # landmark  outlier removal 
                beta = 0.6,                 # variance tuning parameter
                knn = 1) {                  # nearest neighbors used for knn
  
  # Set seed
  set.seed(seed)
  
  # Make sure sparse storage is used if applicable
  data <- Matrix(data)
  
  # Convert alpha1 and alpha2 from percentage to proportion
  alpha1 <- alpha1/100
  alpha2 <- alpha2/100
  
  
  # Some rules of thumb to pick p and r
  if (is.null(p)) {
    p <- floor(sqrt(nrow(data) * k))
  }
  
  if (is.null(r)) {
    r <- floor(p / 10)
  }
  
  if (method == "random") {
    # Sampling p landmarks randomly
    pindex <- sample(nrow(data), p)
    pmatrix <- data[pindex, ]
    # Cosine similarity by default. Gaussian otherwise.
    if (similarity == "cosine") {
      # Compute affinity matrix with cosine similarity
      A <- (data / sqrt(rowSums(data ^ 2))) %*% t(pmatrix / sqrt(rowSums(pmatrix ^ 2)))
    } else {
      # Compute distance matrix for p landmarks for estimating sigma
      p_dist <- as.matrix(dist(pmatrix, method = "euclidean"))
      keepr.small <- function(x, m) {
        return(sort(x)[1:m + 1])
      }
      m <- 2
      p_dist <- t(apply(p_dist, 1, keepr.small, m = m)) / 2
      # Estimate sigma using mean distance between two closest landmarks
      g_sigma <- beta * mean(p_dist)
      # Calculate euclidean distance between points and landmarks
      A <- as.matrix(pdist(data, pmatrix))
      A <- 1 / exp((A ^ 2) / (2 * (g_sigma ^ 2)))
    }
    
  } else {
    # Sampling p landmarks with kmeans
    if (similarity == "cosine") {
      # Suppress warnings because failure to converge is irrelevant
      suppressWarnings(pmatrix <- kmeans(data / sqrt(rowSums(data ^ 2)), centers = p)$centers)
    } else {
      suppressWarnings(pmatrix <- kmeans(data, centers = p)$centers)
    }
    if (similarity == "cosine") {
      # Compute affinity matrix with cosine similarity
      A <- (data / sqrt(rowSums(data ^ 2))) %*% t(pmatrix / sqrt(rowSums(pmatrix ^ 2)))
    } else {
      # Compute affinity matrix with gaussian similarity
      pfit <- kmeans(data, centers = p)
      pmatrix <- pfit$centers
      # Use total within sum of squares to estimate sigma
      g_sigma <- beta * sqrt(pfit$tot.withinss / (nrow(data) - p))
      A <- as.matrix(pdist(data, pmatrix))
      A <- 1 / exp((A ^ 2) / (2 * (g_sigma ^ 2)))
    }
  }
  
  # Landmark outlier removal. Outliers defined as having small column sums.
  if ((alpha2 * p) >= 1) {
    colsums <- colSums(A)
    landmark_quantile <- quantile(colsums, alpha2)
    A <- A[, colsums >= landmark_quantile]
    if (method == "random") {
      pindex <- pindex[colsums >= landmark_quantile]
    }
  }
  
  # Check for rowSums == 0. Give 0 label if there are any.
  zerooutliers <- FALSE
  if (any(rowSums(A) == 0)) {
    zerooutliers <- TRUE
    zoutliers <- rowSums(A) == 0
    znum <- sum(zoutliers)
    warning(znum, " row(s) of zeros in affinity matrix. Will give zero label.")
    A <- A[!zoutliers,]
  }
  
  # Data outlier removal. Outliers defined as having small row sums.
  outliers <- FALSE
  if ((alpha1 * nrow(A)) >= 1) {
    outliers <- TRUE
    rowsums <- rowSums(A)
    obs_quantile <- quantile(rowsums, alpha1)
    outlier_index <- rowsums < obs_quantile
    
    if (method == "random" & clustering != "traditional") {
      # If random sampling and landmark clustering, outliers can't be landmarks
      outlier_index[which(outlier_index)[which(outlier_index) %in% pindex]] <- FALSE
      # Shift pindex for splitting into test and training sets later
      porder <- order(order(pindex))
      pindex <- which((1:length(outlier_index) %in% pindex)[!outlier_index])
      pindex <- pindex[porder]
    }
    
    outlier_similarity <- A[outlier_index, ]
    A <- A[!outlier_index,]
    
    # Can give outliers zero label or reclassify with knn
    if (reclassify == TRUE) {
      knnindex <- function(x, length, knn) {
        if (knn == 1) {
          return(which.max(x))
        } else {
          return(order(x, decreasing = TRUE)[1:knn])
        }
      }
      # Construct nearest neighbor matrix for classification of outliers
      oknnmatrix <- Matrix(t(apply(outlier_similarity, 1, knnindex, length = ncol(A), knn = knn)))
    }
  }
  
  if (clustering == "traditional") {
    keepr <- function(x, length, r) {
      # Keep r largest entries in each row
      x[order(x)[1:(length - r)]] <- 0L
      return(x)
    }
    # Keep r nearest landmark distances
    A <- Matrix(t(apply(A, 1, keepr, length = ncol(A), r = r)))
  } else {
    keeprknn <- function(x, length, r, knn) {
      # Keep r largest entries in each row
      sorted <- order(x)
      x[sorted[1:(length - r)]] <- 0L
      return(c(sorted[length - ((knn - 1):0)], x))
    }
    A <- Matrix(t(apply(A, 1, keeprknn, length = ncol(A), r = r, knn = knn)))
    knnmatrix <- A[, 1:knn]
    A <- A[, -(1:knn)]
  }
  
  # Remove columns with colSums == 0
  if (any(colSums(A) == 0)) {
    colzeroindex <- colSums(A) != 0
    A <- A[, colzeroindex]
    if (method == "random") {
      pindex <- pindex[colzeroindex]
    }
  }
  
  if (clustering == "traditional") {
    if (t > 0) {
      A1  <- A / rowSums(A)  # Calculate A1
      A2 <- t(t(A1) / sqrt(rowSums(t(A1))))  # Calculate A2
      svdresult <- svds(A2, k = k, nu = k, nv = 0)
      sigma <- svdresult$d
      U <- svdresult$u
      U <- t(t(U) * sigma^t)
      U <- U / sqrt(rowSums(U ^ 2))
      labelsout <- kmeans(U, centers = k, iter.max = iter.max, nstart = nstart)$cluster
    } else {
      A1  <- A / rowSums(A)  # Calculate A1
      A2 <- t(t(A1) / sqrt(rowSums(t(A1))))  # Calculate A2
      # Take top k left singular vectors of A2
      U <- svds(A2, k = k, nu = k, nv = 0)$u
      # L2 normalize rows of U
      U <- U / sqrt(rowSums(U ^ 2))
      # Kmeans on U
      labelsout <- kmeans(U, centers = k, iter.max = iter.max, nstart = nstart)$cluster
    }
  } else {
    A1  <- t(t(A) / (rowSums(t(A))))  # Calculate A1
    A2  <- A1 / sqrt(rowSums(A1))  # Calculate A2
    # Take top k right singular vectors of A2
    V <- svds(A2, k = k, nu = 0, nv = k)$v
    # L2 normalize rows of V
    V <- V / sqrt(rowSums(V ^ 2))
    # Kmeans on V
    landmarkfit <- kmeans(V, centers = k, iter.max = iter.max, nstart = nstart)$cluster
    
    if (method == "random") {
      #Split test and train sets
      fullindex <- 1:nrow(A2)
      trainindex <- pindex
      testindex <- fullindex[-trainindex]
      if (knn == 1) {
        knnmatrix <- knnmatrix[-trainindex]
      } else {
        knnmatrix <- knnmatrix[-trainindex, ]
      }
      votematrix <- matrix(landmarkfit[as.vector(knnmatrix)], nrow = length(testindex), ncol = knn)
      #Do knn on data
      if (knn == 1) {
        knnresults <- votematrix
      } else {
        getmode <- function(x) {
          tab <- table(x)
          max <- names(which(tab == max(tab)))
          if (length(max) > 1) {
            max <- sample(max, 1)
          }
          return(as.integer(max))
        }
        knnresults <- apply(votematrix, 1, getmode)
      }
      # Recombine labels
      predictedlabels <- rep(NA, length.out = nrow(A2))
      predictedlabels[testindex] <- knnresults
      predictedlabels[trainindex] <- landmarkfit
      labelsout <- predictedlabels
    } else {
      votematrix <- matrix(landmarkfit[as.vector(knnmatrix)], nrow = nrow(A), ncol = knn)
      if (knn == 1) {
        labelsout <- votematrix
      } else {
        getmode <- function(x) {
          tab <- table(x)
          max <- names(which(tab == max(tab)))
          if (length(max) > 1) {
            max <- sample(max, 1)
          }
          return(as.integer(max))
        }
        labelsout <- apply(votematrix, 1, getmode)
      }
    }
  }
  
  # Reclassify outliers or give zero label
  if (outliers == TRUE) {
    if (reclassify == TRUE) {
      ovotematrix <- matrix(labelsout[as.vector(oknnmatrix)], nrow = nrow(oknnmatrix), ncol = knn)
      if (knn == 1) {
        predicted.o.label <- ovotematrix
      } else {
        getmode <- function(x) {
          tab <- table(x)
          max <- names(which(tab == max(tab)))
          if (length(max) > 1) {
            max <- sample(max, 1)
          }
          return(as.integer(max))
        }
        predicted.o.label <- apply(ovotematrix, 1, getmode)
      }
      finallabels <- rep(NA, length.out = length(rowsums))
      finallabels[outlier_index] <- predicted.o.label
      finallabels[!outlier_index] <- labelsout
      labelsout <- finallabels
    } else {
      finallabels <- rep(NA, length.out = length(rowsums))
      finallabels[outlier_index] <- 0
      finallabels[!outlier_index] <- labelsout
      labelsout <- finallabels
    }
  }
  
  if (zerooutliers == TRUE) {
    finallabels <- rep(NA, length.out = nrow(data))
    finallabels[zoutliers] <- 0
    finallabels[!zoutliers] <- labelsout
    labelsout <- finallabels
  }
  
  # Check for missing clusters. Can only occur with kmeans landmark selection and 
  # landmark clustering. 
  if (length(unique(labelsout[labelsout != 0])) != k) {
    warning("The final number of clusters is not equal to k")
  }
  
  return(as.vector(labelsout))
}

accuracy <- function(truelabels, clusters) {
  # Hungarian algorithm
  
  # Remove zeros from labels
  if (any(clusters == 0)) {
    print("Zero labels will be removed from accuracy computation")
    zeroindex <- clusters == 0
    truelabels <- truelabels[!zeroindex]
    clusters <- clusters[!zeroindex]
  }
  
  # Labels from cluster A will be matched on the labels from cluster B
  minWeightBipartiteMatching <- function(clusteringA, clusteringB) {
    require(clue)
    idsA <- unique(clusteringA)  # distinct cluster ids in a
    idsB <- unique(clusteringB)  # distinct cluster ids in b
    nA <- length(clusteringA)  # number of instances in a
    nB <- length(clusteringB)  # number of instances in b
    if (length(idsA) != length(idsB) || nA != nB) {
      stop("The number of clusters or lengths do not match")
    }
    nC <- length(idsA)
    tupel <- c(1:nA)
    # Computing the assignment matrix
    assignmentMatrix <- matrix(rep(-1, nC * nC), nrow = nC)
    for (i in 1:nC) {
      tupelClusterI <- tupel[clusteringA == i]
      solRowI <- sapply(1:nC, function(i, clusterIDsB, tupelA_I) {
        nA_I <- length(tupelA_I)  # number of elements in cluster I
        tupelB_I <- tupel[clusterIDsB == i]
        nB_I <- length(tupelB_I)
        nTupelIntersect <- length(intersect(tupelA_I, tupelB_I))
        return((nA_I - nTupelIntersect) + (nB_I - nTupelIntersect))
      }, clusteringB, tupelClusterI)
      assignmentMatrix[i,] <- solRowI
    }
    # Optimization
    result <- solve_LSAP(assignmentMatrix, maximum = FALSE)
    attr(result, "assignmentMatrix") <- assignmentMatrix
    return(result)
  }
  test <- minWeightBipartiteMatching(clusters, truelabels)
  predicted = NULL
  predicted <- rep(NA, length(clusters))
  for (i in 1:length(test)) {
    predicted[which(clusters == i)] <- test[i]
  }
  table <- table(predicted, truelabels)
  accuracy <- (sum(diag(table)) / length(truelabels))
  classaccuracy <- vector()
  colsums <- colSums(table)
  for (i in 1:length(test)) {
    classaccuracy[i] <- table[i, i] / colsums[i]
  }
  return(
    list(
      "accuracy" = accuracy,
      "classaccuracy" = classaccuracy,
      "table" = table,
      "mapping" = test,
      "mappedlabels" = predicted
    )
  )
}