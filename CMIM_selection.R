#install.packages("infotheo")
library(mlbench)
library(infotheo)
data("BreastCancer") # Przykladowe dane
obserwacje_bez_NA <- complete.cases(BreastCancer) 
X <- BreastCancer[obserwacje_bez_NA , -c(1, 11)]
y <- BreastCancer[obserwacje_bez_NA , 11]


CMIM_selection <- function(X,y,kmax){
  stopifnot(is.numeric(kmax),kmax <= ncol(X),is.matrix(X) | is.data.frame(X))
  if (kmax == 1){ ### Specjalny przypadek dla kmax=1
    wybrane_zmienne <- list()
    zmienne <- 1:ncol(X)
    entropia <- unlist(lapply(X,function(x) entropy(x) - condentropy(X = x,
                                                                     Y = y,
                                                                     method = "emp")))
    wybrane_zmienne[["S"]] <- which.max(entropia)
    wybrane_zmienne[["score"]] <- max(entropia)
    names(wybrane_zmienne[["S"]]) <- NULL
    wybrane_zmienne
  }
  else{
    ####KROk 1
    wybrane_zmienne <- list()
    zmienne <- 1:ncol(X)
    entropia <- unlist(lapply(X,function(x) entropy(x) - condentropy(X = x,
                                                                     Y = y,
                                                                     method = "emp")))
    wybrane_zmienne[["S"]] <- which.max(entropia)
    wybrane_zmienne[["score"]] <- max(entropia)
    zmienne <- zmienne[-wybrane_zmienne$S]
    ####KROK DRUGi i POZOSTALE
    mat <- matrix(nrow = length(zmienne)
                  ,ncol = 1)
    colnames(mat) <- wybrane_zmienne$S
    rownames(mat) <- zmienne
    for (zm in 1:length(zmienne)) {
      mat[zm,1] <- condinformation(X = y,Y = X[,zmienne[zm]],
                                   S = X[,as.numeric(
                                     wybrane_zmienne$S[length(wybrane_zmienne$S)])]
                                   ,method = "emp")
    }
    entropia <- apply(mat,MARGIN = 1,min)
    wybrane_zmienne[["S"]] <- c(wybrane_zmienne[["S"]],
                                as.numeric(names(which.max(entropia))))
    wybrane_zmienne[["score"]] <- c(wybrane_zmienne[["score"]],
                                    max(entropia))
    zmienne <- c(1:ncol(X))[-as.numeric(wybrane_zmienne$S)]
    mat <- mat[rownames(mat) %in% zmienne ]
    if (kmax != 2){ # jezeli kmax=2 to nie musimy tego wykonwywac (nie powinnismy)
      for (k in 1:(kmax-2)){
        column <- matrix(nrow = length(zmienne)
                         ,ncol = 1)
        for (zm in 1:length(zmienne)) {
          column[zm,1] <- condinformation(X = y,Y = X[,zmienne[zm]],
                                          S = X[,as.numeric(
                                            wybrane_zmienne$S[length(wybrane_zmienne$S)])]
                                          ,method = "emp")
          
        }
        mat <- cbind(as.matrix(mat),column)
        rownames(mat) <- zmienne
        colnames(mat) <- wybrane_zmienne$S
        entropia <- apply(mat,MARGIN = 1,min)
        wybrane_zmienne[["S"]] <- c(wybrane_zmienne[["S"]],
                                    as.numeric(names(which.max(entropia))))
        wybrane_zmienne[["score"]] <- c(wybrane_zmienne[["score"]],
                                        max(entropia))
        zmienne <- c(1:ncol(X))[-as.numeric(wybrane_zmienne$S)]
        mat <- mat[rownames(mat) %in% zmienne,]
        mat <- matrix(mat,nrow = length(zmienne))
      }
    }
    wybrane_zmienne$S <- unname(wybrane_zmienne$S)
    wybrane_zmienne
  }
}

CMIM_sel <- CMIM_selection(X,y,9)

round(CMIM_sel$score,4)

CMIM_sel$S
