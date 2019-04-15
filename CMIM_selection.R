install.packages("infotheo")
library(mlbench)
library(infotheo)
data("BreastCancer")
obserwacje_bez_NA <- complete.cases(BreastCancer) 
X <- BreastCancer[obserwacje_bez_NA , -c(1, 11)]
y <- BreastCancer[obserwacje_bez_NA , 11]



CMIM_selection <- function(X,y,kmax){
  stopifnot(is.numeric(kmax),kmax < ncol(X),is.matrix(X) | is.data.frame(X))
  if (kmax == 1){
    wybrane_zmienne <- list()
    zmienne <- 1:ncol(X)
    entropia <- unlist(lapply(X,function(x) entropy(x) - condentropy(X = x,
                                                                     Y = y,
                                                                     method = "emp")))
    wybrane_zmienne[["zmienne"]] <- which.max(entropia)
    wybrane_zmienne[["wynik"]] <- max(entropia)
    wybrane_zmienne
  }
  else{
  ####KROk 1
  wybrane_zmienne <- list()
  zmienne <- 1:ncol(X)
  entropia <- unlist(lapply(X,function(x) entropy(x) - condentropy(X = x,
                                                      Y = y,
                                                      method = "emp")))
  wybrane_zmienne[["zmienne"]] <- which.max(entropia)
  wybrane_zmienne[["wynik"]] <- max(entropia)
  zmienne <- zmienne[-wybrane_zmienne$zmienne]
  ####KROK k-ty
  for (k in 1:(kmax-1)){
  mat <- matrix(nrow = length(zmienne)
                ,ncol = length(wybrane_zmienne$zmienne))
  colnames(mat) <- wybrane_zmienne$zmienne
  rownames(mat) <- zmienne
  for (zm in 1:length(zmienne)) {
    for (zm2 in 1:length(wybrane_zmienne$zmienne)){
     mat[zm,zm2] <- condinformation(X = y,Y = X[,zmienne[zm]],
                                    S = X[,as.numeric(
                                      wybrane_zmienne$zmienne[zm2])]
                                    ,method = "emp")
       
    }
  }
  entropia <- apply(mat,MARGIN = 1,min)
  wybrane_zmienne[["zmienne"]] <- c(wybrane_zmienne[["zmienne"]],
                                    as.numeric(names(which.max(entropia))))
  wybrane_zmienne[["wynik"]] <- c(wybrane_zmienne[["wynik"]],
                                  max(entropia))
  zmienne <- c(1:ncol(X))[-as.numeric(wybrane_zmienne$zmienne)]
  }
  wybrane_zmienne$zmienne <- unname(wybrane_zmienne$zmienne)
  wybrane_zmienne
    }
  }


CMIM_sel <- CMIM_selection(X = X,y = y,kmax  = 5)

round(CMIM_sel$wynik,4)

CMIM_sel <- CMIM_selection(X = X,y = y,kmax  = 1)

round(CMIM_sel$wynik,4)

