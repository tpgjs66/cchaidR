
findsplit <- function(response,
                      vars,
                      data,
                      weights,
                      parent = parent,
                      minbucket,
                      alpha_split = 0.05,
                      alpha_merge = 0.05) {

  ## SUBSET the data for recursion
  ## weight is deployed for subsetting original data for recursion purpose.
  ## However, still it is needed to drop parent predictor(s).
  data_for_recursion <- data[as.logical(weights),]
  data_for_recursion <- as.data.frame(data_for_recursion)
  ## For the first split, there is no parent predictor.
  if(is.null(parent)){
    data_for_recursion <- data_for_recursion
    ## After the first split, DROP parent predictor(s) from "data_for_recursion"
  }else{
    data_for_recursion <- data_for_recursion[,-c(parent)]
  }
  data_for_recursion <- as.data.frame(data_for_recursion)
  if (ncol(data_for_recursion) != 1){
  data_for_recursion <- droplevels(data_for_recursion)
  }

  ## MERGE insignificant categories for each predictor
  merged <- ccmerge(response, vars, data_for_recursion, minbucket, alpha_merge)

  ## name of chosen predictor
  chosen_predictor_name <- names((which.min(merged[[1]])))
  ## chosen predictor in integer
  chosen_predictor <- match(chosen_predictor_name,names(data))

  ## CHOOSE the most significant predictor
  xselect <- chosen_predictor
  x <- merged[[2]][[chosen_predictor_name]]

  if(is.null(x)){
    return()
  }

  ## Define splitpoint with the optimally merged category
  splitpoint <- levels(x[drop = TRUE])
  splitpoint <- (strsplit(splitpoint,"-"))
  splitpoint <- lapply(splitpoint, as.integer)

  splitindex <- (1:length(levels(data[,chosen_predictor_name])))
  ## Data-oriented : because all categories start from 0
  splitindex <- splitindex - 1L

  ## splitindex editing
  ## http://stackoverflow.com/questions/11002391/
  ## fast-way-of-getting-index-of-match-in-list
  for (i in 1:length(splitindex)){
    if(length(which(sapply(splitpoint,
                           FUN = function(X) splitindex[i] %in% X)))==0){
      splitindex[i] <- NA
    }else{
      splitindex[i] <- which(sapply(splitpoint,
                                    FUN=function(X) splitindex[i] %in% X))
    }
  }

  return(list(partysplit(varid = as.integer(xselect),
                         index = splitindex,
                         info = (list(p.value=unlist(merged[[1]])))),
              data_for_recursion,
              splitpoint))
}
