
growtree <- function(id = 1L,
                     response,
                     vars,
                     data,
                     weights,
                     parent,
                     minbucket,
                     minsplit,
                     alpha_split,
                     alpha_merge,
                     max_depth) {

  ## for less than MINBUCKET, stop here
  ## !!CHECK MINBUCKET VS MINSPLIT!!
  if (sum(weights) < minsplit) return(partynode(id = id))

  ## stop when max_depth is reached
  if (max_depth == length(parent) ) return(partynode(id = id))

  ## Find best split ###########
  sp <- findsplit(response,
                  vars,
                  data,
                  weights,
                  parent = parent,
                  minbucket,
                  alpha_split,
                  alpha_merge)

  ## no split found, stop here
  if (is.null(sp[[1]])) {return(partynode(id = id))}
  else{
    sp[[2]] <- droplevels(sp[[2]])
  }

  ## NEEDED to be fixed: GOTO findsplit
  ## stop when min_bucket is reached
  split_bucket <- c()
  for(i in 1:length(sp[[3]])){
    bucket <- table(sp[[2]][names(data[sp[[1]][[1]]])])
    split_bucket[i]<- sum(bucket[as.character(unlist(sp[[3]][[i]]))])
  }

  if (any(split_bucket < minbucket)) return(partynode(id = id))

  ## if alpha_split is reached, stop here
  if (min(unlist(sp[[1]][[6]])) > alpha_split) return(partynode(id = id))

  ## actually split the data
  kidids <- kidids_split(sp[[1]], data = data)
  ## if(any(is.na(kidids))){warning("Check split index in partysplit")}

  ## set up all daugther nodes
  kids <- vector(mode = "list", length = max(kidids, na.rm = TRUE))


  for (kidid in 1:length(kids)) {
    ## select observations for current node
    w <- weights
    w[kidids != kidid] <- 0

    ## get next node id
    if (kidid > 1) {
      myid <- max(nodeids(kids[[kidid - 1]]))
      parent <- parent

    } else {
      myid <- id
      parent <- c(parent,sp[[1]][[1]])
    }

    ## start recursion on this daugther node
    kids[[kidid]] <- growtree(id = as.integer(myid + 1),
                              response,
                              vars,
                              data,
                              w,
                              parent,
                              minbucket,
                              minsplit,
                              alpha_split,
                              alpha_merge,
                              max_depth)
  }

  ## return nodes
  return(partynode(id = as.integer(id),
                   split = sp[[1]],
                   kids = kids,
                   info = list(p.value = min(info_split(sp[[1]])$p.value,
                                             na.rm = TRUE))))
}
