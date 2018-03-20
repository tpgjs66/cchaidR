## This is ccmerge
ccmerge <- function(response, vars, data, minbucket, alpha_merge = 0.05) {

  ## Check data is dataframe and alpha_merge <= 1
  if (!is.data.frame(data))
    stop("data is not a dataframe", call. = TRUE)
  if (alpha_merge > 1)
    warning("p value provided for alpha_merge should be <=1", call. = TRUE)

  ## Detect each column on its data_type -> function to be used is class;
  ## and send it to its appropriate function for merging loop
    #n <- colnames(data[, -c(ncol(data)),drop=F])
    n <- colnames(data[, -which(names(data) %in% response), drop = F])
    #t <- vars
    #n <- vars

  ## Initialization of lists to store adjusted p-value and merged category data
  p_adj_list<-list()
  merged<-list()

  ## Initialization of a list to get the result of merging_loop
  merging_result<-vector("list")

  ## Start merging loop for each predictor i
  for (i in n) {

    # Get adj. p-values and merged data[,i] with Bonferroni correction,
    # which is an output of "merging_loop" function
    merging_result[[i]]<-merging_loop(response, data, data[,i], i, minbucket, alpha_merge)

    # A list to store p_adj for respective predictor i
    p_adj_list[i]<-merging_result[[i]][1]

    # A list to store merged category for respective predictor i
    merged[i]<-merging_result[[i]][2]
  }

  # Binding merged category data with response variable
  # NOTE: Response variable should be placed at the end column of dataframe.
  merged<-cbind(merged,data[length(data)])

  return(list(p_adj_list,merged))
}


##--------------------------------------------------------------------------##
# FUNCTION BELOW MERGES INSIGNIFICANT CATEGORIES WITHIN EACH nominal COLUMN
# data[,i] OF DATAFRAME data AND RETURNS ADJ. P-VALUE list and merged categories
##--------------------------------------------------------------------------##


merging_loop <- function(response, data, y, i, minbucket, alpha_merge) {

  # Number of categories for predictor i
  l = length(levels(data[,i]))
  merged_category<-vector("list")

  ## If X has 1 category only (i.e., pure),
  ## Stop merging process and set the adjusted p-value to be 1
  if (l == 1){

    ## The results of merging_loop: p_aov, merged_category
    ## These lines of code always have to stick together
    ## NOTE: Do not need bonferroni correction when there is no merging
    p_aov <- 1
    merged_category <- data[,i]

  ## If X has 2 categories, adjusted p-value is computed for the category
  ## by applying Bonferroni adjustments (c=2, r={1,2}).
  }else if (l == 2){

    ## Check minbucket here!!
    ## If the split of a node results in a child node whose node size is
    ## less than "minbucket", corresponding child nodes will be merged with
    ## the most similar child node as measured by the largest of the p-values.
    ## NOTE: If X has only 2 categories, as in this case, child nodes will be
    ## merged each other and results single merged category.
    if (any(table(data[,i]) < minbucket)){
      r <- names(table(data[,i])[2])
      c <- names(table(data[,i])[1])
      nameofMergedCategory <- paste(r, c, sep = "-")
      data[,i] <- merging_function(data, y, r, c, nameofMergedCategory, i)

      p_aov <- 1
      merged_category <- data[,i]

      ## Bonferroni correction
      number_original_category <- 2
      number_merged_category <- 1
      b <- Bonferroni_correction(type = class(data[,i])[1],
                                 o = number_original_category,
                                 m = number_merged_category)
      if (p_aov * b > 1) {
        p_aov = 1
      } else {
        p_aov = p_aov * b
      }

    }else{

    ## The results of merging_loop: p_aov, merged_category
    ## These lines of code always have to stick together
    p_aov = as.matrix(summary((aov(data[, response]~data[,i]
                                   ,data)))[[1]][,5])[1]
    merged_category <- data[,i]
    }

  ## If X has more than 2 categories, find the allowable pair of categories of X
  ## (an allowable pair of categories for ordinal predictor is two adjacent
  ## categories, and for nominal predictor is any two categories) that is
  ## least significantly different (i.e., most similar).
  ## The most similar pair is the pair whose test statistic gives the
  ## largest p-value with respect to the response variable Y.
  }else{
    k<-c(0)    # index for merging loop

    ## If X has more than two categories with a single observeation,
    ## which prevents to do pairwise t-test,
    ## merge one of them with the most similar category

    if(length(table(data[,i])[table(data[,i])==1]) >= 1) {
      category_single_obs <- names(which(table(data[,i]) == 1))

      ## NOTE: An allowable pair of categories for "nominal" predictor is any
      ##       categories.
      if(class(data[,i])[1]=="factor"){
        n <- length(category_single_obs)
        for (x in n:2) {
          r <- category_single_obs[x]
          r_mean <- mean(data[,response][data[,i]== r])

          c_list <- names(table(data[,i]))[!names(table(data[,i])) %in% r]
          c_list_mean <- c()

            for (z in 1:length(c_list)){
              c_list_mean[z] <- mean(data[,response][data[,i]== c_list[z]])
            }

          names(c_list_mean)<- c_list
          c <- names(which.min(abs(c_list_mean - r_mean)))
          nameofMergedCategory <- paste(r, c, sep = "-")
          data[,i] <- merging_function(data, y, r, c, nameofMergedCategory, i)

          if(length(table(data[,i])[table(data[,i])==1]) <= 1){break}
        }
      ## NOTE: An allowable pair of categories for "ordered" predictor is two
      ##       adjacent categories.
      } else if (class(data[,i])[1]=="ordered") {
        n <- length(category_single_obs)
        for (x in n:1) {
          r <- category_single_obs[x]
          r_mean <- mean(data[,response][data[,i]== r])
          c_list <- names(table(data[,i]))[!names(table(data[,i])) %in% r]

          ## "NAs introduced by coercion" warning message may appear here.
          ## Find adjacent categories (tt==1)
          suppressWarnings(
            tt <- abs(as.numeric(names(table(data[,i]))
                                 [!names(table(data[,i])) %in% r])
                      - as.numeric(r))
            )

          names(tt) <- c_list

          ## A list only contains adjacent categories
          c_list_adjacent <- names(tt[tt %in% 1])

          c_list_adjacent_mean <- c()

          ## Calculate mean difference to examine similarity
            if(length(c_list_adjacent) >= 2) {
              for (z in 1:length(c_list)){
                c_list_adjacent_mean[z] <-
                  mean(data[,response][data[,i]== c_list_adjacent[z]])
              }
              names(c_list_adjacent_mean)<- c_list_adjacent
              c <- names(which.min(abs(c_list_adjacent_mean - r_mean)))

            }else if(length(c_list_adjacent) == 1){
              c <- c_list_adjacent

            }else {
              next
            }

          nameofMergedCategory <- paste(r, c, sep = "-")
          data[,i] <- merging_function(data, y, r, c, nameofMergedCategory, i)

          ## Break if there is less than 2 category with a single observation
          if(length(table(data[,i])[table(data[,i]) == 1]) < 2){break}
        }
          ## Though there are more than 2 categories with single observation,
          ## if the function cannot find allowable pair of categories to be
          ## merged for ordered predictor,
          ## there will be no merging and return p-value as 1.
          if(!exists("nameofMergedCategory")){

              p_aov <- 1
              merged_category <- data[,i]

              return(list(p_aov,merged_category))
          }
      } else {warning("Variable is not factor type")}
    }

      repeat{
      k <- k+1  # ADD  1 to loop index
      l <- length(levels(data[,i])) # Number of categories for predictor i
      
      ## Break if there is a single observation on category  
      if(length(table(data[,i])[table(data[,i])==1])!=0){break}  
        
      ## Pairwise t.test to get p.value matrix
      p =(pairwise.t.test(data[, response],
                         data[,i],p.adjust.method = "none",
                         paired = FALSE, pool.sd = FALSE,
                         var.equal = TRUE))$p.value

      if(any(dim(p)==0)) {break}

      if(class(data[,i])[1]=="factor"){
        ## This picks the max p value for nominal predictor,
        ## which means the least significant pair.
        p_max_value = max(p, na.rm = TRUE)
      }else if(class(data[,i])[1]=="ordered"){
        ## This picks the max p value ordinal predictor,
        ## which means the least significant pair.
        p_max_value = max(diag(p), na.rm = TRUE)
      }else {
        warning("Variable is not factor type")
      }

      p_max = which(p == max(p_max_value, na.rm = TRUE), arr.ind = TRUE)
      r <- rownames(p)[p_max[, 1]]
      c <- colnames(p)[p_max[, 2]]

      ## Do Wilcox test if Normality assumption is violated (i.e., # obs. < 30)
      if((length(data[data[,i]==r,response]) < 30 |
          length(data[data[,i]==c,response]) < 30)){

        wilcox_p_max_value <-
          round(wilcox.test(data[data[,i]==r,response],
                            data[data[,i]==c,response],
                            exact=FALSE)$p.value,digits=2)

        ## Break the loop if stopping criteria is reached
        if(l==2 | (wilcox_p_max_value < alpha_merge)) {break}

        #------- MERGED HERE WITH A HYPHEN
        else(nameofMergedCategory <- paste(r, c, sep = "-"))

        # CALL MERGING FUNCTION here for actual merging
        data[,i]<-merging_function(data, y, r, c, nameofMergedCategory, i)

      ## Break the loop if stopping criteria is reached
      } else if(l==2 | p_max_value < alpha_merge) {break}

      ## Merging categories if both category obs. > 30
      else{
        p_max = which(p == max(p_max_value, na.rm = TRUE), arr.ind = TRUE)
        r <- rownames(p)[p_max[, 1]]
        c <- colnames(p)[p_max[, 2]]

        #------- MERGED HERE WITH A HYPHEN
        nameofMergedCategory <- paste(r, c, sep = "-")

        # CALL MERGING FUNCTION here for actual merging
        data[,i] <- merging_function(data, y, r, c, nameofMergedCategory, i)

      }

    }

    ##### This lines of code always have to stick together
    p_aov =
      as.matrix(summary((aov(data[, response]~data[,i],data)))[[1]][,5])[1]
    merged_category <- data[,i]

    ## Bonferroni correction
    number_original_category <- length(levels(y))
    number_merged_category <- length(levels(data[,i]))
    b <- Bonferroni_correction(type = class(data[,i])[1],
                               o = number_original_category,
                               m = number_merged_category)
    if (p_aov * b > 1) {
      p_aov = 1
    } else {
      p_aov = p_aov * b
    }

  ## Check minbucket here
  ## Any category having too few observations (as compared with a user-specified
  ## minimum segment size) is merged with the most similar other category
  ## as measured by the largest of the p-values.
  while(any(table(data[,i]) < minbucket)){

    if(any(table(data[,i]) == 1)){
      p_aov = 1
      merged_category <- data[,i]

      break

    } else if(length(levels(data[,i])) == 2){
      c <- names(table(data[,i])[1])
      r <- names(table(data[,i])[2])

      #------- MERGED HERE WITH A HYPHEN
      nameofMergedCategory <- paste(r, c, sep = "-")

      # CALL MERGING FUNCTION here for actual merging
      data[,i] <- merging_function(data, y, r, c, nameofMergedCategory, i)

      # Set p_aov as 1 and stop merging
      p_aov <- 1
      merged_category <- data[,i]

    }else {

      m <- names(which.min(table(data[,i])))

      if (!any((c(colnames(p),rownames(p)))==m)) {
        ## Pairwise t.test to get p.value matrix
        p = (pairwise.t.test(data[, response],
                             data[,i],
                             p.adjust.method = "none",
                             paired = FALSE,
                             pool.sd = FALSE,
                             var.equal = TRUE))$p.value
      }

      if (any(colnames(p) == m) & any(rownames(p) == m)) {

        if (class(data[,i])[1]=="factor"){
          p_max_value = max(p[m,],p[,m], na.rm=TRUE)
          p_max = which(p == max(p_max_value, na.rm = TRUE), arr.ind = TRUE)
          r <- rownames(p)[p_max[, 1]]
          c <- colnames(p)[p_max[, 2]]

        }else if (class(data[,i])[1]=="ordered") {
          x <- c(p[m,],p[,m])
          x2 <- x[diag(p) %in% x]

          p_max_value = max(x2, na.rm=TRUE)
          p_max = which(p == max(p_max_value, na.rm = TRUE), arr.ind = TRUE)
          r <- rownames(p)[p_max[, 1]]
          c <- colnames(p)[p_max[, 2]]

        }else {warning("Variable is not factor type")}

      } else if (m %in% rownames(p)){

        if (class(data[,i])[1]=="factor"){
          p_max_value = max(p[m,], na.rm=TRUE)
          p_max = which(p == max(p_max_value, na.rm = TRUE), arr.ind = TRUE)
          r <- rownames(p)[p_max[, 1]]
          c <- colnames(p)[p_max[, 2]]

        }else if (class(data[,i])[1]=="ordered") {
          x <- c(p[m,])
          x2 <- x[diag(p) %in% x]

          p_max_value = max(x2, na.rm=TRUE)
          p_max = which(p == max(p_max_value, na.rm = TRUE), arr.ind = TRUE)
          r <- rownames(p)[p_max[, 1]]
          c <- colnames(p)[p_max[, 2]]

        }else {warning("Variable is not factor type")}

      } else {

        if (class(data[,i])[1]=="factor"){
          p_max_value = max(p[,m], na.rm=TRUE)
          p_max = which(p == max(p_max_value, na.rm = TRUE), arr.ind = TRUE)
          r <- rownames(p)[p_max[, 1]]
          c <- colnames(p)[p_max[, 2]]

        }else if (class(data[,i])[1]=="ordered") {
          x <- c(p[,m])
          x2 <- x[diag(p) %in% x]

          p_max_value = max(x2, na.rm=TRUE)
          p_max = which(p == max(p_max_value, na.rm = TRUE), arr.ind = TRUE)
          r <- rownames(p)[p_max[, 1]]
          c <- colnames(p)[p_max[, 2]]

        }else {warning("Variable is not factor type")}
       }

        #------- MERGED HERE WITH A HYPHEN
        nameofMergedCategory <- paste(r, c, sep = "-")
        # CALL MERGING FUNCTION here for actual merging
        data[,i] <- merging_function(data, y, r, c, nameofMergedCategory, i)

        p_aov =
          as.matrix(summary((aov(data[, response]~data[,i],
                                 data)))[[1]][,5])[1]
        merged_category <- data[,i]

        ## Bonferroni correction
        number_original_category <- length(levels(y))
        number_merged_category <- length(levels(data[,i]))
        b <- Bonferroni_correction(type = class(data[,i])[1],
                                   o = number_original_category,
                                   m = number_merged_category)
        if (p_aov * b > 1) {
          p_aov = 1
        } else {
          p_aov = p_aov * b
        }
    }
  }

}
  return(list(p_aov,merged_category))
}


merging_function <- function(data, y, r, c, nameofMergedCategory, i) {

  list_of_levels = levels(data[,i])

  for (q in seq_along(list_of_levels)) {
    if (list_of_levels[q] == r || list_of_levels[q] == c) {
      list_of_levels[q] = nameofMergedCategory
    }
  }
  levels(data[,i]) <- list_of_levels
  return(data[,i])
}

Bonferroni_correction <- function(type,o,m){
  if(type=="ordered"){
    b <- factorial(o-1)/((factorial(m-1)*factorial((o-1)-(m-1))))

  }else if(type=="factor"){
    b <- 0
    for(i in 0:(m-1)){
      b <- b + (-1)^(i)*(((m-i)^(o))/(factorial(i)*factorial(m-i)))
    }
  }
  return(b)
}

