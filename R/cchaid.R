#' @title Chi-squared Automated Interaction Detection for continuous response
#'  variable and tree
#'
#' @description Fits a classification tree by the CHAID algorithm for
#' continuous response variable
#'
#' @param formula an object of class formula (or one that can be coerced to
#' that class): a symbolic description of the model to be fitted. Response
#' variable should be continuous and all predictors should be categorical
#' (either ordered or not).
#' @param weights an optional vector of weights to be used in the fitting
#' process. Should be NULL or a numeric vector.
#' @param minbucket Minimum number of observations in terminal nodes.
#' @param minsplit Number of observations in splitted response at which
#' no further split is desired.
#' @param alpha_split Level of significance used for splitting of a node
#' in the most significant predictor
#' @param alpha_merge Level of significance used for merging of predictor
#' categories
#' @param max_depth Maximum depths for the tree
#' @details The CHAID algorithm is originally proposed by Kass (1980) which
#' allow multiple splits of a node. The current implementation only accepts
#' continuous response variable and categorical predictors (nominal or ordinal).
#' If response variable is categorical, refer to \code{\link[CHAID]{chaid}}.
#' CHAID consist of three steps: merging, splitting and stopping.
#' A tree is grown by repeatedly using the three steps below on each
#' node starting from the root node.
#'\itemize{
#'\item Merging: For each predictor variable X, merge non-significant
#' categories. Each final category of X will result in one child node if X is
#' used to split the node. The merging step also calculates the adjusted p-value
#' that is to be used in the splitting step.
#' \itemize{
#' \item 1. If X has 1 category only, stop and set the adjusted p-value to be 1.
#' \item 2. If X has 2 categories, go to step 8.
#' \item 3. Else, find the allowable pair of categories of X (an allowable pair
#' of categories for ordinal predictor is two adjacent categories, and for
#' nominal predictor is any two categories) that is least significantly
#' different (i.e., most similar). The most similar pair is pair whose test
#' statistic gives the largest p-value with respect to the dependent variable Y.
#' \item 4. For the pair having the largest p-value, check if its p-value is
#' larger than a user-specified alpha-level α merge (alpha_merge). If it does,
#' this pair is merged into a single compound category. Then a new set of
#' categories of X is formed. If it does not, then go to step 7.
#' \item 5. (Optional) If the newly formed compound category consists of three
#' or more original categories, then find the best binary split within the
#' compound category which p-value is the smallest. Perform this binary split
#' if its p-value is not larger than an alpha-level split-merge α
#' (alpha_spli-merge).
#' \item 6. Go to step 2.
#' \item 7. (Optional) Any category having too few observations (as compared
#' with a user-specified minimum segment size) is merged with the most similar
#' other category as measured by the largest of the p-values.
#' \item 8. The adjusted p-value is "mandatorily" computed for the merged
#' categories by applying Bonferroni adjustments.
#'}
#'\item Splitting: The “best” split for each predictor is found in the merging
#' step. The splitting step selects which predictor to be used to best split
#' the node. Selection is accomplished by comparing the adjusted p-value
#' associated with each predictor. The adjusted p-value is obtained in the
#' merging step.
#' \itemize{
#'  \item 1. Select the predictor that has the smallest adjusted p-value
#'  (i.e., most significant).
#'  \item 2. If this adjusted p-value is less than or equal to a user-specified
#'  alpha-level split α (alpha_split), split the node using this predictor.
#'  Else, do not split and the node is considered as a terminal node.
#' }
#'\item Stopping: The stopping step checks if the tree growing process should
#' be stopped according to the following stopping rules.
#'  \itemize{
#'\item 1. If a node becomes pure; that is, all cases in a node have identical
#' values of the dependent variable, the node will not be split.
#'\item 2. If all cases in a node have identical values for each predictor,
#' the node will not be split.
#'\item 3. If the current tree depth reaches the user specified maximum tree
#' depth limit value, the tree growing process will stop.
#'\item 4. If the size of a node is less than the user-specified minimum node
#' size value, the node will not be split.
#'\item 5. If the split of a node results in a child node whose node size is
#' less than the userspecified minimum child node size value, child nodes that
#' have too few cases (as compared with this minimum) will merge with the most
#' similar child node as measured by the largest of the p-values. However, if
#' the resulting number of child nodes is 1, the node will not be split.
#' }
#' }
#' @return An object of class constparty, see package
#'  \code{\link[partykit]{party}}.
#' @export
#' @references Kass, G. V. (1980). An exploratory technique for investigating
#' large quantities of categorical data. Applied statistics, 119-127.\cr
#' Hothorn T, Zeileis A (2015). partykit: A Modular Toolkit for
#' Recursive Partytioning in R. Journal of Machine Learning Research,
#' 16, 3905–3909.
#' @seealso \code{\link[CHAID]{chaid}} \code{\link[partykit]{ctree}}
#' \code{\link[partykit]{glmtree}}
#'
#' @examples
#' library("cchaid")
#' require("partykit")
#'
#' set.seed(100)
#' WorkDurationS <- WorkDuration[sample(1:nrow(WorkDuration), 1000),]
#' formula <- (Dur ~ Urb + Comp + Child + Day + pAge + SEC + Ncar + Gend +
#'             Driver + wstat + Pwstat + Xdag + Xn.dag + Xarb + Xpop + Ddag +
#'             Dn.dag + Darb + Dpop)
#' mytree <- cchaid(formula,data = WorkDurationS, weights = NULL, minbucket = 57,
#'                  minsplit = 114, alpha_split=0.05, alpha_merge=0.05,
#'                  max_depth = 8)
#' mytree
#' plot(mytree)



##--------------------------------------------------------------------------##
#                                                                            #
#                       MAIN FUNCTION CALLED BY THE USER                     #
#                                                                            #
##--------------------------------------------------------------------------##

cchaid <- function(formula,
                   data,
                   weights = NULL,
                   minbucket = 100,
                   minsplit = 200,
                   alpha_split = 0.05,
                   alpha_merge = 0.05,
                   max_depth = 10) {

  ## max_depth check
  if(isTRUE(max_depth == -1) | is.null(max_depth)){
    max_depth <- .Machine$integer.max
    }
    
  ## name of the response variable and independent variables
  allvars<- all.vars(formula)
  response <- all.vars(formula[[2]])
  vars<- all.vars(formula[[3]])

  ## Drop the variable not specifieds in formula
  data <- subset(data,select=allvars)

  ## data without missing values, response comes last
  #data <- data[complete.cases(data), c(all.vars(formula)[-1], response)]

  ## data is factors only
  ## stopifnot(all(sapply(data, is.factor)))

  if(is.null(weights)) weights <- rep(1L, nrow(data))
  ## weights are case weights, i.e., integers
  stopifnot(length(weights) == nrow(data) &
              max(abs(weights - floor(weights))) < .Machine$double.eps)

  ## Initialize a list to store parent nodes
  parent <- c()

  ## grow tree
  nodes <- growtree(id = 1L,
                    response,
                    vars,
                    data,
                    weights,
                    parent,
                    minbucket,
                    minsplit,
                    alpha_split,
                    alpha_merge,
                    max_depth)

  ## compute terminal node number for each observation
  fitted <- fitted_node(nodes, data = data)
  ## return rich constparty object
  ret <- party(nodes, data = data,
               fitted = data.frame("(fitted)" = fitted,
                                   "(response)" = data[[response]],
                                   "(weights)" = weights,
                                   check.names = FALSE),
               terms = terms(formula))
  as.constparty(ret)
}

## minbucket calculation
#
# int detmincas(int d,int n)
# {  float y;
#
#   if (d==0) { //discreet
#     y=(float)n/150;
#     if (y<75) y=75;
#     if (y>150) y=150;
#   }
#   else {   //continue
#     y=(float)n/200;
#     if (y<50) y=50;
#     if (y>100) y=100;
#   }
#   return (int)y;
# }

#data<-WorkDuration

#formula <- (Dur ~ Urb + Comp + Child + Day + pAge + SEC + Ncar + Gend + Driver + wstat + Pwstat + Xdag + Xn.dag + Xarb + Xpop + Ddag + Dn.dag + Darb + Dpop)
#formula <- (Dur.ratio ~ Urb + Comp + Child + Day + pAge + SEC + Ncar + Gend + Driver + wstat + Pwstat + Xdag + Xn.dag + Xarb + Xpop + Ddag + Dn.dag + Darb + Dpop + Dur)
#formula <- (Stfix ~ Urb + Comp + Child + Day + pAge + SEC + Ncar + Gend + Driver + wstat + Pwstat + Xdag + Xn.dag + Xarb + Xpop + Ddag + Dn.dag + Darb + Dpop + Act + yWo + Dur + yNep + Ratio + Inter + BT + T1 + T2 + T3 + T4 + T5 + T6 + T7 + T8)
#formula <- (Dur.fix ~ Urb + Comp + Child + Day + pAge + SEC + Ncar + Gend + Driver + wstat + Pwstat + Xdag + Xn.dag + Xarb + Xpop + Ddag + Dn.dag + Darb + Dpop + Act + yWo + Dur + yNep + Ratio + Inter + BT + T1 + T2 + T3 + T4 + T5 + T6 + T7 + T8)

#mytree <- cchaid(formula, data = data, weights = NULL, minbucket = 57, minsplit = 114, alpha_split=0.05, alpha_merge=0.05, max_depth = 8)
#mytree <- cchaid(formula, data = aidwork4, weights = NULL, minbucket = 50, minsplit = 100, alpha_split=0.05, alpha_merge=0.05, max_depth = 8)
#mytree <- cchaid(formula, data = actsec3, weights = NULL, minbucket = 50, minsplit = 100, alpha_split=0.05, alpha_merge=0.05, max_depth = 8)
#mytree <- cchaid(formula, data = actsec5, weights = NULL, minbucket = 50, minsplit = 100, alpha_split=0.05, alpha_merge=0.05, max_depth = 8)



