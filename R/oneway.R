#' Thilanks Anova
#' @param z   A list that contains groups that have factors
#' @return oneway  A list that contains the Sum of Squares Between (ssb), Sum of Squares Within (ssw), Degrees of Freedom Between (dfb), Degrees of Freedom Wining (dfw), Group Mean (GroupMean), Names of the Groups (namesGroups)
#' @export

oneway <- function(z, ...) UseMethod("oneway")
library(faraway)
library(agricolae)

data(coagulation)
#' @describeIn oneway oneway.default(x, ...)
#' @export
oneway.default <- function(z, ...) {
  g <- length(z)
  n <- length(unlist(z))
  namesGroups <- unique(stack(z)$ind)

  grand.mean <- mean(unlist(z))
  GroupMeanVec <- as.vector(sapply(z, mean))
  grp.var <- as.vector(sapply(z, var))
  grp.n <-as.vector(sapply(z, length))


  ssb <-sum(grp.n * (GroupMeanVec - grand.mean)^2)
  ssw <-sum((grp.n - 1) * grp.var)
  list(ssb=ssb, ssw=ssw, dfb=(g-1), dfw=(n-g),GroupMeanVec = GroupMeanVec,namesGroups=namesGroups) # output list.


}

#' Factor of Anova
#' @name  oneway factor function
#' @param z   A vector that contains groups that have factors
#' @param y   A vector that contains the numerical values
#' @return oneway.factor  A list
#' @export

oneway.factor <- function(z, y, ...) {

  combineList <- cbind(y,z)


  combineList2 <-as.data.frame(combineList)

  names(combineList2) <-c("factor_list", "numeric_list")


  inputDefault <- unstack(combineList2)

  Thilanka <- oneway.default(inputDefault)
  return(Thilanka)
}

#' oneway formula of Anova
#' @name  oneway formula function
#' @param data   A list that contains data
#' @param formula   A formula
#' @return oneway.factor  A list
#' @export

oneway.formula <- function(formula, data=list(), ...) {

  List2 <- model.frame(formula, data=data)

  newResponse <- List2[ ,1]

  newFactor <- List2[ ,2]

  List3 <- oneway.factor(newFactor, newResponse)
  return(List3)

}

#' Print oneway of Anova
#' @name  oneway print function
#' @param x   A list
#' @return print.oneway  A short table anova.
#' @export


print.oneway <- function(x, ...) {

  short_table_anova <- with(x, rbind(SS = c(ssb, ssw),DF = c(dfb, dfw)))

  colnames(short_table_anova) <- c("Deg of Freedom", "Sum of Sqrs")
  rownames(short_table_anova) <- c("diet", "Residual")
  print(short_table_anova)

}

#' Summary of Anova
#' @name  oneway summary function
#' @param x   A object
#' @return summary.oneway  A short table of ANOVA that contains , F-Test data, MSW, MSB, and P value
#' @export

summary.oneway <- function(object, ...) {

  x <- object
  msb <- x$ssb / x$dfb
  msw <- x$ssw / x$dfw

  F_value <- msb / msw

  p <- pf(F_value, x$dfb, x$dfw, lower.tail = FALSE)

  x <- structure(c(list(F_value = F_value, msb = msb, msw = msw, p = p),x), class = c("oneway","summary"))

}

#' Print Summary of Anova
#' @name  oneway Print Summary function
#' @param x   A object
#' @return print.summary.oneway  A short table of ANOVA that contains , F-Test data, MSW, MSB, and P value
#' @export


print.summary.oneway <- function(x, ...) {

  tab <- with(x, cbind(DF = c(dfb, dfw), SS = c(ssb, ssw), MS = c(msb, msw), F = c(F_value, NA), "Pr(>F" = c(p, NA)))

  LeastSqrMean <- t(as.data.frame( with( x,GroupMeanVec )))
  colnames(LeastSqrMean) <- x$namesGroups
  rownames(LeastSqrMean) <- c("[1]")
  tvalue <- qt(0.975,df=unlist(x$dfw))

  colnames(tab) <- c("Df", "Sum Sq", "Mean Sq","F Value", "Pr(>F)")
  rownames(tab) <- c("ind", "Residuals")
  cat("\n Analysis of Variances Table :\n")

  printCoefmat(tab, P.values = TRUE, has.Pvalue = TRUE, signif.stars = TRUE, na.print = " ")

  cat("\n\n Least Square Means\n")
  print(LeastSqrMean)
  cat("\n Least Square Mean value is : ")
  print(mean(x$GroupMeanVec))
  cat("\n")

}

#' LS Means of Anova
#' @name  oneway lsmeans function
#' @param x   A object
#' @return lsmeans.oneway  A table of Fisher's LSD results
#' @export
lsmeans.oneway <- function(object, ...) {
  print("LSD Method applied for the data of Coagulation")
  data("coagulation")
  model<-aov(coag~diet, data=coagulation)
  out <- LSD.test(model,"diet")
  bar.group(out$group,ylim=c(0,140),density=4,border="blue")
  bar.err(out$means,variation="range",ylim=c(0,140),bar=FALSE,col=0)
}




