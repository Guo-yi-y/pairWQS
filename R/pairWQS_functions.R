utils::globalVariables(c("final.weights", "weights_avg", "n_covars", "time"))

#' @import dplyr
#' @import reticulate
#' @import glue
#' @import survival
#' @import Rsolnp
#' @import magrittr
#' @import magrittr
#' @importFrom stats as.formula
#' @importFrom survival Surv
#' @importFrom utils globalVariables

NULL

#' Weighted Quantile Sum Regression with Bootstrapping
#'
#' @param wqsdata The input data
#' @param col_vars Exposure variable names
#' @param col_covars Covariates
#' @param id Matching ID variable
#' @param event Event indicator
#' @param q Number of quantiles
#' @param B Number of bootstrap
#' @return A list containing WQS weights and the regression coefficient of WQS index
#' @export

pairwqs_boot = function(wqsdata, col_vars, col_covars, id = "studyid", event = "event", q=4, B=10){

  res_list = lapply(1:B, function(i){

    return(data.frame(t(py$pairwqs_noboot(sample_n(wqsdata, nrow(wqsdata), replace = T), col_vars, col_covars,id, event, q)$final.weights)))

  })
  weights_avg = res_list %>% bind_rows() %>% colMeans()

  wholedata = wqsdata
  wholedata$id = wholedata[[id]]
  wholedata$event= wholedata[[event]]

  design_mat = as.matrix(wholedata[col_vars])

  wholedata$wqs = design_mat %*% weights_avg

  if (length(col_covars) == 0){
    fit = coxph(Surv(time, event) ~ wqs +strata(id), data = wholedata)
  } else{
    fit = coxph(as.formula( glue("Surv(time, event) ~ wqs+{paste(col_vars, collapse = '+')}+strata(id)") ), data = wholedata)
  }

  newlist<-list (vars = col_vars,
                 "final.weights" = weights_avg, "wqs_beta" = summary(fit)$coefficient)

  return(newlist)
}


#' Weighted Quantile Sum Regression for design with paired data structure
#'
#' @param train_data The training data
#' @param valid_data The validation data
#' @param col_vars Exposure variable names
#' @param col_covars Covariates
#' @param id Matching ID variable
#' @param event Event indicator
#' @param q Number of quantiles
#' @param boot Whether to perform bootstrap
#' @param B Number of bootstrap
#' @return A list containing WQS weights and the regression coefficient of WQS index
#' @export
#'
#'
pairwqs = function(train_data, valid_data = NULL, col_vars, col_covars, id = "studyid", event = "event", q=4, boot = FALSE, B=10){

  q = as.integer(q)
  if (is.null(valid_data)) {
    valid_data <- train_data
  }

  if(boot == FALSE){
    train_res = py$pairwqs_noboot(train_data, col_vars, col_covars, id, event, q)
  } else {
    train_res = pairwqs_boot(train_data, col_vars, col_covars, id, event, q, B)
  }

  wholedata = valid_data
  wholedata$id = wholedata[[id]]
  wholedata$event= wholedata[[event]]

  weights = train_res$final.weights

  wqs_beta = train_res$wqs_beta

  valid_design_mat = as.matrix(wholedata[col_vars])

  wholedata$wqs = valid_design_mat %*% weights

  if (length(col_covars) == 0){
    fit = coxph(Surv(time, event) ~ wqs +strata(id), data = wholedata)
  } else{
    fit = coxph(as.formula( glue("Surv(time, event) ~ wqs+{paste(col_covars, collapse = '+')}+strata(id)") ), data = wholedata)
  }

  newlist<-list (vars = col_vars,
                 wqs_beta = wqs_beta,
                 "final.weights" = weights, valid_mod = fit)
  return(newlist)

}
