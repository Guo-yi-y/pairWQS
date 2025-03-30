utils::globalVariables(c("final.weights", "weights_avg", "n_covars", "time"))

#' @import dplyr
#' @import glue
#' @import survival
#' @import Rsolnp
#' @import magrittr
#' @import magrittr
#' @importFrom stats as.formula
#' @importFrom survival Surv
#' @importFrom utils globalVariables

NULL

#' Weighted Quantile Sum Regression without Bootstrapping
#'
#' @param wqsdata The input data
#' @param col_vars Exposure variable names
#' @param col_covars Covariates
#' @param id Matching ID variable
#' @param event Event indicator
#' @param q Number of quantiles
#' @return A list containing WQS weights and the regression coefficient of WQS index
#' @export


pairwqs_noboot <- function(wqsdata, col_vars, col_covars, id = "studyid", event = "event", q=4){
  wholedata <- wqsdata %>%
    mutate(
      across(
        all_of(col_vars),
        ~ as.numeric(
          cut(
            .,
            breaks = quantile(., probs = 0:q/q, na.rm = TRUE),
            include.lowest = TRUE,
            labels = 1:q
          )
        ),
        .names = "{.col}"
      )
    )

  wholedata$id = wholedata[[id]]
  wholedata$event= wholedata[[event]]

  design_mat = as.matrix(wholedata[col_vars])

  if (!is.null(col_covars)) {
    covar_mat = as.matrix(wholedata[col_covars])
    n_covars = length(col_covars)
  } else {
    covar_mat = NULL
    n_covars = 0
  }

  index.cont = which(wholedata$event==0);index.case = which(wholedata$event ==1)
  n_exps = length(col_vars)
  n_covars = length(col_covars)

  if (n_covars == 0){
    fn=function(x){
      -sum((x[1] + x[2] * (design_mat[index.case,] %*% x[3:(2 + n_exps)])-
              log(exp(x[1] + x[2] * (design_mat[index.case,] %*% x[3:(2 + n_exps)])) +
                    sum(exp(
                      x[1] + x[2] * (design_mat[index.cont,] %*% x[3:(2 + n_exps)])))
              )
      ))
    }
  } else {
    fn=function(x){
      -sum((x[1] + x[2] * (design_mat[index.case,] %*% x[3:(2 + n_exps)]) + (covar_mat[index.case,] %*% x[(3 + n_exps):(2+n_exps+n_covars)])-
              log(exp(x[1] + x[2] * (design_mat[index.case,] %*% x[3:(2 + n_exps)]) + (covar_mat[index.case,] %*% x[(3 + n_exps):(2+n_exps+n_covars)])) +
                    sum(exp(
                      x[1] + x[2] * (design_mat[index.cont,] %*% x[3:(2 + n_exps)]) + (covar_mat[index.cont,] %*% x[(3 + n_exps):(2+n_exps+n_covars)])))
              )
      ))
    }
  }


  eqn1=function(x){z1=sum(x[3:(2 + n_exps)]); return(z1)}
  eqn2=function(x){vec = x[3:(2 + n_exps)];return(vec)}
  x0 = c(rep(0.5,2 + n_exps), rep(0, n_covars))#initial values
  pars = solnp(x0, fun = fn, eqfun = eqn1 , eqB = 1, ineqfun = eqn2,
               ineqLB = rep(0.001,n_exps), ineqUB = rep(0.999,n_exps), control = list(maxiter = 1000, tol = 1e-10, trace = 1))



  wholedata$wqs = design_mat %*% pars$pars[3: (2 + n_exps) ]

  if (n_covars == 0){
    fit = coxph(Surv(time, event) ~ wqs +strata(id), data = wholedata)
  } else{
    fit = coxph(as.formula( glue("Surv(time, event) ~ wqs+{paste(col_vars, collapse = '+')}+strata(id)") ), data = wholedata)
  }




  newlist<-list (vars = col_vars,
                 "final.weights" = pars$pars[3:(2 + n_exps)], "wqs_beta" = summary(fit)$coefficient)
  return(newlist)
}


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

    return(data.frame(t(pairwqs_noboot(sample_n(wqsdata, nrow(wqsdata), replace = T), col_vars, col_covars,id, event, q)$final.weights)))

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

  if (is.null(valid_data)) {
    valid_data <- train_data
  }

  if(boot == FALSE){
    train_res = pairwqs_noboot(train_data, col_vars, col_covars, id, event, q)
  } else {
    train_res = pairwqs_boot(train_data, col_vars, col_covars, id, event, q, B)
  }

  wholedata = valid_data
  wholedata$id = wholedata[[id]]
  wholedata$event= wholedata[[event]]

  weights = train_res$final.weights

  valid_design_mat = as.matrix(wholedata[col_vars])

  wholedata$wqs = valid_design_mat %*% final.weights

  if (length(col_covars) == 0){
    fit = coxph(Surv(time, event) ~ wqs +strata(id), data = wholedata)
  } else{
    fit = coxph(as.formula( glue("Surv(time, event) ~ wqs+{paste(col_vars, collapse = '+')}+strata(id)") ), data = wholedata)
  }

  newlist<-list (vars = col_vars,
                 "final.weights" = weights_avg, valid_mod = fit)
  return(newlist)

}
