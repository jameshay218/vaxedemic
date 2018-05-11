#' linear function of vaccine production
#' 
#' no vaccine produced until time vax_production_params[["detection_delay"]] + vax_production_params[["production_delay"]], then constant production rate until max number of doses ever made reached, then no production
#' @param vax_production_params list of parameters for vaccine production
#' @param t scalar time
#' @return needs to return a scalar: the number of vaccines ever produced up to time t
#' @export

produce_vax_linear_with_delay <- function(vax_production_params, t) {
  t_since_production <- t - (vax_production_params[["detection_delay"]] + 
                               vax_production_params[["production_delay"]])
  if(t_since_production < 0) {
    0
  } else {
    min(vax_production_params[["max_vax"]],
        t_since_production * vax_production_params[["production_rate"]])
  }
}