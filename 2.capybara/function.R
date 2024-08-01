## functions

### https://math.stackexchange.com/questions/453113/how-to-merge-two-gaussians
calculate.params <- function(fitted.params) {
  mu.est <- sum(fitted.params$parameters$mu * fitted.params$parameters$pi)
  var.num <- (fitted.params$parameters$sigma)^2 
  var.est <- sum(var.num * (fitted.params$parameters$pi)^2)
  
  return(list(mu.est, var.est))
}

model.estimation <- function(background) {
  d <- density(background)
  his <- hist(background, breaks = 100)
  
  df <- data.frame(mid=his$mids, cou=his$counts)
  ts_y <- ts(d$y)
  tp <- turnpoints(ts_y)
  guemea <- d$x[tp$peaks]
  guesig <- rep((max(df$mid) - min(df$mid))/4, length(guemea))
  
  guedis <- "norm"
  fitpro <- mix(as.mixdata(df), mixparam(mu=guemea, sigma=guesig))
  
  return(fitpro)
}

calculate.deviance.p.val <- function(deviance.df, centers.smu) {
  cell.types <- ncol(deviance.df) - 1
  sd.est <- (1/cell.types)/(cell.types - 1)
  deviance.df$p.val.single <- pnorm(deviance.df$total.deviance, mean = centers.smu$single.center, sd = sd.est)
  deviance.df$p.val.multi <- pnorm(deviance.df$total.deviance, mean = centers.smu$multi.center, sd = sd.est)
  deviance.df$p.val.unknown <- pnorm(deviance.df$total.deviance, mean = centers.smu$unknown.center, sd = sd.est)
  
  return(deviance.df)
}

calculate.deviance <- function(qp.test.mtx) {
  cell.types <- ncol(qp.test.mtx) - 2
  deviance.from.all <- abs(qp.test.mtx[,c(1:cell.types)] - 1/cell.types)
  deviance.from.all$total.deviance <- rowSums(deviance.from.all)
  
  return(deviance.from.all)
}

get.thresholds <- function(qp.test.mtx) {
  cell.types <- ncol(qp.test.mtx) - 2
  exp.score <- 1/cell.types
  ## single id thresholds
  single.id.th.top <- (1-exp.score) + (exp.score * (cell.types - 1))
  if (cell.types >= 3) {
    single.id.th.bottom <- (1/2-exp.score) * 2 + (exp.score * (cell.types - 2))
  } else {
    print("No Solid Evidence for Multi-ID! Check the QP scores for continuous measure!")
    single.id.th.bottom <- single.id.th.top
  }
  
  ## multi id thresholds
  multi.id.th.top <- (1/2-exp.score) * 2 + (exp.score * (cell.types - 2))
  if (cell.types > 3) {
    multi.id.th.bottom <- (1/3-exp.score) * 3 + (exp.score * (cell.types - 3))
  } else {
    multi.id.th.bottom <- multi.id.th.top
  }
  
  ## Unknown thresholds
  unknown.center <- 0
  
  return(list(single.center = mean(c(single.id.th.top)),
              multi.center = mean(c(multi.id.th.top, multi.id.th.bottom)),
              unknown.center = unknown.center))
}

##
