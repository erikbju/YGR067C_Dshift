growth_rate <- function(growth_data, platenum, strains, phases, timespan) {
  
  combination <- expand.grid(platenum, strains, phases)
  
  res <- data.frame(Phase = combination[,3], Strain =  combination[,2], Plate = combination[,1], x0 = 0, mu = 0, R2 = 0)
  
  # exp_labels <- c(paste0(combination[,3], " ", combination[,2], ", Plate ", combination[,1]))
  # res <- data.frame(matrix(ncol = 3, nrow = 0))
  # colnames(res) <- c("x0 [g/L]", "μ [1/h]", "Adjusted R2")
  
  for (i in 1:nrow(combination)) {
    my_data <- growth_data %>% subset(plate == combination[i,1] & Strain == combination[i,2] & Phase == combination[i,3] & hours > timespan[1] & hours < timespan[2])
    model <- lm(log(my_data$OD560) ~ my_data$hours)
    res[i,4:6] <- round(c(exp(model$coefficients[1]), model$coefficients[2], summary(model)$adj.r.squared), digits = 3)
  }
  
  # rownames(res) <- exp_labels
  
  return(res)
  
  # model.df <- data.frame(hours = test$hours, OD560 =exp(fitted(model))) %>% unique() %>% mutate(Strain = strainkind) %>% mutate(Phase = phasekind)
}