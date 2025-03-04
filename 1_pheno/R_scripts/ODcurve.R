
ODcurve <- function(growth_data, plate, k = 6) {
  
  mydata <-growth_data[growth_data$plate == plate, ] %>% group_by(hours, Strain, Phase) %>% 
    summarise(medOD = median(OD560)) %>% ungroup() %>% group_by(Strain, Phase) %>% arrange(Strain, Phase) %>% 
    mutate(mm5 = rollapply(medOD, k, median, fill = NA, align="center"), mad5 = rollapply(medOD, k, mad, fill = NA,align="center")) %>% 
    filter(!(Phase == "Pre-shift" & hours > 11.5))
  
  p <- ggplot(mydata, aes(x = hours, y = mm5, color = Strain, linetype = Phase, fill = Strain, group = interaction(Strain, Phase))) +
    # stat_summary(
    #   fun = mean,
    #   geom='line',
    #   size = 0.7) +
    # stat_summary(
    #   fun.data=mean_cl_boot,
    #   geom='ribbon',
    #   alpha=0.2) +
    geom_line() +
    geom_ribbon(aes(ymin = mm5-mad5, ymax = mm5+mad5, fill = Strain), color = NA, alpha = 0.2) +
    geom_vline(xintercept = 11.5, linetype = "dotted", colour = "black") +
    ylim(0, 0.4) +
    labs(x = "Time (h)", y = paste0("Rolling median of OD560 (k=", k, ")")) + 
    ggtitle(paste0("Plate: ", plate)) +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(size = 20),
      axis.line = element_line(colour = "black"),
      # legend.title = element_blank(),
      legend.text  = element_text(size = 10),
      legend.key = element_rect(fill = "white"),
      axis.title.x = element_text(size = 12),
      axis.text.x = element_text(size = 10),
      axis.title.y = element_text(size = 12),
      axis.text.y = element_text(size = 10)
    ) 
  return(p)
}