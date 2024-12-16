Hypothesis=extinc.plots[[1]]
Nullmodel=w_null
#{    grey88
  #if (class(Hypothesis$sims)[2] == "SimulateExt") {
    NumExt <- sd <- AccSecExt <- AccSecExt_mean <- Lower <- Upper <- nsim <- NULL
    g <- ggplot(Nullmodel$sims, aes(x = NumExt, y = AccSecExt_mean)) + 
      geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "grey88") + 
      geom_line(aes(color = "blue")) + 
      ylab("Acc. Secondary extinctions") + 
      xlab("Primary extinctions") + theme_bw() #+ ggplot2::theme(axis.title.x = ggplot2::element_blank())
    g <- g + geom_line(data = Hypothesis$sims, aes(y = AccSecExt), 
                        color = "black") + geom_line(data = Hypothesis$sims, 
                                                     aes(y = AccSecExt, color = "black")) + scale_color_manual(name = "Comparison", 
                                                                                                               values = c("black", "blue"), label = c("Observed", "Null hypothesis"))
   
                                                                                                                                           
g + geom_line(data = extinc.plots[[1]]$sims, aes(y=AccSecExt),color = "red", linewidth=1.01)+
geom_line(data = extinc.plots[[2]]$sims, aes(y=AccSecExt),color = "blue", linewidth=1.01)+
geom_line(data = extinc.plots[[3]]$sims, aes(y=AccSecExt),color = "gray", linewidth=1.01)+
geom_line(data = extinc.plots[[4]]$sims, aes(y=AccSecExt),color = "black", linewidth=1.01)+
geom_line(data = extinc.plots[[5]]$sims, aes(y=AccSecExt),color = "cyan", linewidth=1.01)+
geom_line(data = extinc.plots[[6]]$sims, aes(y=AccSecExt),color = "orange", linewidth=1.01)+
geom_line(data = extinc.plots[[7]]$sims, aes(y=AccSecExt),color = "purple", linewidth=1.01)+
geom_line(data = extinc.plots[[8]]$sims, aes(y=AccSecExt),color = "darkgreen", linewidth=1.01)+
geom_line(data = extinc.plots[[9]]$sims, aes(y=AccSecExt),color = "darkgreen", linewidth=1.01)+
geom_line(data = extinc.plots[[10]]$sims, aes(y=AccSecExt),color = "darkgreen", linewidth=1.01)



g + geom_line(data = extinc.plots[[1]]$sims, aes(y=AccSecExt),color = plotcol[1], linewidth=1.01)+
  geom_line(data = extinc.plots[[2]]$sims, aes(y=AccSecExt),color = plotcol[2], linewidth=1.01)+
  geom_line(data = extinc.plots[[3]]$sims, aes(y=AccSecExt),color = plotcol[3], linewidth=1.01)+
  geom_line(data = extinc.plots[[4]]$sims, aes(y=AccSecExt),color = plotcol[4], linewidth=1.01)+
  geom_line(data = extinc.plots[[5]]$sims, aes(y=AccSecExt),color = plotcol[5], linewidth=1.01)+
  geom_line(data = extinc.plots[[6]]$sims, aes(y=AccSecExt),color = plotcol[6], linewidth=1.01)+
  geom_line(data = extinc.plots[[7]]$sims, aes(y=AccSecExt),color = plotcol[7], linewidth=1.01)+
geom_line(data = extinc.plots[[8]]$sims, aes(y=AccSecExt),color = plotcol[8], linewidth=1.01)+
geom_line(data = extinc.plots[[9]]$sims, aes(y=AccSecExt),color = plotcol[9], linewidth=1.01)+
geom_line(data = extinc.plots[[10]]$sims, aes(y=AccSecExt),color = plotcol[10], linewidth=1.01)

  g+ geom_line(data=Nullmodel$sims, aes(x = NumExt, y = AccSecExt_mean), color="red")                                                                                                                                                      
                                                                                                                                                                                                                                                                                                         "Null hypothesis"))
    #h <- ggplot(Nullmodel$sims, aes(x = NumExt, y = nsim/max(nsim))) + 
    #  geom_line() + theme_bw() + labs(y = "Prop", x = "Primary extinctions")
    #I <- patchwork::wrap_plots(g, h, ncol = 1, guides = "keep", 
    #                           heights = c(3, 1)) + patchwork::plot_annotation(tag_levels = "A")
    #I
    #return(I)
 # }
 # if (class(Hypothesis$sims)[2] %in% c("Mostconnected", "ExtinctionOrder")) {
    #NumExt <- sd <- AccSecExt <- AccSecExt_mean <- NULL
    #g <- ggplot(Nullmodel$sims, aes(x = NumExt, y = AccSecExt_mean)) + 
    #  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = scales::muted("red")) + 
    #  geom_line(aes(color = "blue")) + ylab("Acc. Secondary extinctions") + 
    #  xlab("Primary extinctions") + theme_bw() + ggplot2::theme(axis.title.x = ggplot2::element_blank())
    #g <- g + geom_point(data = Hypothesis, aes(y = AccSecExt), 
    #                    color = "black") + geom_line(data = Hypothesis, aes(y = AccSecExt, 
    #                                                                        color = "black")) + scale_color_manual(name = "Comparison", 
                                                                                                                   values = c("black", "blue"), label = c("Observed", 
                                                                                                                                                          "Null hypothesis"))
    #I <- patchwork::wrap_plots(g, h, ncol = 1, guides = "keep", 
      #                         heights = c(3, 1))
    #I
    #return(I)
  #}
  ##else {
  #  message("Hypothesis not of class Mostconnected or ExtinctionOrder")
  #}
#}
#<bytecode: 0x000002960bf42748>
#  <environment: namespace:NetworkExtinction>
#  > 