

heatmap <-function(fit, id, ggdata_){
  group1 <-list()
  for (i in 1:length(id)){
    if (abs(fit$U[i,1]) >.01){
      group1[[length(group1)+1]] <-1} else {group1[[length(group1)+1]] <-0}} 
  
  group2 <-list()
  for (i in 1:length(id)){
    if (abs(fit$U[i,2]) > .01){
      group2[[length(group2)+1]] <-1} else {group2[[length(group2)+1]] <-0}} 
  
  group3 <-list()
  for (i in 1:length(id)){
    if (abs(fit$U[i,3]) >.01){
      group3[[length(group3)+1]] <-1} else {group3[[length(group3)+1]] <-0}}
  
  layerdata <- data.frame(ID = unlist(id), layer1= unlist(group1), layer2= unlist(group2), layer3= unlist(group3))
  gglayerdat1 <-merge(ggdata_, layerdata, by = "ID")
  gglayerdat1$layer1 = as.factor(gglayerdat1$layer1)
  gglayerdat1$layer2 = as.factor(gglayerdat1$layer2)
  gglayerdat1$layer3 = as.factor(gglayerdat1$layer3)
  colnames(gglayerdat1)[2:4] = c("PC1", "PC2", "PC3")
  
  x <- ggplot(subset(gglayerdat1, layer1 == "1"), aes(x=PC1, y=PC2)) +
    geom_hex(bins = 20)+ theme(aspect.ratio=1, axis.title.x = element_blank(), axis.ticks = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank())+
    scale_fill_gradientn(colours=c("khaki","red"),name = "Frequency", limit = range(c(1,33))) + ylim(-3.3,4.8) +xlim(-5,5)+labs(title = "Humic-like 1")
  y <- ggplot(subset(gglayerdat1, layer2 == "1"), aes(x=PC1, y=PC2)) +
    geom_hex(bins = 20)+ theme(aspect.ratio=1, axis.title.x = element_blank(), axis.ticks = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position =  "none" )+
    scale_fill_gradientn(colours=c("khaki","red"),name = "Frequency", limit = range(c(1,33)))+ ylim(-3.3,4.8) +xlim(-5,5)+labs(title = "Protein-like")
  z <- ggplot(subset(gglayerdat1, layer3 == "1"), aes(x=PC1, y=PC2)) +
    geom_hex(bins = 20)+ theme(aspect.ratio=1, axis.title.x = element_blank(), axis.ticks = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position =  "none")+
    scale_fill_gradientn(colours=c("khaki","red"),name = "Frequency", limit = range(c(1,33)))+ ylim(-3.3,4.8) +xlim(-5,5)+labs(title = "Humic-like 2")
  legend <- get_legend(x)
  l <- get_legend(x)
  x <- x + theme(legend.position = "none")
  print(plot_grid(x, y, z, l, nrow = 1, rel_widths = c(1, 1, 1, .5)))
  
}


pca_heatmap <- function(ggdata_, loads){
  ggplot(ggdata_, aes(x=Comp.1, y=Comp.2)) + 
    geom_hex(alpha = .7) +
    scale_fill_gradientn(colours=c("darkgrey","darkgrey"), name = "Frequency") +                     
    annotate("segment", x = 0, y = 0, xend = 4 * loads$Comp.1[1], yend = 4 * loads$Comp.2[1], colour = '#2C5B94', size = 1.5, arrow = arrow(length = unit(0.15, "inches"))) +
    annotate("text", x = 4 * loads$Comp.1[1], y = 4 * loads$Comp.2[1], label = 'O/C', colour = "black", hjust = -0.1, vjust = -0.5) +
    annotate("segment", x = 0, y = 0, xend = 4 * loads$Comp.1[2], yend = 4 * loads$Comp.2[2], colour = "#2C5B94", size = 1.5, arrow = arrow(length = unit(0.15, "inches")))  +
    annotate("text", x = 4 * loads$Comp.1[2], y = 4 * loads$Comp.2[2], label = 'H/C', colour = "black", hjust = .8, vjust = -2.7) +
    annotate("segment", x = 0, y = 0, xend = 4 * loads$Comp.1[3], yend = 4 * loads$Comp.2[3], colour = "#2C5B94", size = 1.5, arrow = arrow(length = unit(0.15, "inches")))  +
    annotate("text", x = 4 * loads$Comp.1[3], y = 4 * loads$Comp.2[3], label = 'N/C', colour = "black", hjust = -.5, vjust = 1.1) +
    annotate("segment", x = 0, y = 0, xend = 4 * loads$Comp.1[4], yend = 4 * loads$Comp.2[4], colour = "#2C5B94", size = 1.5, arrow = arrow(length = unit(0.15, "inches")))  +
    annotate("text", x = 4 * loads$Comp.1[4], y = 4 * loads$Comp.2[4], label = 'KMD', colour = "black", hjust = -.2, vjust = -.7) +
    annotate("segment", x = 0, y = 0, xend = 4 * loads$Comp.1[6], yend = 4 * loads$Comp.2[6], colour = "#2C5B94", size = 1.5, arrow = arrow(length = unit(0.15, "inches")))  +
    annotate("text", x = 4 * loads$Comp.1[6], y = 4 * loads$Comp.2[6], label = 'Neutral Mass', colour ="black", hjust = .6, vjust = 1.3) +
    annotate("segment", x = 0, y = 0, xend = 4 * loads$Comp.1[7], yend = 4 * loads$Comp.2[7], colour = "#2C5B94", size = 1.5, arrow = arrow(length = unit(0.15, "inches"))) +
    annotate("text", x = 4 * loads$Comp.1[7], y = 4 * loads$Comp.2[7], label = 'DBE', colour = "black", hjust = 1.1, vjust = -0.5)+
    theme(legend.position = "none", axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), aspect.ratio=1)
}