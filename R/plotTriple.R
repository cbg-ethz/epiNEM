color_list = unique(ggpdat$colorz)
labels_list = unique(ggpdat$logic)

ggplot(data=ggpdat, aes(index, log_liklihood, colour=colorz)) + 
  geom_point() +
  scale_colour_manual(values=setNames(ggpdat$colorz, ggpdat$colorz), labels=setNames(ggpdat$logic, ggpdat$logic)) +
  labs(title= paste(parents2[1], " and ", parents2[2], sep="")) +
  scale_size_manual(values = c(2,3,3)) +
  guides(size = F) +
  coord_cartesian(xlim = c(-1, 164), ylim = c(min(dat)-50, max(dat)+200)) +
  scale_y_continuous(#limits = c(-5,0),
    name = "Log Likelihood") +
  scale_x_discrete(breaks = c(1, index[dim(ggpdat)[1]]),
                   labels= c(1, dim(ggpdat)[1]), name = "Modulators")  +
  #annotate(geom = "rect", xmin = -1, xmax = ggpdat$index[dim(ggpdat)[1]], ymin =  dat[length(dat)-4]-30, ymax = 0, alpha = .15 ) +
  geom_text(aes(label= annotate), size=6, angle = 45, hjust = 0) +
  theme(legend.position = c(.77,.75), legend.key.size = unit(2, "cm"),
        legend.text = element_text(size = 15), legend.title = element_text(size=16),
        plot.title = element_text(size = rel(2)))
  
  
  ggplot(data=ggpdat, aes(index, log_liklihood, colour=colorz)) + 
    geom_point() +
    scale_colour_manual(values=setNames(ggpdat$colorz, ggpdat$colorz)) +
    labs(title= paste(parents2[1], " and ", parents2[2], sep="")) +
    scale_size_manual(values = c(2,3,3)) +
    guides(size = F) +
    coord_cartesian(xlim = c(-1, 164), ylim = c(min(dat)-50, max(dat)+200)) +
    scale_y_continuous(#limits = c(-5,0),
      name = "Log Likelihood") +
    scale_x_discrete(breaks = c(1, index[dim(ggpdat)[1]]),
                     labels= c(1, dim(ggpdat)[1]), name = "Modulators") 