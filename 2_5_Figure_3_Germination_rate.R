library(ggplot2)
library(patchwork)

load("output/rda/plot_alnus.rda")

alnus <- plot +
  ggtitle("Alnus glutinosa") +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "italic", size = 12))

load("output/rda/plot_betula_pub.rda")

betula <- plot +
  ggtitle("Betula pubescens") +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "italic", size = 12))

load("output/rda/plot_pinus.rda")

pinus <- plot +
  ggtitle("Pinus sylvestris") +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "italic", size = 12))


final_plot <- (alnus / betula / pinus) +
  plot_layout(guides = "collect", axis_titles = "collect_x") &
  theme(legend.position = "bottom")

final_plot

ggsave("output/imgs/Fig_3_cardinal_temperatures.png", final_plot, width = 9.5, height = 8.5, units = "in")

# ggsave("output/imgs/production/Fig_3_cardinal_temperatures.pdf", final_plot, width = 9.5, height = 8.5, units = "in")
# 
# ggsave("output/imgs/production/Fig_3_cardinal_temperatures.tif", final_plot, width = 9.5, height = 8.5, units = "in", dpi = 600, compression = "lzw")

