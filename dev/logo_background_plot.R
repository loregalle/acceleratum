library("tidyverse")

load_all()

brst <- burst(muskox, "accelerations_raw", ts_col = "timestamp")

acc <- b2a(brst[1:1000,])
acc <- downsample(acc, ds_factor = 8, FUN = mean)
acc <- as.data.frame(acc)

p <- acc |>
  mutate(timestamp = row_number()) |>
  pivot_longer(!timestamp) |>
  ggplot(aes(timestamp, value, colour = name)) +
  geom_line(lwd = .5) +
  scale_color_brewer(palette = "Set2") +
  xlab("Time") +
  ylab("Acceleration")

p +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",
                                        colour = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        legend.position = "none",
        text = element_text(colour = "white"),
        # axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_line(arrow = grid::arrow(length = unit(2, "mm"),
                                                     angle = 25,
                                                     type = "closed"),
                                 colour = "white"))

p +
  theme_void() +
  theme(panel.background = element_rect(fill = "transparent",
                                        colour = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        legend.position = "none")

ggsave("logo.svg",
       bg = "transparent",
       width = 160, height = 50, units = "mm")
