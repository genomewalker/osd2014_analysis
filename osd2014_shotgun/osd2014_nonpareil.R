library(tidyverse)
library(Nonpareil)
library(ggpubr)
source("https://gist.githubusercontent.com/genomewalker/8abc47a044f0ac98b7392bbef8afcffa/raw/f4cc661bc3c7db5998b89d01006f710b1eb936d0/geom_flat_violin.R")
source("https://gist.githubusercontent.com/genomewalker/001b93a7415f2652ae0e9d6447a72065/raw/48698cc1de6d8b619747f4fe55d50b2ea523cec7/format_bp.R")
load("osd2014_shotgun/data/osd2014_nonpareil.Rdata")
load("osd2014_shotgun/data/tara_nonpareil.Rdata")

np_parse <- function(X, data = data){
  np <- data[[X]]
  xlim = c(1000, 1e+13)
  ylim = c(1e-06, 1)
  model.x <- exp(seq(log(xlim[1]), log(xlim[2]), length.out = 1000))
  model.y <- predict(np, lr = model.x)
  tibble(label = np@label,x.adj.m = model.x, y.cov.m = model.y, LR = np@LR, model = predict(np), C = np@C, diversity = np@diversity)
}

np_curves_osd <- bind_rows(lapply(1:length(results), np_parse, data = results))

np_plot <- ggplot(np_curves_osd, aes(x.adj, y.cov, group = label)) +
  geom_line(aes(x.adj.m, y.cov.m), alpha = 0.5, color = "#103C54") +
  #geom_line() +
  geom_point(aes(LR,model), shape = 21, fill = "white", color = "#103C54") +
  scale_x_log10(labels = format_bp(), limits = c(1e5, 1e13), breaks = c(1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12, 1e13)) +
  scale_y_continuous(labels = scales::percent)+
  theme_light() +
  theme(legend.position = "none") +
  annotation_logticks(sides = "b") +
  xlab("Sequencing effort") +
  ylab("Estimated average coverage")

cov_plot <- ggplot(np_curves_osd %>% select(C) %>% unique %>% mutate(var = "Coverage"), aes(var, C)) +
  geom_flat_violin(scale = "count", trim = FALSE) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
               geom = "pointrange", position = position_nudge(0.05)) +
  geom_dotplot(binaxis = "y", dotsize = 0.5, stackdir = "down",
               position = position_nudge(-0.025)) +
  theme_light() +
  scale_y_continuous(labels = scales::percent) +
  theme(legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  xlab("") +
  ylab("Estimated average coverage")

div_plot <- ggplot(np_curves_osd %>% select(diversity) %>% unique %>% mutate(var = "Diversity"), aes(var, diversity)) +
  geom_flat_violin(scale = "count", trim = FALSE) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
               geom = "pointrange", position = position_nudge(0.05)) +
  geom_dotplot(binaxis = "y", dotsize = 0.5, stackdir = "down",
               position = position_nudge(-0.025)) +
  theme_light() +
  theme(legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  xlab("") +
  ylab("Nonpareil diversity")

np_comb <- ggpubr::ggarrange(np_plot, cov_plot, div_plot, nrow = 1, ncol = 3, widths = c(0.5, 0.2, 0.2))
osd_np_plot <- np_plot
ggsave(np_comb, filename = "osd2014_shotgun/figures/osd2014_nonpareil_plot.pdf", width = 11.69, height = 4)


# TARA nonpareil ----------------------------------------------------------
load("osd2014_shotgun/data/tara_nonpareil.Rdata")


np_curves_tara <- bind_rows(lapply(1:length(results_np_tara), np_parse, data = results_np_tara))

np_plot <- ggplot(np_curves_tara, aes(x.adj, y.cov, group = label)) +
  geom_line(aes(x.adj.m, y.cov.m), alpha = 0.5, color = "#103C54") +
  #geom_line() +
  geom_point(aes(LR,model), shape = 21, fill = "white", color = "#103C54") +
  scale_x_log10(labels = format_bp(), limits = c(1e5, 1e13), breaks = c(1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12, 1e13)) +
  scale_y_continuous(labels = scales::percent)+
  theme_light() +
  theme(legend.position = "none") +
  annotation_logticks(sides = "b") +
  xlab("Sequencing effort") +
  ylab("Estimated average coverage")

cov_plot <- ggplot(np_curves_tara %>% select(C) %>% unique %>% mutate(var = "Coverage"), aes(var, C)) +
  geom_flat_violin(scale = "count", trim = FALSE) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
               geom = "pointrange", position = position_nudge(0.05)) +
  geom_dotplot(binaxis = "y", dotsize = 0.5, stackdir = "down",
               position = position_nudge(-0.025)) +
  theme_light() +
  scale_y_continuous(labels = scales::percent) +
  theme(legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  xlab("") +
  ylab("Estimated average coverage")

div_plot <- ggplot(np_curves_tara %>% select(diversity) %>% unique %>% mutate(var = "Diversity"), aes(var, diversity)) +
  geom_flat_violin(scale = "count", trim = FALSE) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
               geom = "pointrange", position = position_nudge(0.05)) +
  geom_dotplot(binaxis = "y", dotsize = 0.5, stackdir = "down",
               position = position_nudge(-0.025)) +
  theme_light() +
  theme(legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  xlab("") +
  ylab("Nonpareil diversity")

np_comb <- ggpubr::ggarrange(np_plot, cov_plot, div_plot, nrow = 1, ncol = 3, widths = c(0.5, 0.2, 0.2))
tara_np_plot <- np_plot
ggsave(np_comb, filename = "osd2014_shotgun/figures/tara_nonpareil_plot.pdf", width = 11.69, height = 4)


bind_rows(np_curves_osd %>% mutate(study = "OSD"),
          np_curves_tara %>% mutate(study = "TARA")) %>%
  ggplot(aes(study, diversity)) +
  geom_boxplot() +
  ggpubr::stat_compare_means() +
  theme_light()

ggpubr::ggarrange(osd_np_plot, tara_np_plot, ncol = 1, nrow = 2, align = "hv")



# nonpareil for tara randomised samples -----------------------------------

library(Nonpareil)
library(tidyverse)
f <- list.files(path="~/Downloads/test_nonpareil", pattern = "npo", full.names = FALSE)

sample.names <- sapply(strsplit(f, ".npo"), `[`, 1)

np_curves_osd <- bind_rows(lapply(1:length(results), np_parse, data = results_np_tara))


results_np_tara <- lapply(file.path("~/Downloads/test_nonpareil",f), Nonpareil.curve)


np_plot <- ggplot(np_curves_osd, aes(x.adj, y.cov, group = label)) +
  geom_line(aes(x.adj.m, y.cov.m), alpha = 0.5, color = "#103C54") +
  #geom_line() +
  geom_point(aes(LR,model), shape = 21, fill = "white", color = "#103C54") +
  scale_x_log10(labels = format_bp(), limits = c(1e5, 1e13), breaks = c(1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12, 1e13)) +
  scale_y_continuous(labels = scales::percent)+
  theme_light() +
  theme(legend.position = "none") +
  annotation_logticks(sides = "b") +
  xlab("Sequencing effort") +
  ylab("Estimated average coverage")


names(results_np_tara) <- sample.names

purrr::map_df(results_np_tara, "C") %>% gather() %>%
  ggplot(aes(1,value)) +
  geom_boxplot()

purrr::map_df(results_np_tara, "diversity") %>% gather() %>%
  ggplot(aes("1", value)) +
  geom_boxplot()
