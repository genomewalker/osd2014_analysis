library(tidyverse)
library(Biostrings)
library("mixtools")

# Learn distribution from data with shorter read fragment

short_length <- data.table::fread(input = "/scratch/antonio/tara_rand_reads/tara_1M_stats.tsv", header = FALSE) %>%
  tbl_df() %>%
  dplyr::rename(label = V1, len = V2)

# Calculate proportion merged/not merged
d <- short_length %>%
  mutate(dist = ifelse(len > 101, TRUE, FALSE)) %>%
  group_by(label) %>%
  count(dist) %>%
  mutate(P = n/sum(n)) %>%
  ungroup() %>%
  group_by(dist) %>%
  summarise(min = min(P), mean = mean(P), median = median(P), max = max(P))


# Use mixture models to fit both distributions
plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}

# We take a random sample of 10M lengths
set.seed(1)
sam <- short_length$len
sam10M <- sam[sample(length(sam), 10000000)]
hist(sam10M)
mixmdl <- normalmixEM(sam10M, k = 2, fast = TRUE, verb = TRUE)

lplot<-data.frame(x = mixmdl$x) %>%
  ggplot() +
  geom_histogram(aes(x, ..density..), binwidth = 1, color = "black", fill = "white", size = 0.1) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
                colour = "firebrick2", lwd = 0.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
                colour = "deepskyblue4", lwd = 0.5) +
  ylab("Density") +
  xlab("Sequence length (bp)") +
  theme_light()


# We have the length distribution, now let's generate the new distribution
setwd("/scratch/antonio/tara_rand_reads")

file_list <- list.files(path = "OSD", pattern = "*fasta", full.names = TRUE)

chop_seqs <- function(X){

  for (iter in 1:100) {

    infile <- X

    fname <- tools::file_path_sans_ext(basename(infile))
    file_out <- file.path("OSD","new",paste(fname, ".", iter, ".new.fasta", sep = ""))

    read_length <- round(mean(short_length$len))

    seqs<-readDNAStringSet(infile)

    # Total number of nucleotides
    nt <- sum(width(seqs))

    # get summary sequence length

    # How many reads of 100nt I need to get closer to nt
    total_reads <- round(nt/(read_length))
    reads_per_seq <- total_reads/length(seqs)

    # Sample random reads
    positions <- sample(1:length(seqs), total_reads, replace=TRUE)


    seqs.new <- seqs[positions]

    # First the larger than 101 bp
    probs_small <- mixmdl$lambda[1]*dnorm(x = min(short_length$len):max(short_length$len), mean = mixmdl$mu[1], sd = mixmdl$sigma[1])
    nlen_small <- sample(min(short_length$len):max(short_length$len), size = round(length(seqs.new)*d$mean[[1]]), replace=T, prob = probs_small)


    probs_large <- mixmdl$lambda[2]*dnorm(x = min(short_length$len):max(short_length$len), mean = mixmdl$mu[2], sd = mixmdl$sigma[2])
    nlen_large <- sample(min(short_length$len):max(short_length$len), size = round(length(seqs.new)*d$mean[[2]]), replace=T, prob = probs_large)

    # bind_rows(data.frame(x = nlen_large) %>% mutate(class = "Large"), data.frame(x = nlen_small) %>% mutate(class = "Small")) %>%
    #   ggplot(aes(x, fill = class)) +
    #   geom_histogram(binwidth = 1, color = "black", size = 0.1, alpha = 0.6) +
    #   ylab("Counts") +
    #   xlab("Sequence length (bp)") +
    #   scale_fill_manual(values = c("deepskyblue4", "firebrick2"), name = element_blank()) +
    #   theme_light()

    r.len <- c(nlen_small, nlen_large)

    seq_df <- tibble(seq_name = names(seqs.new), seq_len = width(seqs.new)) %>% arrange(desc(seq_len)) %>%
      mutate(new_len = sort(r.len, decreasing = TRUE)) %>%
      mutate(diff = ifelse(seq_len <= new_len, TRUE, FALSE),
             diff_n = new_len - seq_len,
             new_len = ifelse(diff == TRUE, new_len - diff_n, new_len)) %>%
      select(-contains("diff")) %>%
      rowwise() %>%
      mutate(seq_ini = sample(1:abs(seq_len-new_len), size = 1)) %>%
      ungroup() %>%
      mutate(seq_end = seq_ini + new_len - 1,
             seq_ini = ifelse(seq_ini == 0, 1, seq_ini)) %>%
      unite(seq_name_export, seq_name, seq_ini, seq_end, sep = "_", remove = FALSE) %>%
      unique()

    seqs.new<-subseq(seqs.new[seq_df$seq_name],start=seq_df$seq_ini,end=seq_df$seq_end)
    names(seqs.new)<-seq_df$seq_name_export
    writeXStringSet(seqs.new,filepath=file_out)
    write_tsv(seq_df,path = paste(file_out,"info",sep="."))
  }
}
l <- pbmcapply::pbmclapply(file_list, chop_seqs, mc.cores = 16)







l1 <- c("OSD/OSD116_2014-06-21_0m_NPL022.me.fasta",
        "OSD/OSD141_2014-06-20_5m_NPL022.me.fasta",
        "OSD/OSD155_2014-06-21_1m_NPL022.me.fasta",
        "OSD/OSD174_2014-06-25_3m_NPL022.me.fasta",
        "OSD/OSD26_2014-06-21_0m_NPL022.me.fasta",
        "OSD/OSD49_2014-06-21_2m_NPL022.me.fasta",
        "OSD/OSD63_2014-06-20_0m_NPL022.me.fasta",
        "OSD/OSD90_2014-06-21_2m_NPL022.me.fasta")






# Use mixture models to fit both distributions
plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}

# We take a random sample of 10M lengths
set.seed(1)
sam <- short_length$len
sam10M <- sam[sample(length(sam), 10000000)]
hist(sam10M)
mixmdl_osd <- normalmixEM(osd_reads$X1, k = 2, fast = TRUE, verb = TRUE)

lplot1<-data.frame(x = mixmdl$x) %>%
  ggplot() +
  geom_histogram(aes(x, ..density..), binwidth = 1, color = "black", fill = "white", size = 0.1) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
                colour = "firebrick2", lwd = 0.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
                colour = "deepskyblue4", lwd = 0.5) +
  ylab("Density") +
  xlab("Sequence length (bp)") +
  theme_light()




