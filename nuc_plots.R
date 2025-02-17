library(ggplot2)
library(rtracklayer)
args = commandArgs(trailingOnly=TRUE)
inFile = args[1]
print(inFile)

test <- read.csv(inFile, sep="", header = F)
chr_l <- 250000000 # 102000000

create_tfbs <- function(test, len) {
  test_seq <- rep(0, len)
  for (i in 1:nrow(test)) {
    s <- min(as.numeric(test[i,2]), as.numeric(test[i,3]))
    e <- max(as.numeric(test[i,2]), as.numeric(test[i,3]))
    test_seq[s:e] <- test_seq[s:e] + rep(1, ((e + 1) - s))
  }
  return(test_seq)
}

tfbs_l <- create_tfbs(test, chr_l)

plot_centered <- function(bed, bed_name, window_len) {
 # test2 <- read.csv(bed, sep ="", header=F)
  test2 <- readGFF(bed)
  test2 <- test2[which(test2$class == "W"),]
  get_windows <- function(regions, window_len) {
    starts <- regions$start
    ends <- regions$end
    mids <- round(abs(ends - starts)/2)
    new_starts <- starts + mids - window_len
    new_ends <- starts + mids + window_len
    positions <- data.frame(start = new_starts,
                            end = new_ends)
    return(positions)
  }
  
  windows <- get_windows(test2, window_len)
  windows <- windows[complete.cases(windows),]
  res <- sapply(1:nrow(windows), function(X){
    start <- windows[X,1]
    end <- windows[X,2]
    tfbs_l[start:end] 
  })
  
  line_plot <- rowSums(res)
  write.csv(line_plot, "nuc_vector.csv")
  
  p <- ggplot() +
    geom_smooth(aes(c(-window_len:(window_len)),line_plot),se=F, size=2) +
    xlab("Nuc Dyad") +
    ylab(bed_name) +
    theme_classic()
  ggsave(filename="Rplot2.pdf", plot=p)
  dev.off()
  
}

plot_centered("../../../human_nucs/all_nucs/GSM907784/03_nucleR/GSM907784.60.gff", paste0("TTS Density"), 500)
