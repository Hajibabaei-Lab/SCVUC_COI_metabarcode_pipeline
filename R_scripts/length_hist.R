# By Teresita M. Porter, Sept. 6, 2019

library(Biostrings)
library(ggplot2)

args=(commandArgs(TRUE))

# grab seq lengths from FASTA file
lengths<-fasta.seqlengths(args[1], seqtype="DNA")

# create df for ggplot
df<-as.data.frame(lengths)

# Get outlier values
#outlier_values <- boxplot.stats(df$lengths)$out  # outlier values.

# Calculate CDS outliers, values that are lower/higher than 1.5*IQR
# IQR = inter quartile range

# calculate percentiles
percentile25th = quantile(df$lengths, probs=0.25)
# 309
percentile75th = quantile(df$lengths, probs=0.75)
# 321

IQR = percentile75th-percentile25th
# 12

lowerCutoff = percentile25th-(1.5*IQR)
# 291

upperCutoff = percentile75th+(1.5*IQR)
# 339
                        
# plot histogram
h <- ggplot(df, aes(df$lengths)) + 
  geom_histogram(col="black", fill="white", binwidth=1) +
  geom_vline(xintercept = lowerCutoff, col="red4") +
  geom_vline(xintercept = upperCutoff, col="red4") +
  labs(title="cds.fasta", x="Length (bp)", y="Number of CDS (ESVs)") +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))

ggsave(args[2], h)

# combine lower and upper cutoffs
range<-as.data.frame(rbind(lowerCutoff, upperCutoff))

# write to csv
write.table(range, file = args[3], row.names = FALSE, col.names = FALSE)


