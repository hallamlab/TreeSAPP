library(ggplot2)
library(RColorBrewer)
library(Rmisc)
library(dplyr)
library(optparse)
library(stringr)

option_list = list(make_option(c("-i", "--input_table"), type="character", default=NULL,
                               help="Tab-separated value file output by MCC_calculator.py", metavar="character"),
                   make_option(c("-d", "--output_dir"), type="character", default="./",
                               help="Output directory to write output files.", metavar="character"));

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input_table)){
  print_help(opt_parser)
  stop("--input_table must be provided!", call.=FALSE)
}


mcc_dat <- read.table(opt$input_table,
                      sep="\t", header=TRUE)
# header <- c("Tool", "Tax.dist", "MCC", "True.Pos", "True.Neg"," "False.Pos" "False.Neg")
# names(mcc_dat) <- header

mcc_plt <- mcc_dat %>%
  ggplot(aes(x=Tax.dist, y=MCC, fill=Tool)) +
  geom_line() +
  geom_point(colour="black", shape=21, size=5, stroke=1, alpha=0.7) +
  scale_fill_brewer(palette = "Paired") +
  scale_y_continuous(limits=c(0.3,1),
                     breaks=seq(0.3,1,0.2)) +
  scale_x_continuous(breaks=seq(0,7,1)) +
  xlab("Distance from Optimal Rank") +
  ylab("Matthews correlation coefficient") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#bdbdbd", linetype = "dotted"),
        panel.grid.minor = element_blank())

pr_curve <- mcc_dat %>% 
  mutate("Precision" = True.Pos/(True.Pos+False.Pos)) %>% 
  mutate("Recall" = True.Pos/(True.Pos+False.Neg)) %>% 
  ggplot(aes(x=Recall, y=Precision, fill=Tool)) +
  geom_line() +
  geom_point(colour="black", shape=21, size=3, stroke=1) +
  scale_fill_brewer(palette = "Paired") +
  scale_y_continuous(limits=c(0,1)) +
  scale_x_continuous(limits=c(0,1)) +
  xlab("Recall") +
  ylab("Precision") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#bdbdbd", linetype = "dotted"),
        panel.grid.minor = element_blank())

# Save all the files
ggsave(plot= mcc_plt, filename = paste(opt$output_dir, "MCC_Rplot.png", sep='/'),
       width = 8, height = 6, dpi = 400)
ggsave(plot= mcc_plt, filename = paste(opt$output_dir, "MCC_Rplot.svg", sep='/'),
       width = 8, height = 6)
ggsave(plot= pr_curve, filename = paste(opt$output_dir, "Prec-Recall_curve.png", sep='/'),
       width = 8, height = 6, dpi = 400)
ggsave(plot= pr_curve, filename = paste(opt$output_dir, "Prec-Recall_curve.pdf", sep='/'),
       width = 8, height = 6)

