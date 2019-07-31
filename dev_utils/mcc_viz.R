library(ggplot2)
library(RColorBrewer)
library(Rmisc)
library(dplyr)
library(optparse)
library(stringr)

option_list = list(make_option(c("-i", "--input_table"), type="character", default=NULL,
                               help="Tab-separated value file output by MCC_calculator.py", metavar="character"),
                   make_option(c("-p", "--prefix"), type="character", default="~/MCC_Rplot",
                               help="Prefix for the output file (no extension).", metavar="character"));

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

ggplot(mcc_dat, aes(x=Tax.dist, y=MCC, fill=Tool)) +
  geom_line() +
  geom_point(colour="black", shape=21, size=3, stroke=1) +
  scale_fill_brewer(palette = "PuOr") +
  scale_y_continuous(limits=c(0.3,1),
                     breaks=seq(0.3,1,0.2)) +
  scale_x_continuous(breaks=seq(0,7,1)) +
  xlab("Distance from Optimal Rank") +
  ylab("Matthews correlation coefficient") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#bdbdbd", linetype = "dotted"),
        panel.grid.minor = element_blank())

ggsave(filename = spec_out <- paste(opt$prefix, ".eps", sep=''), width = 8, height = 6, dpi = 400)
