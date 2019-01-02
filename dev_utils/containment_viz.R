library(ggplot2)
library(Rmisc)
library(dplyr)
library(RColorBrewer)
library(optparse)

option_list = list(make_option(c("-i", "--input_table"), type="character", default=NULL,
                               help="Tab-separated value file output by Clade_exclusion_analyzer.py", metavar="character"),
                   make_option(c("-p", "--prefix"), type="character", default="out.txt",
                               help="Prefix for the output files.", metavar="character"));

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input_table)){
  print_help(opt_parser)
  stop("--input_table must be provided!", call.=FALSE)
}

contain_out <- paste(opt$prefix, "Containment_bars.png", sep='_')

contain_dat <- read.table(opt$input_table,
                          sep="\t", header=TRUE)
acc_header <- c("Trial", "Marker", "Software", "Rank", "Contained", "Sensitivity")
names(contain_dat) <- acc_header

contain_se <- summarySE(contain_dat, measurevar="Contained", groupvars=c("Trial", "Rank", "Software")) %>% 
  mutate(Rank = reorder(Rank, Contained))

ggplot(contain_se, aes(x=Rank, y=Contained, fill=Rank)) +
  geom_bar(stat="identity") +
  facet_wrap(~Trial+Software) +
  geom_errorbar(aes(ymin=Contained-se, ymax=Contained+se),
                width=0.2) +
  scale_fill_brewer(palette = "PuOr") +
  scale_y_continuous(breaks=seq(0,100,10)) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#bdbdbd", linetype = "dotted"),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x = element_blank())

ggsave(filename=contain_out, width = 8, height = 5)
