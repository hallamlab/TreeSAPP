library(ggplot2)
library(RColorBrewer)
library(Rmisc)
library(dplyr)
library(optparse)

option_list = list(make_option(c("-i", "--input_table"), type="character", default=NULL,
                               help="Tab-separated value file output by Clade_exclusion_analyzer.py", metavar="character"),
                   make_option(c("-p", "--prefix"), type="character", default="~/Clade_Exclusion",
                               help="Prefix for the output files.", metavar="character"),
                   make_option(c("-r", "--build_params"), type="character", default="data/tree_data/ref_build_parameters.tsv",
                               help="Path to the TreeSAPP reference package build parameters table.", metavar="character"));

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input_table)){
  print_help(opt_parser)
  stop("--input_table must be provided!", call.=FALSE)
}

##
# For debugging
##
# prefix <- "~/Desktop/Trimming"
# input_table <- "manuscript/alignment_trimming_clade_exclusion_performance_compilation.tsv"
# build_params <- "~/Bioinformatics/Hallam_projects/TreeSAPP/data/tree_data/ref_build_parameters.tsv"
# opt <- data.frame(input_table, prefix, build_params, stringsAsFactors = FALSE)

spec_out <- paste(opt$prefix, "Classification_specificity_bars.png", sep='_')
sens_out <- paste(opt$prefix, "Classification_sensitivity.png", sep='_')
pdist_out <- paste(opt$prefix, "Classification_WeightedDistance.png", sep='_')

acc_dat <- read.table(opt$input_table,
                      sep="\t", header=TRUE)
clade_exclusion_header <- c("Trial", "Tool", "Marker", "Rank", "Sequences", "Classified", "Distance", "Proportion")
names(acc_dat) <- clade_exclusion_header
acc_dat <- filter(acc_dat, Rank %in% c("Strain", "Species", "Genus", "Family", "Order", "Class"))
acc_dat$GraftM <- grepl(pattern = "graftm", acc_dat$Tool, ignore.case = T)
acc_dat$Software <- ifelse(acc_dat$GraftM == TRUE, "GraftM", "TreeSAPP")

build_params <- read.table(opt$build_params, sep="\t", header=TRUE)

acc_dat$Distance <- as.character(acc_dat$Distance)
Ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")
Position <- c(1,2,3,4,5,6,7,8)
taxonomic_hierarchy <- data.frame(Ranks, Position)

acc_dat <- merge(acc_dat, taxonomic_hierarchy, by.x = "Rank", by.y = "Ranks")

acc_dat <- acc_dat %>% 
  mutate("Sensitivity" = (Classified/Sequences))


##
# Figure 1: Rank-exclusion specificity in relation to optimal placement distance
##
spec_se <- summarySE(acc_dat, measurevar="Proportion", groupvars=c("Rank", "Trial", "Software", "Distance", "Position")) %>% 
  select(c("Rank", "Trial", "Software", "Distance", "Position", "Proportion", "se")) %>% 
  filter(Distance <= 6) %>% 
  mutate(Rank = reorder(Rank, Position))
filter(spec_se, Distance >= 4) %>% 
  filter(Proportion > 0)
pd <- position_dodge(width = 0.75)

ggplot(spec_se, aes(x=Distance, y=Proportion, fill=Rank)) +
  geom_bar(stat="identity", position=pd, colour="black", width = 0.75) +
  facet_wrap(~Trial + Software) +
  geom_errorbar(aes(ymin=Proportion-se, ymax=Proportion+se),
                width=0.5, position=pd) +
  scale_fill_brewer(palette = "PuOr") +
  scale_y_continuous(breaks=seq(0,100,10)) +
  xlab("Distance from Optimal Rank") +
  ylab("Percentage of Queries") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#bdbdbd", linetype = "dotted"),
        panel.grid.minor = element_blank())

ggsave(filename = spec_out, width = 10, height = 6, dpi = 400)


##
# Figure 2: Rank-exclusion sensitivity in relation to taxonomic rank 
##
ggplot(acc_dat, aes(x=Rank, y=Sensitivity, fill=Rank)) +
  geom_jitter(colour="black", shape=21, size=3, stroke=1) +
  facet_wrap(~Trial) +
  scale_fill_brewer(palette = "PuOr") +
  scale_y_continuous(breaks=seq(0.8,1,0.05), limits = c(0.8,1)) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#bdbdbd", linetype = "dotted"),
        panel.grid.minor = element_blank())

ggsave(filename = sens_out, width = 8, height = 5, dpi = 400)


##
# Figure 3: The taxonomic-hazard-weighted classification scores
## 
acc_dat$Distance <- as.numeric(acc_dat$Distance)
harm_dist_dat <- acc_dat %>% 
  filter(Proportion > 0) %>% 
  merge(build_params, by.x = "Marker", by.y = "name") %>% 
  mutate(Marker = reorder(Marker, ref_sequences)) %>% 
  mutate(Rank = reorder(Rank, Position)) %>%
  mutate(PlaceDist = Distance*Proportion) %>% 
  group_by(Trial, Software, Marker, Rank, ref_sequences) %>% 
  summarise_at(vars(PlaceDist), funs(sum))

ggplot(harm_dist_dat, aes(x=Marker, y=PlaceDist)) +
  geom_point(aes(fill=Rank), colour="black",
             shape=21, size=3, alpha=2/3) +
  scale_fill_brewer(palette = "PuOr") +
  facet_wrap(~Trial) +
  ylab("Cumulative Taxonomic Distance") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#bdbdbd", linetype = "dotted"),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = pdist_out, width = 8, height = 5, dpi = 400)

harm_dist_dat %>%
  group_by(Trial) %>% 
  summarise_at(vars(PlaceDist), funs(mean, median))

##
# Welch Two Sample t-test between the two trials
##
cat("Trials evaluated: ")
unlist(unique(harm_dist_dat$Trial)[1:2])
t.test(filter(harm_dist_dat, Trial == unique(harm_dist_dat$Trial)[1])$PlaceDist,
       filter(harm_dist_dat, Trial == unique(harm_dist_dat$Trial)[2])$PlaceDist)
