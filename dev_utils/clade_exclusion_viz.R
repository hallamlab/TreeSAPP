library(ggplot2)
library(RColorBrewer)
library(Rmisc)
library(dplyr)
library(tidyr)
library(optparse)
library(stringr)

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
# prefix <- "~/Desktop/treesapp_evaluate_test"
# input_table <- "~/Bioinformatics/Hallam_projects/TreeSAPP_manuscript/clade_exclusion_performance_compilation.tsv"
# build_params <- "~/Bioinformatics/Hallam_projects/TreeSAPP/treesapp/data/ref_build_parameters.tsv"
# opt <- data.frame(input_table, prefix, build_params, stringsAsFactors = FALSE)

##
# Set figure names here
##
spec_out <- paste(opt$prefix, "Classification_specificity_bars.eps", sep='_')
sens_out <- paste(opt$prefix, "Classification_sensitivity.eps", sep='_')
pdist_out <- paste(opt$prefix, "Classification_MeanDistance.eps", sep='_')
f1_out <- paste(opt$prefix, "F1-score.eps", sep='_')
o_u_bars <- paste(opt$prefix, "Prediction_bars.eps", sep='_')

##
# Function definitions
##
f1_score <- function(p, r) {
  numer <- p*r
  denom <- p+r
  return(2*(numer/denom))
}
## end

##
# Loading data, cleaning, merging, mutating, etc.
##
acc_dat <- read.table(opt$input_table,
                      sep="\t", header=TRUE) %>% 
  mutate(Proportion = (Correct*100)/Queries)
acc_dat <- filter(acc_dat, Rank %in% c("Strain", "Species", "Genus", "Family", "Order", "Class"))
acc_dat$Software <- str_replace_all(acc_dat$Tool, c("treesapp" = "TreeSAPP", "graftm" = "GraftM", "diamond" = "DIAMOND"))

build_params <- read.table(opt$build_params, sep="\t", header=TRUE)

acc_dat$TaxDist <- as.character(acc_dat$TaxDist)
Ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")
Position <- c(1,2,3,4,5,6,7,8)
taxonomic_hierarchy <- data.frame(Ranks, Position)

acc_dat <- merge(acc_dat, taxonomic_hierarchy, by.x = "Rank", by.y = "Ranks")
## end

##
# Figure 1: Rank-exclusion specificity in relation to optimal placement distance
##
spec_se <- summarySE(acc_dat, measurevar="Proportion", groupvars=c("Rank", "Software", "TaxDist", "Position")) %>% 
  select(c("Rank", "Software", "TaxDist", "Position", "Proportion", "se")) %>% 
  filter(TaxDist <= 6) %>% 
  mutate(Rank = reorder(Rank, Position))
filter(spec_se, TaxDist >= 4) %>% 
  filter(Proportion > 0)
pd <- position_dodge(width = 0.75)

ggplot(spec_se, aes(x=TaxDist, y=Proportion, fill=Rank)) +
  geom_bar(stat="identity", position=pd, colour="black", width = 0.75) +
  facet_wrap(~Software) +
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
## end


##
# Figure 2: The average taxonomic distance across marker genes
## 
acc_dat$TaxDist <- as.numeric(acc_dat$TaxDist)
harm_dist_dat <- acc_dat %>% 
  filter(Proportion > 0) %>% 
  merge(build_params, by.x = "RefPkg", by.y = "name") %>% 
  mutate(RefPkg = reorder(RefPkg, ref_sequences)) %>% 
  mutate(Rank = reorder(Rank, Position)) %>%
  mutate(PlaceDist = TaxDist*Proportion) %>% 
  group_by(Software, RefPkg, Rank, ref_sequences) %>% 
  summarise_at(vars(PlaceDist), funs(sum)) %>% 
  mutate(AvgDist = PlaceDist/100)

ggplot(harm_dist_dat, aes(x=RefPkg, y=AvgDist)) +
  geom_point(aes(fill=Rank), colour="black",
             shape=21, size=3, alpha=2/3) +
  geom_rug(sides="l") +
  scale_fill_brewer(palette = "PuOr") +
  facet_wrap(~Software) +
  ylab("Average Taxonomic Distance") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#bdbdbd", linetype = "dotted"),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = pdist_out, width = 8, height = 5, dpi = 400)
## end


##
# Classification performance data frame for the following two figures for Recall and F1-score
##
summary_dat <- acc_dat %>% 
  filter(TaxDist == 2) %>% 
  group_by(Software, RefPkg, Rank, Position) %>% 
  summarize(Queries = sum(Queries),
            Cumulative = sum(Cumulative),
            Over = sum(Over),
            Under = sum(Under)) %>% 
  mutate("Precision" = (Cumulative/(Cumulative+Over))) %>% 
  mutate("Recall" = (Cumulative/Queries)) %>%
  ungroup() %>% 
  mutate(Rank = reorder(Rank, Position)) 

##
# Figure 3: Rank-exclusion sensitivity in relation to taxonomic rank 
##
ggplot(summary_dat, aes(x=RefPkg, y=Recall, fill=Rank)) +
  geom_point(colour="black", shape=21, size=3, stroke=1) +
  facet_wrap(~Software) +
  scale_fill_brewer(palette = "PuOr") +
  scale_y_continuous(breaks=seq(0,1.05,0.1), limits = c(0,1)) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#bdbdbd", linetype = "dotted"),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = sens_out, width = 8, height = 5, dpi = 400)
## end


##
# Figure 4: Plot of F1-score across taxonomic distance
##
f1_dat <- summary_dat %>% 
  group_by(Software, Rank, Position) %>% 
  summarize(Queries = sum(Queries),
            Cumulative = sum(Cumulative),
            Over = sum(Over),
            Under = sum(Under)) %>% 
  ungroup() %>% 
  mutate("Precision" = (Cumulative/(Cumulative+Over))) %>% 
  mutate("Recall" = (Cumulative/Queries)) %>% 
  mutate("F1_score" = f1_score(Precision, Recall)) %>% 
  mutate(Rank = reorder(Rank, Position))

ggplot(f1_dat, aes(x=Rank, y=F1_score, group=Software, colour=Software)) +
  geom_point(colour="black", size=3) +
  geom_line(size=2) +
  scale_colour_brewer(palette = "Set2") +
  scale_y_continuous(breaks=seq(0,1,0.1), 
                     limits = c(0,1)) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#bdbdbd", linetype = "dotted"),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = f1_out, width = 7, height = 7, dpi = 400)
## end


##
# Welch Two Sample t-test between the two software
##
harm_dist_dat %>%
  group_by(Software) %>% 
  summarise_at(vars(AvgDist), funs(mean, median))

cat("Softwares evaluated: ")
unlist(unique(harm_dist_dat$Software)[2:3])
t.test(filter(harm_dist_dat, Software == unique(harm_dist_dat$Software)[2])$AvgDist,
       filter(harm_dist_dat, Software == unique(harm_dist_dat$Software)[3])$AvgDist)
## end


##
# Figure 5: Characterizing the proportion of over- and under-predictions
##
over_under_dat <- acc_dat %>% 
  filter(TaxDist == 2) %>% 
  group_by(Software, RefPkg) %>% 
  summarize(Queries = sum(Queries),
            Cumulative = sum(Cumulative),
            Over = sum(Over),
            Under = sum(Under)) %>% 
  ungroup() %>% 
  mutate("P_over" = Over/Queries) %>% 
  mutate("P_under" = (Under/Queries)*-1) %>% 
  gather("P_over", "P_under",
         key="Predict_type", value = "Predict_count")

ggplot(over_under_dat, aes(x=RefPkg, y=Predict_count, fill=Software, colour=Predict_type)) +
  geom_bar(stat = "identity",
           position = "dodge",
           size = 1) +
  xlab("Reference Package") +
  ylab("Proportion of Over- and Under-predictions") +
  scale_fill_brewer(palette = "Set2") +
  scale_colour_manual(values = c("#252525", "#d9d9d9")) +
  guides(colour = FALSE) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#bdbdbd", linetype = "dotted"),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = o_u_bars, width = 8, height = 5, dpi = 400)
## end
