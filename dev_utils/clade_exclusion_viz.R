deps <- c("ggplot2", "RColorBrewer", "Rmisc", "dplyr", "tidyr", "optparse", "stringr", "Cairo")

new.packages <- deps[!(deps %in% installed.packages()[,"Package"])]
if(length(new.packages))
  install.packages(new.packages)

library(ggplot2)
library(RColorBrewer)
library(Cairo)
library(tidyr)
library(optparse)
library(stringr)
library(Rmisc)
library(dplyr)

setHook(packageEvent("grDevices", "onLoad"),
        function(...) grDevices::X11.options(type='cairo'))
options(device='x11')
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
# input_table <- "~/Bioinformatics/Hallam_projects/TreeSAPP_manuscript/Performance_analyses/clade_exclusion_performance_compilation.tsv"
# build_params <- "~/Bioinformatics/Hallam_projects/TreeSAPP/treesapp/data/ref_build_parameters.tsv"
# opt <- data.frame(input_table, prefix, build_params, stringsAsFactors = FALSE)

##
# Set figure names here
##
img_format <- "eps"
spec_out <- paste(paste(opt$prefix, "Classification_specificity_bars",  sep='_'),
                  img_format, sep='.')
sens_out <- paste(paste(opt$prefix, "Classification_sensitivity", sep='_'),
                  img_format, sep='.')
pdist_out <- paste(paste(opt$prefix, "Classification_MeanDistance", sep='_'),
                   img_format, sep='.')
f1_out <- paste(paste(opt$prefix, "F1-score", sep='_'),
                img_format, sep='.')
o_u_bars <- paste(paste(opt$prefix, "Prediction_bars", sep='_'),
                  img_format, sep='.')

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
taxonomic_hierarchy <- data.frame(Ranks = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"),
                                  Depth = c(1,2,3,4,5,6,7,8))

refpkg_plt_dat <- data.frame(Cycle = c("Carbon", "Carbon", "Carbon", "Carbon",
                                       "Nitrogen", "Nitrogen", "Nitrogen", "Nitrogen", "Nitrogen", "Nitrogen",
                                       "Sulphur"),
                             RefPkg = c("McrA", "McrB", "McrG", "p_amoA",
                                        "napA", "nirK", "nirS", "nifD", "NorB", "NxrB",
                                        "DsrAB"),
                             Position = as.numeric(seq(1, 11)))

acc_dat <- merge(acc_dat, taxonomic_hierarchy, by.x = "Rank", by.y = "Ranks") %>% 
  merge(refpkg_plt_dat, by="RefPkg")

## end

##
# Figure 1: Rank-exclusion specificity in relation to optimal placement distance
##
spec_se <- summarySE(acc_dat, measurevar="Proportion", groupvars=c("Rank", "Software", "TaxDist", "Depth")) %>% 
  select(c("Rank", "Software", "TaxDist", "Depth", "Proportion", "se")) %>% 
  filter(TaxDist <= 6) %>% 
  mutate(Rank = reorder(Rank, Depth))
filter(spec_se, TaxDist >= 4) %>% 
  filter(Proportion > 0)
pd <- position_dodge(width = 0.75)

spec_plot <- ggplot(spec_se, aes(x=TaxDist, y=Proportion, fill=Rank)) +
  geom_bar(stat="identity", position=pd, colour="black", width = 0.75, alpha = 2/3) +
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

ggsave(plot = spec_plot, filename = spec_out,
       width = 10, height = 6, dpi = 400, device = cairo_ps)
ggsave(plot = spec_plot, filename = gsub(img_format, "png", spec_out),
       width = 10, height = 6, dpi = 400)
## end


##
# Figure 2: The average taxonomic distance across marker genes
## 
acc_dat$TaxDist <- as.numeric(acc_dat$TaxDist)
harm_dist_dat <- acc_dat %>% 
  filter(Proportion > 0) %>% 
  merge(build_params, by.x = "RefPkg", by.y = "name") %>% 
  mutate(RefPkg = reorder(RefPkg, Position)) %>% 
  mutate(Rank = reorder(Rank, Depth)) %>%
  mutate(PlaceDist = TaxDist*Proportion) %>% 
  group_by(Software, RefPkg, Rank, ref_sequences) %>% 
  summarise_at(vars(PlaceDist), funs(sum)) %>% 
  mutate(AvgDist = PlaceDist/100) %>% 
  ungroup()

avg_dist_plot <- harm_dist_dat %>% 
  ggplot(aes(x=RefPkg, y=AvgDist)) +
  geom_point(aes(fill=Rank), colour="black",
             shape=21, size=3, alpha=2/3) +
  geom_rug(sides="l") +
  scale_fill_brewer(palette = "PuOr") +
  facet_wrap(~ Software) +
  xlab("Reference Package") +
  ylab("Average Taxonomic Distance") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#bdbdbd", linetype = "dotted"),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(plot = avg_dist_plot, filename = pdist_out,
       width = 8, height = 5, dpi = 400, device = cairo_ps)
ggsave(plot = avg_dist_plot, filename = gsub(img_format, "png", pdist_out),
       width = 8, height = 5, dpi = 400)
## end


##
# Classification performance data frame for the following two figures for Recall and F1-score
##
summary_dat <- acc_dat %>% 
  filter(TaxDist == 2) %>% 
  group_by(Software, RefPkg, Rank, Depth, Position) %>% 
  summarize(Queries = sum(Queries),
            Cumulative = sum(Cumulative),
            Over = sum(Over),
            Under = sum(Under)) %>% 
  mutate("Precision" = (Cumulative/(Cumulative+Over))) %>% 
  mutate("Recall" = (Cumulative/Queries)) %>%
  ungroup() %>% 
  mutate(RefPkg = reorder(RefPkg, Position)) %>% 
  mutate(Rank = reorder(Rank, Depth)) 

##
# Figure 3: Rank-exclusion sensitivity in relation to taxonomic rank 
##
sens_plot <- summary_dat %>% 
  ggplot(aes(x=RefPkg, y=Recall, fill=Rank)) +
  geom_point(colour="black", shape=21, size=3, stroke=1, alpha=2/3) +
  facet_wrap(~Software) +
  scale_fill_brewer(palette = "PuOr") +
  scale_y_continuous(breaks=seq(0,1.05,0.1), limits = c(0,1)) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#bdbdbd", linetype = "dotted"),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(plot = sens_plot, filename = sens_out,
       width = 8, height = 5, dpi = 400, device = cairo_ps)
ggsave(plot = sens_plot, filename = gsub(img_format, "png", sens_out),
       width = 8, height = 5, dpi = 400)
## end


##
# Figure 4: Plot of F1-score across taxonomic distance
##
f1_dat <- summary_dat %>% 
  group_by(Software, Rank, Depth) %>% 
  summarize(Queries = sum(Queries),
            Cumulative = sum(Cumulative),
            Over = sum(Over),
            Under = sum(Under)) %>% 
  ungroup() %>% 
  mutate("Precision" = (Cumulative/(Cumulative+Over))) %>% 
  mutate("Recall" = (Cumulative/Queries)) %>% 
  mutate("F1_score" = f1_score(Precision, Recall)) %>% 
  mutate(Rank = reorder(Rank, Depth))

fone_plot <- f1_dat %>% 
  ggplot(aes(x=Rank, y=F1_score, group=Software, colour=Software)) +
  geom_point(colour="black", size=3) +
  geom_line(size=2) +
  scale_colour_brewer(palette = "Set2") +
  scale_y_continuous(breaks=seq(0,1,0.1), 
                     limits = c(0,1)) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#bdbdbd", linetype = "dotted"),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(plot = fone_plot, filename = f1_out,
       width = 7, height = 7, dpi = 400, device = cairo_ps)
ggsave(plot = fone_plot, filename = gsub(img_format, "png", f1_out),
       width = 7, height = 7, dpi = 400)
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
  group_by(Software, RefPkg, Position) %>% 
  summarize(Queries = sum(Queries),
            Cumulative = sum(Cumulative),
            Over = sum(Over),
            Under = sum(Under)) %>% 
  ungroup() %>% 
  mutate(RefPkg = reorder(RefPkg, Position)) %>% 
  mutate("P_over" = Over/Queries) %>% 
  mutate("P_under" = (Under/Queries)*-1) %>% 
  gather("P_over", "P_under",
         key="Predict_type", value = "Predict_count")

o_u_plot <- over_under_dat %>% 
  ggplot(aes(x=RefPkg, y=Predict_count, fill=Software, colour=Predict_type)) +
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
ggsave(plot = o_u_plot, filename = o_u_bars,
       width = 8, height = 5, dpi = 400, device = cairo_ps)
ggsave(plot = o_u_plot, filename = gsub(img_format, "png", o_u_bars),
       width = 8, height = 5, dpi = 400)
## end
