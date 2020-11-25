#!/usr/bin/env Rscript

##
# Boilerplate
##
r = getOption("repos")
r["CRAN"] = "https://cloud.r-project.org/"
options(repos = r)

# Ensure the required packages are installed
deps <- c("optparse", "logging", "ggplot2", "RColorBrewer", "dplyr", "tidyr", "purrr", "stringr")
new.packages <- deps[!(deps %in% installed.packages()[,"Package"])]
if(length(new.packages))
  install.packages(new.packages)

# Load the R packages
library(optparse)
library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(dplyr, warn.conflicts = F)
library(stringr)
library(purrr)

# Set up logging
log_file <- "treesapp_summariser.log.txt"
if (file.exists(log_file))
  #Delete file if it exists
  file.remove(log_file)

library(logging, quietly = TRUE)
logging::basicConfig()
addHandler(writeToFile, logger="ts", file=log_file)
setLevel("INFO")

##

##
# Load functions
##
get_args <- function() {
  option_list = list(make_option(c("-t", "--treesapp_table"), type="character", default=NULL,
                                 help="Path to the classification table made by TreeSAPP.", metavar="character"),
                     make_option(c("-o", "--output_dir"), type="character", default="./",
                                 help="Path to a directory where the figures and tables should be saved. [ DEFAULT = './' ]", metavar="character"))
  
  opt_parser = OptionParser(option_list=option_list);
  opt = parse_args(opt_parser);

  return(opt)
}

# List the package dependency versions
list_dep_versions <- function(dependencies) {
  for (d in deps) {
    logging::loginfo(paste(d, "version:", packageVersion(d)), logger="ts")
  }
}

# A function for generating n number of hex colour codes
MoreColours <- function(n_o, palette="BrBG") {
  ramp_func <- colorRampPalette(brewer.pal(9, palette))
  ramp_func(n_o)
}

# Function for determining the rank a lineage was resolved to
resolved_rank <- function(lineage_str) {
  taxa <- str_split(lineage_str, pattern="; ")
  len <- length(unlist(taxa))
  if (len == 0) {
    return("Unclassified")
  }
  rank <- as.character(with(rank_depth_map, Rank[Depth %in% len])) %>% 
    unlist()
  return(rank)
}

# Load and format the TreeSAPP classification table, marker_contig_map.tsv
load_classification_table <- function(ts_tbl_path) {
  dat <- read.table(file=ts_tbl_path,
                    header=TRUE,
                    sep="\t", quote="", comment.char = "") %>% 
    mutate(Length=End_pos-Start_pos) %>% 
    mutate(Taxonomy = as.character(Taxonomy)) %>% 
    mutate(Resolved = as.character(lapply(Taxonomy, resolved_rank))) %>% 
    separate(col = Taxonomy, into = c('r', 'k', 'p', 'c', 'o', 'f', 'g', 's'), sep = "; ", fill = "right", remove = T) %>% 
    select(Sample, Marker, Length, EvoDist, Resolved, k, p, c, o, f, g)
  return(dat)
}

# Function for plotting evolutionary distance
evodist_density <- function(ts_tbl_df, output_dir, min_length=30) {
  evodist_dat <- ts_tbl_df %>%
    filter(Length >= min_length) %>% 
    select(Marker, EvoDist)
  
  dist_density <- ggplot(evodist_dat, aes(x=EvoDist)) +
    geom_density(colour="black", fill="grey", alpha=0.5) +
    xlim(0,3) +
    ylab("Density") +
    xlab("Evolutionary Distance") +
    facet_wrap(~Marker) +
    treesapp_theme +
    theme(axis.text.x = element_text(angle=0, hjust=0.5))
  ggsave(plot=dist_density,
         filename=paste0(output_dir, "/EvoDist_Density.png"),
         height = 6, width = 8)
}

# Function for making a bubble plot of sample_name by taxonomy
taxa_bubbles <- function(ts_tbl_df, output_dir, max_evodist=2, rank=c) {
  ranks <- list("d"="Domain", "p"="Phylum", "c"="Class", "o"="Order", "f"="Family", "g"="Genus", "s"="Species")
  rank_var <- enquo(rank)
  rank_name = ranks[rlang::as_name(rank_var)]
  
  # Prep the data frame by tallying the different groups
  census_dat <- ts_tbl_df %>%
    filter(!is.na(!! rank_var)) %>% 
    group_by(Sample, Marker, !! rank_var) %>% 
    dplyr::count() %>% 
    ungroup()
  # Generate the plot
  bub_plt <- ggplot(census_dat, aes(x=Sample, y=!! rank_var, fill=!! rank_var)) +
    geom_point(aes(size=n), pch=21, colour="black") +
    facet_wrap(~Marker) +
    scale_fill_viridis_d() +
    guides(fill=F,
           size=guide_legend(title="Relative abundance")) +
    ylab(rank_name) +
    treesapp_theme
  ggsave(plot=bub_plt,
         filename=paste0(output_dir, "/Taxon_bubble_plot.png"),
         height=8,
         width=8)
}

# Function for plotting the taxonomic resolution of the classifications
taxa_res_bar <- function(ts_tbl_df, output_dir) {
  # Make the input data frame
  known_df <- ts_tbl_df %>% 
    group_by(Resolved) %>% 
    dplyr::count() %>% 
    ungroup() %>% 
    merge(rank_depth_map, by.x="Resolved", by.y="Rank") %>% 
    mutate(Proportion = n/sum(n)) %>% 
    mutate(Resolved = reorder(Resolved, Depth))
  scalar <- sum(known_df$n)
  # Generate the plot
  res_bar <- known_df %>% 
    ggplot(aes(fill=Resolved, y=n, x=1)) +
    geom_bar(stat = "identity", colour="black") +
    ylab("Assigned Sequences") +
    rank_palette +
    scale_y_continuous(limits = c(0,scalar),
                       sec.axis = sec_axis(~./scalar, name = "Proportion")) +
    coord_flip() +
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text())
  ggsave(plot=res_bar,
         file=paste0(output_dir, "/classification_resolution.png"),
         height = 3, width = 8)
}

# Function for plotting abundance of different taxa across samples
taxa_bars <- function(ts_tbl_df, output_dir, rank=c) {
  ranks <- list("d"="Domain", "p"="Phylum", "c"="Class", "o"="Order", "f"="Family", "g"="Genus", "s"="Species")
  rank_var <- enquo(rank)
  rank_name = ranks[rlang::as_name(rank_var)]

  # Prep the data frame by tallying the different groups
  census_dat <- ts_tbl_df %>%
    filter(!is.na(!! rank_var)) %>% 
    group_by(Sample, Marker, !! rank_var) %>% 
    dplyr::count() %>% 
    ungroup()
  # Generate the plot
  taxa_hist <- ggplot(census_dat, aes(x=Marker, y=n, fill=!! rank_var)) +
    geom_bar(stat="identity") +
    scale_fill_viridis_d() +
    guides(fill=guide_legend(title=rank_name)) +
    xlab("Protein family") +
    ylab("Count") +
    facet_wrap(~Sample) +
    treesapp_theme +
    theme(plot.margin = unit(c(1, 2, 1, 2), "cm"))
  # Save the output
  ggsave(plot=taxa_hist,
         filename=paste0(output_dir, "/Taxon_abundances.png"),
         height=6,
         width=10)
}

##

##
# Load global variables
##
# A ggplot theme used in all plots
treesapp_theme <- theme(panel.background = element_blank(),
                        panel.grid.major = element_line(colour = "#bdbdbd", linetype = "dotted"),
                        panel.grid.minor = element_blank(),
                        axis.text.x = element_text(angle = 45, 
                                                   hjust = 1, vjust = 1))
# A mapping table between the taxonomic rank and position in a lineage
Rank <- c("Root", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
Depth <- c(1, 2, 3, 4, 5, 6, 7, 8)
rank_depth_map <- data.frame(Rank, Depth,
                             stringsAsFactors = FALSE)


rank_palette <- scale_fill_manual(name="Resolved", drop=TRUE,
                                  values=setNames(MoreColours(length(rank_depth_map$Rank), "PuOr"), rank_depth_map$Rank))

##


##
# Main workflow
##
# Set up the options parser
opt <- get_args()


if (is.null(opt$treesapp_table)){
  print_help(opt_parser)
  stop("--treesapp_table must be provided!", call.=FALSE)
}
# List the dependency versions
list_dep_versions(deps)
# Load the classification table
ts_dat <- load_classification_table(opt$treesapp_table)

# Summarize taxonomic classification resolution
taxa_res_bar(ts_dat, opt$output_dir)

# Make a bubble plot of taxonomy by sample name
taxa_bubbles(ts_dat, opt$output_dir)

# Make a faceted bar plot of taxonomy by sample name
taxa_bars(ts_dat, opt$output_dir, f)

# Make a density plot of evolutionary distance faceted by RefPkg
evodist_density(ts_dat, opt$output_dir)

## End
