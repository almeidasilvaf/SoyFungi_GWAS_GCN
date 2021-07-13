

#---Load packages
library(here)
library(ggplot2)
library(tidyverse)
library(ggpubr)

#----Plot number of SNPs per species----

# Load SNP and atlases metadata info
load(here("data", "stress_samples_for_DGE_analysis.RData"))
stress_info <- Reduce(rbind, biotic_samplelist)
stress_info$Pathogen[stress_info$Pathogen == "SBA"] <- "Aglycines"
stress_info$Pathogen[stress_info$Pathogen == "Pgregata"] <- "Cgregata"
snp_info <- read.csv(here("data", "snp_list.csv"), header=TRUE, sep=";")
snp_info$Organism[snp_info$Organism == "Pgregata"] <- "Cgregata"

# Create a stress class ontology for each pathogen/pest
stress_ontology <- data.frame(
    Pathogen = c("Aglycines", "SBA", "Fgraminearum", "Foxysporum",
                 "Fvirguliforme", "Hglycines", "MAMP", "Mphaseolina",
                 "Pgregata", "Ppachyrhizi", "Psojae", "Psojae_glucan_elicitor",
                 "Rreniformis", "Rsolani", "Slitura", "SMV",
                 "Ssclerotiorum", "Evarivestis", "Efabae", "Pincludens",
                 "Agemmatalis", "Xaxonopodis", "Cgregata", "Dphaseolorum",
                 "BPMV", "PMV", "Mincognita", "TRSV", 
                 "Psylvaticum", "Fequiseti"),
    Class = c("insect", "insect", "fungus", "fungus",
              "fungus", "nematode", "MAMP", "fungus", 
              "fungus", "fungus", "oomycete", "oomycete",
              "nematode", "fungus", "insect", "virus",
              "fungus", "insect", "insect", "insect",
              "insect", "bacterium", "fungus", "fungus",
              "virus", "virus", "nematode", "virus",
              "oomycete", "fungus")
)

# Count occurrences for SNP and transcriptome
st1 <- stress_ontology %>%
    inner_join(snp_info, by = c("Pathogen" = "Organism")) %>%
    filter(Class == "fungus")
write_tsv(st1, here("products", "tables", "sup_table1.tsv"))

count_snp <- stress_ontology %>%
    inner_join(snp_info, by = c("Pathogen" = "Organism")) %>%
    select(Pathogen, Class) %>%
    group_by(Pathogen) %>%
    summarise(n = n())

st2 <- stress_ontology %>%
    inner_join(stress_info)
write_tsv(st2, here("products", "tables", "sup_table2.tsv"))

count_transcriptome <- stress_ontology %>%
    inner_join(stress_info) %>%
    select(Pathogen, Class) %>%
    replace("SBA", "Aglycines") %>%
    group_by(Pathogen) %>%
    summarise(n = n())

# Merge data frames by pathogen 
count_snp_transcriptome <- full_join(count_snp, count_transcriptome,
                                     by = c("Pathogen" = "Pathogen")) %>%
    dplyr::rename(SNP = n.x, 
           Transcriptome = n.y) %>%
    replace(is.na(.), 0) %>%
    filter(SNP > 5 | Transcriptome > 5) %>%
    left_join(stress_ontology)

# Plot data
# Fungi
fungi <- count_snp_transcriptome %>%
    pivot_longer(!c(Pathogen, Class)) %>%
    filter(Class == "fungus") %>%
    arrange(Pathogen) %>%
    mutate(Pathogen = fct_relevel(Pathogen)) %>%
    dplyr::rename(Species = Pathogen, Frequency = value) %>%
    ggbarplot(., x="Species", y="Frequency", fill="Species",
              palette=ggsci::pal_d3("category20")(20),
              facet.by="name", orientation="horiz",
              legend="none")

fungi_overlap <- count_snp_transcriptome %>%
    filter(SNP > 0 & Transcriptome > 0) %>%
    pivot_longer(!c(Pathogen, Class)) %>%
    filter(Class == "fungus") %>%
    arrange(Pathogen) %>%
    mutate(Pathogen = fct_relevel(Pathogen)) %>%
    dplyr::rename(Species = Pathogen, Frequency = value) %>%
    ggbarplot(., x="Species", y="Frequency", fill="Species", color="Species",
              palette=ggsci::pal_aaas()(10),
              facet.by="name", orientation="horiz",
              legend="none") +
    theme_bw() +
    ggtitle("SNPs and RNA-seq samples") +
    theme(plot.title = element_text(face="bold", size=12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none")
    
ggsave(fungi_overlap, width=6, height=6,
       filename = here("products", "plots", 
                       "frequency_of_snps_and_transcriptome_samples_overlap.png"))

#----Plot SNP positions in the genome----
library(GenomicRanges)
library(GenomeInfoDb)

# Load SNP ranges
load(here("products", "result_files", "snp_granges.rda"))

# Load gff and remove scaffolds
R.utils::gunzip(here("data", "PLAZA_selected.transcripts.gff.gz"), remove=FALSE)
gff <- gread::read_gff(here("data", "PLAZA_selected.transcripts.gff"))
file.remove(here("data", "PLAZA_selected.transcripts.gff"))
scaffolds <- grep("scaffold", seqlevels(gff), value=TRUE)
gff <- dropSeqlevels(gff, scaffolds, pruning.mode = "tidy")

# Create a different object for each genome location
table(gff$feature)


granges_all <- gread::construct_introns(gff, update=TRUE)
intergenic <- gaps(granges_all)
introns <- granges_all[granges_all$feature == "intron", ]
exons <- gff[gff$feature == "exon", ]
fivep_utr <- gff[gff$feature == "five_prime_UTR", ]
threep_utr <- gff[gff$feature == "three_prime_UTR", ]


# Create a data frame containing the position of each SNP in the genome
snps_intergenic <- as.data.frame(sapply(snp_grangeslist, function(x) 
    length(subsetByOverlaps(x, intergenic))))
snps_intron <- as.data.frame(sapply(snp_grangeslist, function(x) 
    length(subsetByOverlaps(x, introns))))
snps_exon <- as.data.frame(sapply(snp_grangeslist, function(x) 
    length(subsetByOverlaps(x, exons))))
snps_fivep <- as.data.frame(sapply(snp_grangeslist, function(x) 
    length(subsetByOverlaps(x, fivep_utr))))
snps_threep <- as.data.frame(sapply(snp_grangeslist, function(x) 
    length(subsetByOverlaps(x, threep_utr))))
location_of_snps <- cbind(snps_intergenic, 
                          snps_intron, 
                          snps_exon, 
                          snps_fivep,
                          snps_threep)
names(location_of_snps) <- c("Intergenic", "Intron", "Exon", "Fiveprime_UTR", 
                             "Threeprime_UTR")
head(location_of_snps)

# Melt data frame
location_of_snps$Organism <- rownames(snps_exon)
location_of_snps_plotdata <- reshape2::melt(location_of_snps)
location_of_snps_plotdata$variable <- factor(location_of_snps_plotdata$variable, 
                                             levels=c("Intergenic", "Exon", 
                                                      "Intron",
                                                      "Fiveprime_UTR", 
                                                      "Threeprime_UTR"))

snp_positions <- location_of_snps_plotdata %>%
    ggpubr::ggbarplot(data=.,
                      x="variable", y="value", 
                      facet.by = "Organism",
                      fill="Organism", color="Organism", palette = "aaas",
                      orientation="horiz", ncol=5, legend="none",
                      xlab="Genome location", ylab="SNP frequency",
                      label=TRUE, lab.pos="out", lab.hjust = -0.3, lab.vjust = 0.1,
                      title="Location of SNPs in the genome", font.title="bold") +
    ggplot2::expand_limits(y = 100) + 
    ggplot2::scale_y_continuous(expand = c(0,0)) +
    ggplot2::theme(axis.text.x = element_text(size=11)) +
    scale_x_discrete(labels=c("Threeprime_UTR" = "3'-UTR", 
                              "Fiveprime_UTR" = "5'-UTR",
                              "Intron" = "Intron",
                              "Exon" = "Exon",
                              "Intergenic" = "Intergenic")) +
    theme_bw() +
    theme(plot.title = element_text(face="bold", size=12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none",
          axis.text.x = element_text(angle=45, hjust=1))
    

snp_positions


#----Circos plot of SNP positions across chromosomes per stress class----
library(cageminer)
load(here("data", "soybean_genome_ranges.rda"))
gff_gene <- gff[gff$feature == "gene", ]
snp_circos <- plot_snp_circos(genome.ranges, gff_gene, snp_grangeslist)


#----Plot SNP distribution across chromosomes----
library(ggsci)
snp_granges$trait <- snp_granges$Organism
snp_dist <- plot_snp_distribution(snp_granges, "trait") +
    scale_x_discrete(labels = rev(c(
        "Gm20", "Gm19", "Gm18", "Gm17", "Gm16", 
        "Gm15", "Gm14", "Gm13", "Gm12", "Gm11",
        "Gm10", "Gm09", "Gm08", "Gm07", "Gm06", 
        "Gm05", "Gm04", "Gm03", "Gm02", "Gm01" 
    ))) +
    scale_fill_aaas()

# Save plots individually
ggsave(snp_circos, filename = here("products", "plots", "snp_circos.png"),
       width=8, height=8)
ggsave(snp_dist, filename = here("products", "plots", "snp_dist.png"),
       width=8, height=4)
ggsave(snp_positions, filename = here("products", "plots", "snp_positions.png"),
       width=8, height=4)


library(ggpubr)
panel1 <- ggarrange(fungi_overlap, snp_circos, ncol=2, widths = c(2,3))
panel2 <- ggarrange(snp_dist, snp_positions, nrow=2, labels=c("B", "C"))
final_plot <- ggarrange(panel1, panel2, nrow=2, labels=c("A", ""))
ggsave(final_plot, filename = here("products", "plots",
                                   "main_stats_SNPs.pdf"),
       width=9, height=12)










