
# DATA ACQUISITION

library(here)
library(GenomicRanges)

#----Download .vcf file for genome version Wm82.a2.v1----
options(timeout = 200)
download.file(url="https://soybase.org/snps/soysnp50k_wm82.a2_41317.vcf.gz", 
              destfile = here("data", "snp_positions_a2.v1.vcf.gz"))
system2("ls")
system2("cat", "data/snp_positions_a2.v1.vcf.gz | zcat | grep -v '^##' | cut -f1-3", 
        stdout = "data/snp_positions_a2.v1.txt")
system2("rm", "data/snp_positions_a2.v1.vcf.gz")


#----Download positions of SNPs starting with BARC*----
download.file(url="https://www.soybase.org/SeqMapSearch/searchSeqMap4.php?xsome=all&category=marker&version=Glyma2.0&sid=0.0030980404742696477", 
              destfile = here("data", "BARC_positions.txt"))


#----Download .gff file from PLAZA----
download.file(url="ftp://ftp.psb.ugent.be/pub/plaza/plaza_public_dicots_04/GFF/gma/annotation.selected_transcript.all_features.gma.gff3.gz", 
              destfile = here("data", "PLAZA_selected.transcripts.gff.gz"))

#----Create GRanges object with chromosome lengths----
options(timeout=300)
genome <- Biostrings::readDNAStringSet("ftp://ftp.psb.ugent.be/pub/plaza/plaza_public_dicots_04/Genomes/gma.con.gz")
chr_size <- as.data.frame(seqlengths(genome))[1:20, , drop=FALSE]

genome.ranges = GRanges(
    seqnames = rownames(chr_size),
    strand = "*",
    ranges = IRanges(start = 1, width = chr_size[,1])
)
seqlengths(genome.ranges) <- width(genome.ranges)
save(genome.ranges, file = here("data", "soybean_genome_ranges.rda"), compress="xz")


#----Soybase annotation----
# Go to https://soybase.org/genomeannotation/ and click on: 
# "Download All Genome Annotations for Version Wm82.a2.v1"
soybase <- read.csv("~/Documents/soybase_genome_annotation_v2.0_04-21-2020.txt",
                    header=TRUE, sep="\t", skip=10)
ath_best_hit <- data.frame(
    Gene = rownames(soybase),
    Ath_ortho = soybase[,8]
)
ath_best_hit$Ath_ortho <- gsub("\\.[0-9]\\|.*", "", ath_best_hit$Ath_ortho)
ath_best_hit <- ath_best_hit[!is.na(ath_best_hit$Ath_ortho), ]
save(ath_best_hit, file = here("data", "ath_best_hit.rda"), compress="xz")

