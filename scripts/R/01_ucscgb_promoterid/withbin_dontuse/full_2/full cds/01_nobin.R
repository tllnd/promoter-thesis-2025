#setwd('~/Library/Mobile Documents/com~apple~CloudDocs/UMassPhD/UMassChan/Tai_Gao/Canavan/UCSC genome browser analyses/2025-08_mAspa figure redo/')
library(ggplot2)
library(ggfx)
library(patchwork)
library(dplyr)

# load in current tracks file from UCSC genome browser, and parse separate tracks
file <- "./2025-08-19_mm9_mAspaP_AllTracks_Full cds.tracks.tsv"
lines <- readLines(file)
# record start and end line of each section/track
starts <- grep("^track\\s+name=", lines)
stops <- c(starts[-1]-1L, length(lines))
# init empty vector to catch the start/stop separated sections
catch <- vector('list', length(starts))
# do the separation
for(i in seq_along(starts)){
  catch[[i]] <- lines[starts[i]:stops[i]]
}
# init empty vector to fill with each section converted to a df
dfs <- vector("list", length(catch))
# get the trackname line, extract just the name from it
for (j in seq_along(catch)) {
  h  <- catch[[j]][1]
  nm <- strsplit(strsplit(h, 'name="', fixed=TRUE)[[1]][2], '"', fixed=TRUE)[[1]][1]
  nm <- chartr(" ", "_", nm)
  
  # ignore first (name) line from each section, split out the text, remove tab separators and brackets
  x <- unlist(strsplit(catch[[j]][-1], "\t", fixed=TRUE), use.names=FALSE)
  x <- trimws(x)
  x <- x[nzchar(x)]
  # as long as there are lines with text, keep going
  if (!length(x)) next
  
  # if line starts with a bracket, its an encode data section, so process it
  if (substr(x[1],1,1)=="[") {
    x <- substring(x, 2, nchar(x)-1)
    # make sure there are 3 entries per data point (start, end, value)
    parts <- strsplit(x,",", fixed = TRUE)
    parts <- Filter(function(v) length(v)==3, parts)
    # make them numeric
    nums <- as.numeric(unlist(parts))
    # reshape and stick in a df, name the df the trackname
    m <- matrix(nums, ncol=3, byrow=TRUE)
    df <- as.data.frame(m)
    colnames(df) <- c("start","end","value")
    } else { #if the line does not start with a bracket, it is something else like a BED block, so capture it in a generic table here
      dat <- paste(catch[[j]][-1], collapse="\n")
      df <- read.table(text=dat, sep="\t", header=FALSE, quote="\"", comment.char="", stringsAsFactors=FALSE)
    }
  # save df in dfs vector 
  dfs[[j]] <- df
  names(dfs)[j] <- nm
}

# Keep only 3-column tracks, this restricts data to start, end, value generally. 
#removes nondata tracks in this case. eg. chr11 and ensGene tracks.
dfs3 <- Filter(function(d) is.data.frame(d) && ncol(d)==3, dfs)

# add name column to each entry in dfs3 so that name of the dataset is present as a column for each row, then bind all dfs into one

out <- bind_rows(lapply(names(dfs3), function(nm){
  d <- dfs3[[nm]]
  d$name <- nm
  d
}))

# clean up column names
# pull out the dfs3 names assigned in the loop, ignoring start and end column names
badnames <- (out[,4])

# trim off the text before the name of each assay, keep everything after
goodcore <- sub(".*([Dd]nase.*|[Hh]istone.*)", "\\1", badnames)
# trim text off the ends, starting from 'std' or 'sig'
# this leaves a name like AssayTissueMarkerAgeStrain
goodcore <- sub("(?i)(Std|Sig).*", "", goodcore, perl=TRUE)
# assign the cleaned names back to the df
out[,4] <- goodcore
# separate the information in the name with pipes
# get the names
nm <- out[,4]

# to all names, put a pipe | before every capital letter
nm <- gsub("([A-Z])", "|\\1", nm, perl = TRUE)
#reassign names with pipes to name column
out[,4] <- nm

# separate by data type...
# pull out dnase
dnase   <- out[grepl('Dnase', out[,4], ignore.case = TRUE),]
# add metadata columns
tmp <- do.call(rbind, strsplit(dnase$name, "|", fixed = TRUE))
dnase[c("assay","tissue","strain","sex","age")] <- tmp[, 2:6]
rm(tmp)
# adults only
dnase_adult <- dnase[grepl('adult', dnase[,4], ignore.case = TRUE),]
# filter for standard tissue set
dnase_adult <- dnase_adult[grepl("Wbrain|Cerebrum|Cerebellum|Heart|Skmuscle|Kidney|Liver", dnase_adult[,4], ignore.case = TRUE),]

# pull out histones...
histone <- out[grepl('Histone', out[,4], ignore.case = TRUE),]
tmp <- do.call(rbind, strsplit(histone$name, "|", fixed = TRUE))
histone[c("assay","tissue","marker","sex","age","strain")] <- tmp[, 2:7]
# filter for adults only
histone_adult <- histone[grepl('adult', histone[,4], ignore.case = TRUE),]
# filter for specific histone markers
histone_adult <- histone_adult[grepl("H3K4ME1|H3K4ME3|H3K27AC", histone_adult[,4], ignore.case = TRUE),]
# remove extra tissues that dont have full data / aren't in dnase
histone_adult <- histone_adult[!grepl("BAT|Spleen", histone_adult[,4], ignore.case = TRUE),]

# order the tissues similarly for each dataset
ord <- c("Wbrain","Cortex","Cerebrum","Cerebellum","Cbellum","Liver","Kidney","Heart","Skmuscle","Spleen","Fat","Gfat","Bat")
dnase_adult$tissue <- factor(dnase_adult$tissue, levels = ord)
histone_adult$tissue <- factor(histone_adult$tissue, levels = ord)
histone_adult$marker <- factor(
  histone_adult$marker,
  levels = c("H3k27ac", "H3k4me1", "H3k4me3")
)

# one common palette (includes DNase)
pal_all <- c(
  DNase   = "grey",
  H3k27ac = "#FFD000",
  H3k4me1 = "#00C4FF",
  H3k4me3 = "#FF2D95"
)

scale_common <- scale_fill_manual(
  name = "Track",
  values = pal_all,
  breaks = names(pal_all),
  limits = names(pal_all),
  drop = FALSE)

# dnase plot
p_dn <- ggplot(dnase_adult, aes(start, value, group = name, fill = "DNase")) +
  geom_area(fill = "grey", position="identity", alpha=0.25)+
  geom_line(color = 'black', linewidth = 0.2, show.legend = FALSE) +
  ylim(0,50)+
  facet_grid(rows = vars(tissue), scales = "free_y") +
  scale_x_reverse() +
  labs(x = "Genomic position (bp)", y = "ENCODE DNase Signal") +
  theme_minimal(base_size = 11) +
  theme(
    strip.text.y = element_text(angle = 0), 
    panel.spacing.y = unit(3, "mm"),
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank())

# histone marker plot
p_chip <- ggplot(histone_adult, aes(start, value, group = marker, fill = marker)) +
# p_chip <- ggplot(histone_adult[grepl("me3", histone_adult[,4]),], aes(start, value, group = tissue, fill = marker)) +
  as_reference(
    geom_area(position="identity", alpha=0.25), id="areas") +
  with_blend(
    geom_area(position="identity", alpha=0.25),
    bg_layer="areas",
    blend_type="multiply")+
  #geom_area(position="identity", alpha=0.25)+
  geom_line(color = 'black', linewidth = 0.2, show.legend = FALSE) +
  ylim(0,5)+
  scale_x_reverse() +
  facet_grid(rows = vars(tissue), scales = "free_y") +
  labs(x = "Genomic position (bp)", y = "ENCODE ChIP Signal") +
  scale_common + 
  theme_minimal(base_size = 11) +
  theme(
    strip.text.y = element_text(angle = 0), 
    panel.spacing.y = unit(3, "mm")
    )

p_combo <- (p_dn+labs(x=NULL)) / p_chip + 
  plot_layout(guides = 'collect') &
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal")&
  guides(color="none")

p_combo

today <- format(Sys.time(), "%Y-%m-%d_%H-%M")
ggsave(
  filename = paste0(today,"_fullcds_nobin_ucscplot.tiff"),
  plot = p_combo,
  device = "tiff",
  width = 8,
  height = 12,
  units = "in",
  dpi = 600,
  compression = "lzw"
)
