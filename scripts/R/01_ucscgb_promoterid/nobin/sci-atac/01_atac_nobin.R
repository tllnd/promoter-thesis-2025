#setwd('~/Library/Mobile Documents/com~apple~CloudDocs/UMassPhD/UMassChan/Tai_Gao/Canavan/UCSC genome browser analyses/2025-08_mAspa figure redo/')
library(ggplot2)
library(ggfx)
library(patchwork)

# load in current tracks file from UCSC genome browser, and parse separate tracks
file <- "./2025-10-19_ucsc_udub sci-atac_mAspaP_cds.tsv"
lines <- readLines(file)

# load in file containing track w/ mAspaP region coordinates
regions <- read.csv("./ucsc-maspap-regions.tsv", sep = "t")

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
goodcore <- sub(".*(Astro.*|Cereb.*|Ex.*|Inh.*|Olig.*|Purk.*|Micro.*|Inter.*)", "\\1", badnames)
# trim text off the ends, starting from 'std' or 'sig'
# this leaves a name like AssayTissueMarkerAgeStrain
goodcore <- sub("(?i)-clusters.*(_\\d+)$", "\\1", goodcore, perl=TRUE)
# assign the cleaned names back to the df
out[,4] <- goodcore


# order the tissues similarly for each dataset
ord <- c("Oligodendrocytes_1","Oligodendrocytes_2","Astrocytes_1","Astrocytes_2","Astrocytes_3",
         "Astrocytes_4","Microglia_3","Ex_neurons_CPN_1","Ex_neurons_SCPN_1","Ex_neurons_SCPN_2",
         "Ex_neurons_CThPN_3","Ex_neurons_CThPN_4","Inhibitory_neurons_1","Inhibitory_neurons_2","Inhibitory_neurons_5",
         "Interneurons_3" ,"Purkinje_cells_1","Cerebellar_granule_cells_2")
out$name <- factor(out$name, levels = ord)

wid <- out$end-out$start
wid2 <- median(wid)


# plot
p_atac <- ggplot(out, aes(start, value, group = name)) +
  geom_col(width = wid2, fill = "grey40", color = NA) +
  # geom_area(position="identity", fill="grey", alpha=0.25)+
  # geom_line(color = 'black', linewidth = 0.2, show.legend = FALSE) +
  facet_grid(rows = vars(name), scales = "free_y") +
  scale_x_reverse() +
  #scale_common +
  labs(x = "Genomic position (bp)", y = "UW sci-ATAC Signal") +
  theme_minimal(base_size = 11) +
  theme(
    strip.text.y = element_text(angle = 0), 
    panel.spacing.y = unit(3, "mm"))
p_atac

today <- format(Sys.time(), "%Y-%m-%d_%H-%M")
ggsave(
  filename = paste0(today,"_sci-atac_nobin_full-cds.tiff"),
  plot = p_atac,
  device = "tiff",
  width = 8,
  height = 12,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

p_atac <- p_atac +
  geom_vline(
    xintercept = c(73137989,73138599,73139400,73140448,73140702,73141400),
    linetype = "dotted", color = "red"
  )
ggsave(
  filename = paste0(today,"_sci-atac_nobin_full-cds_wregs.tiff"),
  plot = p_atac,
  device = "tiff",
  width = 8,
  height = 12,
  units = "in",
  dpi = 600,
  compression = "lzw"
)
