#setwd('~/Library/Mobile Documents/com~apple~CloudDocs/UMassPhD/UMassChan/Tai_Gao/Canavan/UCSC genome browser analyses/2025-08_mAspa figure redo/')
library(ggplot2)
library(ggfx)
library(patchwork)

file <- "./2025-08-11_mm9_encode dnase and chip.tracks.tsv"
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

# this will make all dfs in dfs into separate dfs in the environment
#list2env(dfs, .GlobalEnv)

# Keep only 3-column tracks
dfs3 <- Filter(function(d) is.data.frame(d) && ncol(d)==3, dfs)
# Normalize bin widths
# Common bin = GCD of per-track MODE widths (robust to odd/partial spans)
mode_width <- function(d){ w <- as.integer(d$end)-as.integer(d$start) 
                           w <- w[w>0]
                           as.integer(names(which.max(table(w)))) 
                           }
gcd <- function(a,b) if (b==0) abs(a) else Recall(b, a %% b)
bin <- Reduce(gcd, vapply(dfs3, mode_width, integer(1)))

# Build a shared grid (same width and phase)
all_start <- min(vapply(dfs3, function(d) min(as.integer(d$start), na.rm=TRUE), integer(1)))
all_end   <- max(vapply(dfs3, function(d) max(as.integer(d$end),   na.rm=TRUE), integer(1)))
g0 <- floor(all_start/bin)*bin
edges <- seq(g0, all_end, by=bin)
out <- data.frame(start=edges[-length(edges)], end=edges[-1])

# Paint each track onto that grid as a step function
for(nm in names(dfs3)){
  d <- dfs3[[nm]]
  s <- as.integer(d$start)
  e <- as.integer(d$end)
  v <- as.numeric(d$value)
  col <- rep(NA_real_, nrow(out))
  for(i in seq_along(s)){ 
    k1 <- max(1, floor((s[i]-g0)/bin)+1)
    k2 <- min(nrow(out), ceiling((e[i]-g0)/bin))
    if(k1<=k2) col[k1:k2] <- v[i] }
  out[[nm]] <- col
}
# 'out' has aligned bins (start,end) and one value column per track
# set NAs to zero
out[is.na(out)] <- 0

# clean up column names
badnames <- names(out[-(1:2)])
goodcore <- sub(".*([Dd]nase.*|[Hh]istone.*)", "\\1", badnames)
goodcore <- sub("(?i)(Std|Sig).*", "", goodcore, perl=TRUE)
colnames(out)[-(1:2)] <- goodcore
# separate the information in the name with pipes
nm <- names(out); i <- -(1:2)
nm[i] <- paste0(substr(nm[i], 1, 1),
                gsub("([A-Z])", "|\\1", substring(nm[i], 2), perl = TRUE))
names(out) <- nm

# function to reshape the data df for plotting, pulling the sample data out of the name as it goes
make_long <- function(x, tissue_i = 2, age_i = 5, track_i = 3, sep = "|", suffix = "_long") {
  stopifnot(all(c("start","end") %in% names(x)))
  nm  <- deparse(substitute(x)); nmL <- paste0(nm, suffix)
  
  meta <- strsplit(names(x)[-(1:2)], sep, fixed = TRUE)
  pick <- function(k) sapply(meta, function(v) if (length(v) >= k) v[k] else NA_character_)
  tissue <- pick(tissue_i) 
  age <- pick(age_i)
  track <- pick(track_i)
  
  m <- as.matrix(x[, -(1:2), drop = FALSE])
  pos <- (x$start + x$end) / 2
  long <- data.frame(pos = rep(pos, ncol(m)),
                     value = as.numeric(m),
                     name = rep(colnames(m), each = nrow(m)),
                     tissue = rep(tissue, each = nrow(m)),
                     age = rep(age, each = nrow(m)),
                     track = rep(track, each = nrow(m)))
  assign(nmL, long, envir = parent.frame())
  invisible(long)
}

# separate by data type...
# pull out dnase
dnase   <- out[, c(1, 2, grep("dnase",   names(out), ignore.case=TRUE)), drop=FALSE]
make_long(dnase)
# adults only
dnase_adult <- dnase[,c(1,2, grep("adult", names(dnase), ignore.case = TRUE)), drop = FALSE]
# filter for standard tissue set
dnase_adult <- dnase_adult [,c(1,2, grep("Wbrain|Cerebrum|Cerebellum|Heart|Skmuscle|Kidney|Liver", names(dnase_adult), ignore.case = TRUE)), drop = FALSE]


make_long(dnase_adult)

# pull out histone data
histone <- out[, c(1, 2, grep("histone", names(out), ignore.case=TRUE)), drop=FALSE]
# filter for adults only
histone_adult <- histone[,c(1,2, grep("adult", names(histone), ignore.case = TRUE)), drop = FALSE]
# filter for specific histone markers
hist_ad_select <- histone_adult[,c(1,2, grep("H3K4ME1|H3K4ME3|H3K27AC", names(histone_adult), ignore.case = TRUE)), drop =FALSE]
# remove extra tissues that dont have full data / aren't in dnase
hist_ad_select <- hist_ad_select[, -grep("BAT|Spleen", names(hist_ad_select), ignore.case = TRUE), drop =FALSE]

make_long(hist_ad_select)

# order the tissues similarly for each dataset
ord <- c("Wbrain","Cortex","Cerebrum","Cerebellum","Cbellum","Liver","Kidney","Heart","Skmuscle","Spleen","Fat","Gfat","Bat")
dnase_adult_long$tissue <- factor(dnase_adult_long$tissue, levels = ord)
hist_ad_select_long$tissue <- factor(hist_ad_select_long$tissue, levels = ord)

pal <- c("#FFD000", "#00C4FF", "#FF2D95")
lv  <- levels(factor(hist_ad_select_long$track))
pal <- setNames(rep(pal, length.out = length(lv)), lv)

# dnase plot
p_dn <- ggplot(dnase_adult_long, aes(pos, value, group = name, fill = 'grey')) +
  geom_area(position="identity", fill = 'grey', alpha=0.25)+
  geom_line(color = 'black', linewidth = 0.2) +
  facet_grid(rows = vars(tissue), scales = "free_y") +
  scale_x_reverse() +
  labs(x = "Genomic position (bp)", y = "ENCODE DNase Signal") +
  theme_minimal(base_size = 11) +
  theme(strip.text.y = element_text(angle = 0), panel.spacing.y = unit(3, "mm"))

# histone marker plot
p_chip <- ggplot(hist_ad_select_long, aes(pos, value, group = name, fill = track)) +
  as_reference(
    geom_area(position="identity", alpha=0.25), 
    id="areas") +
  with_blend(
    geom_area(position="identity", alpha=0.25), 
    bg_layer="areas", 
    blend_type="multiply")+
  geom_line(color = 'black', linewidth = 0.2) +
  scale_color_manual(values = pal) +
  scale_fill_manual(values   = pal) +
  facet_grid(rows = vars(tissue), scales = "free_y") +
  scale_x_reverse() +
  labs(x = "Genomic position (bp)", y = "ENCODE ChIP Signal") +
  theme_minimal(base_size = 11) +
  theme(strip.text.y = element_text(angle = 0), panel.spacing.y = unit(3, "mm"))

p_dn2 <- p_dn + 
  labs(x=NULL)+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_chip2 <- p_chip +
  theme(legend.position = "bottom",
        legend.direction = "vertical")
  

p_combo <- (p_dn2 / p_chip2) + plot_layout(guides = 'collect')
p_combo
ggsave(
  filename = "2025-08-13_DNase and ChIP orig set.tiff",
  plot = p_combo,
  device = "tiff",
  width = 7,
  height = 10,
  units = "in",
  dpi = 600,
  compression = "lzw"
)
