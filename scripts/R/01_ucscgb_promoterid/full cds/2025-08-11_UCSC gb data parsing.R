#setwd('~/Library/Mobile Documents/com~apple~CloudDocs/UMassPhD/UMassChan/Tai_Gao/Canavan/UCSC genome browser analyses/2025-08_mAspa figure redo/')
library(ggplot2)
library(ggfx)
library(patchwork)

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

# this will make all dfs in dfs into separate dfs in the environment
#list2env(dfs, .GlobalEnv)

# Keep only 3-column tracks, this restricts data to start, end, value generally. 
#removes nondata tracks in this case. eg. chr11 and ensGene tracks.
dfs3 <- Filter(function(d) is.data.frame(d) && ncol(d)==3, dfs)

# Goal: Normalize bin widths (each ENCODE track does not necessarily report values for the same binned sequence regions), so plots will be mis-aligned if not normalized
# Plan: Get common bin = GCD of per-track MODE widths (robust to odd/partial spans)
## get most frequent track width (table(w) outputs a table with frequency counts for different widths), then greatest common denominator for all mode widths in dfs3
mode_width <- function(d){ w <- as.integer(d$end)-as.integer(d$start) 
                           w <- w[w>0]
                           as.integer(names(which.max(table(w)))) 
                           }
gcd <- function(a,b) if (b==0) abs(a) else Recall(b, a %% b)

# define the bin - apply mode_width to each 3 column track in dfs3, get values as an integer vector, 
# run gcd on adjacent mode-widths in the vector, left to right. Final value is the GCD of mode widths 
bin <- Reduce(gcd, vapply(dfs3, mode_width, integer(1)))

# Set normalized bin - build a shared grid (same width and phase)
# get the minimum start coordinate in all dfs3
all_start <- min(vapply(dfs3, function(d) min(as.integer(d$start), na.rm=TRUE), integer(1)))
# get the max end coordinate in all dfs3
all_end   <- max(vapply(dfs3, function(d) max(as.integer(d$end),   na.rm=TRUE), integer(1)))
# set start point: divide all_start by bin, use floor to round down to nearest integer. Multiply by bin to get back adjusted start position
g0 <- floor(all_start/bin)*bin
# define bin edges by going from g0 by bin-width to all_end
edges <- seq(g0, all_end, by=bin)
# set a template df that contains start and end columns, start has all but last edge, end has all but first
out <- data.frame(start=edges[-length(edges)], end=edges[-1])

# Paint each track onto that template as a step function
# NOTE, except for the first and last bin of each track, this logic will not increase bin size
# because bin was set with the gcd of mode_width bin widths
# I did not set bin with all bin widths, but each dataset should only have 2-3 bin widths
# the first bin may be shorter (say 720-800bp), most are length 100 (eg 800-900bp), and the last bin could be weird (eg 900-910bp)
# So all but the start and end bin will definitely be decreased in size unless they are already the appropriate size
# this means that for most bins, I create smaller bins which repeat the data value of the original larger bin they are derived from
# start and end bins *may*, not definitely, be inflated in length, and their value assigned to a larger span than the original data reflected (though it will be maximally increased by bin width, which is small)
for(nm in names(dfs3)){
  # load each df loop by loop
  d <- dfs3[[nm]]
  # pull the start column
  s <- as.integer(d$start)
  # pull the end column
  e <- as.integer(d$end)
  # pull the data column
  v <- as.numeric(d$value)
  # make temp df the shape of the 'out' df
  col <- rep(NA_real_, nrow(out))
  # for each integer in the start column, assign to a normalized bin in 'out':
  for(i in seq_along(s)){ 
    # take the sequence coordinate, subtract the start point and divide by bin width (round down with floor), add one, and test if it is >1
    # g0 is the earliest start coordinate in all data, so this logic captures everything from start point to start point + bin width in bin 1
    k1 <- max(1, floor((s[i]-g0)/bin)+1)
    # take end sequence coordinate, apply same logic (but round up with ceiling). test if result is smaller than the total number of bins
    k2 <- min(nrow(out), ceiling((e[i]-g0)/bin))
    # if start bin is <=  end bin, assign value at i to corresponding bin in the df 'col'
    if(k1<=k2) col[k1:k2] <- v[i] }
  # end of the loop, assign a column with the track name from dfs3, give it the binned values in df 'col'
  out[[nm]] <- col
}
# 'out' has aligned bins (start,end) and one value column per track
# set NAs to zero
out[is.na(out)] <- 0

# clean up column names
# pull out the dfs3 names assigned in the loop, ignoring start and end column names
badnames <- names(out[-(1:2)])
# trim off the text before the name of each assay, keep everything after
goodcore <- sub(".*([Dd]nase.*|[Hh]istone.*)", "\\1", badnames)
# trim text off the ends, starting from 'std' or 'sig'
# this leaves a name like AssayTissueMarkerAgeStrain
goodcore <- sub("(?i)(Std|Sig).*", "", goodcore, perl=TRUE)
# assign the cleaned names back to the df
colnames(out)[-(1:2)] <- goodcore
# separate the information in the name with pipes
# get the names
nm <- names(out)
# flag the start and end 
i <- -(1:2)
# to all but the flagged start and end, put a pipe | before every capital letter
nm[i] <- paste0(substr(nm[i], 1, 1),
                gsub("([A-Z])", "|\\1", substring(nm[i], 2), perl = TRUE))
names(out) <- nm

# function to reshape the data df for plotting, pulling the sample data out of the name as it goes
# arguments include location in pipe-separated names of specific info
make_long <- function(x, tissue_i = 2, age_i = 5, track_i = 3, sep = "|", suffix = "_long") {
  # make sure this is being used on a df like 'out', that has the start and end columns
  stopifnot(all(c("start","end") %in% names(x)))
  # pull out the string text of the df 'x' name (ignoring data), add the suffix to the end of it
  nm  <- deparse(substitute(x))
  nmL <- paste0(nm, suffix)
  # get the metadata from the name - split column names, except for start and end, using the 'sep' argument character
  meta <- strsplit(names(x)[-(1:2)], sep, fixed = TRUE)
  # function to get all metadata tagged to the argument defined in k, for all entries in the meta vector
  pick <- function(k) sapply(meta, function(v) if (length(v) >= k) v[k] else NA_character_)
  # use pick to get the vector for each metadata-argument id
  tissue <- pick(tissue_i) 
  age <- pick(age_i)
  track <- pick(track_i)
  
  # make the df data a matrix without the first two start/end columns, this lets it be flattened
  m <- as.matrix(x[, -(1:2), drop = FALSE])
  # define the position that we will plot with as the midpoint of the bin
  pos <- (x$start + x$end) / 2
  # flatten the matrix and assign colums
  # ncol(m) refers to the number of tracks (columns) in m
  # we are repeating the pos vector n times, where n is the number of columns in m
  # this is then the first column of the new df
  long <- data.frame(pos = rep(pos, ncol(m)),
                     # as.numeric takes m and makes it an integer vector ('single column') which is then assigned a column
                     value = as.numeric(m),
                     # add original column name label as a row to each row in the flattened/long df
                     name = rep(colnames(m), each = nrow(m)),
                     # do the same for tissue, age, track
                     tissue = rep(tissue, each = nrow(m)),
                     age = rep(age, each = nrow(m)),
                     track = rep(track, each = nrow(m)))
  # make long a new df with the name nmL, put it in global env
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

# plottable df for dnase seq data for adults only
make_long(dnase_adult)

# pull out histone data
histone <- out[, c(1, 2, grep("histone", names(out), ignore.case=TRUE)), drop=FALSE]
# filter for adults only
histone_adult <- histone[,c(1,2, grep("adult", names(histone), ignore.case = TRUE)), drop = FALSE]
# filter for specific histone markers
hist_ad_select <- histone_adult[,c(1,2, grep("H3K4ME1|H3K4ME3|H3K27AC", names(histone_adult), ignore.case = TRUE)), drop =FALSE]
# remove extra tissues that dont have full data / aren't in dnase
hist_ad_select <- hist_ad_select[, -grep("BAT|Spleen", names(hist_ad_select), ignore.case = TRUE), drop =FALSE]

# plottable df for histone data for adults only with filters
make_long(hist_ad_select)

# order the tissues similarly for each dataset
ord <- c("Wbrain","Cortex","Cerebrum","Cerebellum","Cbellum","Liver","Kidney","Heart","Skmuscle","Spleen","Fat","Gfat","Bat")
dnase_adult_long$tissue <- factor(dnase_adult_long$tissue, levels = ord)
hist_ad_select_long$tissue <- factor(hist_ad_select_long$tissue, levels = ord)

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
  drop = FALSE)

# dnase plot
p_dn <- ggplot(dnase_adult_long, aes(pos, value, group = name, fill = "DNase")) +
  geom_area(position="identity", alpha=0.25)+
  geom_line(color = 'black', linewidth = 0.2, show.legend = FALSE) +
  facet_grid(rows = vars(tissue), scales = "free_y") +
  scale_x_reverse() +
  scale_common +
  labs(x = "Genomic position (bp)", y = "ENCODE DNase Signal") +
  theme_minimal(base_size = 11) +
  theme(
    strip.text.y = element_text(angle = 0), 
    panel.spacing.y = unit(3, "mm"),
    legend.position = "none")

# histone marker plot
p_chip <- ggplot(hist_ad_select_long, aes(pos, value, group = name, fill = track)) +
  as_reference(
    geom_area(position="identity", alpha=0.25), id="areas") +
  with_blend(
    geom_area(position="identity", alpha=0.25), 
    bg_layer="areas", 
    blend_type="multiply")+
  geom_line(color = 'black', linewidth = 0.2, show.legend = FALSE) +
  scale_x_reverse() +
  scale_common + 
  facet_grid(rows = vars(tissue), scales = "free_y") +
  labs(x = "Genomic position (bp)", y = "ENCODE ChIP Signal") +
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
  filename = paste0(today,"_fullcds_ucscplot.tiff"),
  plot = p_combo,
  device = "tiff",
  width = 8,
  height = 12,
  units = "in",
  dpi = 600,
  compression = "lzw"
)
