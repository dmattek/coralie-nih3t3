require(ggplot2)
require(ggExtra)
require(gridExtra)

rhg_cols <- c(
  "#771C19",
  "#AA3929",
  "#E25033",
  "#F27314",
  "#F8A31B",
  "#E2C59F",
  "#B6C5CC",
  "#8E9CA3",
  "#556670",
  "#000000"
)

md_cols <- c(
  "#FFFFFF",
  "#F8A31B",
  "#F27314",
  "#E25033",
  "#AA3929",
  "#FFFFCC",
  "#C2E699",
  "#78C679",
  "#238443"
)

#####
## Custom functions

myCheckDigits <- function(x) {
  grepl('^[-]?[0-9]+[.]?[0-9]*$' , x)
}

myCheckLogical <- function(x) {
  grepl('^TRUE$|^FALSE$' , x)
}

myConvertStringListToTypes <- function(in.l) {
  # convert strings with digits to numeric
  # uses logical indexing: http://stackoverflow.com/questions/42207235/replace-list-elements-by-name-with-another-list
  loc.l = myCheckDigits(in.l)
  in.l[loc.l] = lapply(in.l[loc.l], as.numeric)
  
  # convert strings with TRUE/FALSE to logical
  loc.l = myCheckLogical(in.l)
  in.l[loc.l] = lapply(in.l[loc.l], as.logical)
  
  return(in.l)
}

# f-n to read experimental description
# returns data table with entire experimental data
myExpRead = function(inFname, inCleanRowCol = TRUE, inCleanMissing = TRUE, inStartRow = 1, inSheetName = 1, inRowIndex = NULL) {
  # read the file
  loc.dt.exp = as.data.table(read.xlsx(
    file = inFname,
    sheetName = inSheetName,
    rowIndex = inRowIndex,
    startRow = inStartRow
  ))
  
  if(inCleanRowCol) {
    # sometimes an NA column appears at the end; remove
    loc.dt.exp = loc.dt.exp[, names(loc.dt.exp)[!(names(loc.dt.exp) %like% 'NA')], with = FALSE]
    
    # sometimes an NA row appears at the end; remove
    loc.dt.exp = loc.dt.exp[loc.dt.exp[,!Reduce(`&`, lapply(.SD, is.na))]]
  }
  
  if(inCleanMissing) {
    # replace missing values with ''
    for (i in seq_along(loc.dt.exp))
      set(loc.dt.exp,
          i = which(is.na(loc.dt.exp[[i]])),
          j = i,
          value = '')
  }
  
  return(loc.dt.exp)
}

## f-n to read single Nuclei.csv file
# removes a number of columns; make sure they're indeed unnecessary
myFread = function(fileIn) {
  # Read the first two rows
  dt.head = fread(fileIn, nrows = 2, header = FALSE)
  
  # make a joint single-row header from two rows
  s.head = paste(dt.head[1], dt.head[2], sep = '_')
  
  # read the rest of the output (except first two rows)
  dt.nuc = fread(fileIn, skip = 2)
  
  # set column names
  setnames(dt.nuc, s.head)
  
  # remove duplicated columns
  dt.nuc = dt.nuc[, s.head[!duplicated(s.head)], with = FALSE]
  
  # remove unnecesary columns
  s.cols = c(
    'Image_ImageNumber',
    'Image_Metadata_C',
    'Image_Metadata_ChannelName',
    'Image_Metadata_ColorFormat',
    'Image_Metadata_FileLocation',
    'Image_Metadata_Frame',
    'Image_Metadata_Plate',
    'Image_Metadata_SizeC',
    'Image_Metadata_SizeT',
    'Image_Metadata_SizeX',
    'Image_Metadata_SizeY',
    'Image_Metadata_SizeZ',
    'Image_Metadata_Well',
    'Image_Metadata_Z'
  )
  dt.nuc[, (s.cols) := NULL]
}


# Returns dt with trajectories that last from the first till last frame
# and include at most in.max.break interruptions
myTrajExtr = function(in.dt,
                      in.max.break = 1,
                      in.met.series = 'Metadata_Series',
                      in.met.t = 'Metadata_T',
                      in.met.tracklabel = 'TrackObjects_Label',
                      in.aggr.cols = NULL) {
  loc.dt = copy(in.dt)
  
  # TrackObjects assigns the same label for different cells from IdentifyPrimaryObjects
  # The following aggregation makes sure there's a unique TrackObjects_Label
  # for every Site at every time point: it takes the mean intensity of duplicated labels.
  # Make sure it makes sense!!!
  # Roughly 10% of objects affected
  
  if (!is.null(in.aggr.cols)) {
    loc.dt = loc.dt[, lapply(.SD, function(x) mean(x, na.rm = TRUE)),
                    by = c(in.met.series, in.met.t, in.met.tracklabel), .SDcols = in.aggr.cols]
  }
  
  # cells from different sites have the same TrackObjects_Label
  # make it unique accross the experiment and pad numbers with zeros
  loc.dt[, TrackObjects_Label_uni := paste(sprintf("%02d", get(in.met.series)), sprintf("%04d", get(in.met.tracklabel)), sep = "_")]
  
  ####
  ## Allow for single-breaks in the middle of the trajectory
  ## Single-breaks usually result from a single frame with lost focus
  ## no interpolation
  
  # build vector with unique timepoints for the entire experiment
  loc.t.range = unique(loc.dt[, c(in.met.t), with = FALSE])
  
  # Select cells with number of timepoints equal to
  # the range nrow(loc.t.range)  AND
  # with first and last frame at the beginning and end of the movie, respectively.
  # If no aggregation columns (in.aggr.cols) are provided,
  # tracks with forks will be omitted here because the number of timepoints exceeds nrow(loc.t.range)
  loc.dt.tmp1 = loc.dt[, .(
    Ntpt = .N,
    T.start = first(get(in.met.t)),
    T.end = last(get(in.met.t))
  ),
  by = TrackObjects_Label_uni][Ntpt <=  nrow(loc.t.range) &
                                 T.start == min(loc.t.range) &
                                 T.end == max(loc.t.range)]
  
  # With cells selected above,
  # select tcourses with at most 1 break.
  ## Single-breaks usually result from a single frame with lost focus
  ## no interpolation
  # Calculation based on differences between consecutive Metadata_T.
  # Metadata_T increases by 1, thus a single-frame break would result in a difference of 2.
  
  # select cells recorded at the beginning and end of the experiment
  loc.dt.tmp2 = loc.dt[TrackObjects_Label_uni %in% loc.dt.tmp1$TrackObjects_Label_uni]
  
  # set order
  setkeyv(loc.dt.tmp2, c(in.met.series, in.met.tracklabel, in.met.t))
  
  # calculate difference in consecutive time
  loc.dt.tmp2[, Metadata_T.diff := c(-1, diff(get(in.met.t))), by = TrackObjects_Label_uni]
  
  # identify cells with at least one break longer than 1 frame
  # column Nframes stores the number of instances with break longer than 1 frame
  loc.dt.tmp2 = loc.dt.tmp2[Metadata_T.diff > 1 + in.max.break, .(Nframes = .N), by = TrackObjects_Label_uni]
  
  # Selected trajectories with frames at 1st and last time points AND with at most 1-frame break
  loc.out = loc.dt[TrackObjects_Label_uni %in% setdiff(loc.dt.tmp1$TrackObjects_Label_uni,
                                                       loc.dt.tmp2$TrackObjects_Label_uni)]
  return(loc.out)
}



# Extract full trajectories
# No aggregation
myTrajExtr2 = function(in.dt,
                       in.max.break = 1, 
                       in.met.series = 'Metadata_Series', 
                       in.met.t = 'Metadata_T', 
                       in.met.tracklabel = 'TrackObjects_Label') {
  loc.dt = copy(in.dt)
  
  # build vector with unique timepoints for the entire experiment
  loc.t.range = unique(loc.dt[, c(in.met.t), with = FALSE])
  
  # add unique cell id
  loc.dt[, TrackObjects_Label_uni := paste(sprintf("%02d", get(in.met.series)), sprintf("%04d", get(in.met.tracklabel)), sep = "_")]
  
  # Select cells with number of timepoints equal to
  # the range nrow(loc.t.range)  AND 
  # with first and last frame at the beginning and end of the movie, respectively.
  loc.dt.tmp1 = loc.dt[, .(Ntpt = .N, 
                           T.start = first(Metadata_T), 
                           T.end = last(Metadata_T)), 
                       by = TrackObjects_Label_uni][Ntpt <=  nrow(loc.t.range) & 
                                                      T.start == min(loc.t.range) & 
                                                      T.end == max(loc.t.range)]
  
  
  # With cells selected above,
  # select tcourses with at most 1 break.
  ## Single-breaks usually result from a single frame with lost focus
  ## no interpolation
  # Calculation based on differences between consecutive Metadata_T.
  # Metadata_T increases by 1, thus a single-frame break would result in a difference of 2.
  
  loc.dt.tmp2 = loc.dt[TrackObjects_Label_uni %in% loc.dt.tmp1$TrackObjects_Label_uni]
  loc.dt.tmp2[, Metadata_T.diff := c(NA, diff(get(in.met.t))), by = TrackObjects_Label_uni]
  loc.dt.tmp2 = loc.dt.tmp2[Metadata_T.diff <= in.max.break + 1]
  
  # Selected trajectories with at most 1-frame break
  loc.out = loc.dt[TrackObjects_Label_uni %in% loc.dt.tmp2$TrackObjects_Label_uni]
  return(loc.out)
}


# Returns original dt with RealTime column added
# Real time is based on acquisition frequency
# Input parameters:
# in.dt - data.table
# in.met.t - string with the name of a column in in.dt with frame number, e.g. Metadata_T
# in.acq.freq - acquisition frequency in minutes (integer)
myAddRealTime = function(in.dt, in.met.t, in.acq.freq) {
  
  loc.dt = copy(in.dt)
  
  loc.dt.t.trans = data.table(meta.tmp = unique(loc.dt[[in.met.t]]), 
                              RealTime = seq(min(loc.dt[[in.met.t]])*in.acq.freq, max(loc.dt[[in.met.t]])*in.acq.freq, in.acq.freq))
  
  setnames(loc.dt.t.trans, 'meta.tmp', in.met.t)
  
  loc.dt = merge(loc.dt, loc.dt.t.trans, by = in.met.t)
  
  return(loc.dt)
}

## Custom plotting

myGgplotTraj = function(dt.arg,
                        x.arg,
                        y.arg,
                        group.arg,
                        facet.arg,
                        facet.ncol.arg = 2,
                        line.col.arg = NULL,
                        xlab.arg = "Time",
                        ylab.arg = "Fl. int.",
                        plotlab.arg = "",
                        dt.stim.arg = NULL,
                        tfreq.arg = 1,
                        maxrt.arg = 60,
                        xaxisbreaks.arg = 10,
                        ylim.arg = c(0,1),
                        stim.bar.height.arg = 0.1,
                        stim.bar.width.arg = 0.5) {
  p.tmp = ggplot(dt.arg,
         aes_string(x = x.arg,
                    y = y.arg,
                    group = group.arg))
  
  if (is.null(line.col.arg))
    p.tmp = p.tmp + geom_line(alpha = 0.25, size = 0.25)
  else
    p.tmp = p.tmp + geom_line(aes_string(colour = line.col.arg), alpha = 0.5, size = 0.5)
  
  p.tmp = p.tmp + 
    stat_summary(
      aes_string(y = y.arg, group = 1),
      fun.y = mean,
      colour = 'blue',
      linetype = 'solid',
      size = 1,
      geom = "line",
      group = 1
    ) +
    facet_wrap(as.formula(paste("~", facet.arg)),
               ncol = facet.ncol.arg,
               scales = "free_x")
    # geom_vline(data = dt.stim.arg,
    #   aes(xintercept = Stimulation_time - tfreq.arg),
    #   colour = rhg_cols[[3]],
    #   size = 0.25
    # ) +
  if(!is.null(dt.stim.arg)) {
    p.tmp = p.tmp + geom_segment(data = dt.stim.arg,
                                 aes(x = Stimulation_time - tfreq.arg,
                                     xend = Stimulation_time - tfreq.arg,
                                     y = ylim.arg[1],
                                     yend = ylim.arg[1] + abs(ylim.arg[2] - ylim.arg[1]) * stim.bar.height.arg),
                                 colour = rhg_cols[[3]],
                                 size = stim.bar.width.arg,
                                 group = 1) 
  }
  
  p.tmp = p.tmp + 
    scale_x_continuous(breaks = seq(0, maxrt.arg, xaxisbreaks.arg)) +
    coord_cartesian(ylim = ylim.arg) +
    xlab(paste0(xlab.arg, "\n")) +
    ylab(paste0("\n", ylab.arg)) +
    ggtitle(plotlab.arg) +
    theme_bw(base_size = 18, base_family = "Helvetica") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      axis.line.x = element_line(color = "black", size = 0.25),
      axis.line.y = element_line(color = "black", size = 0.25),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      strip.text.x = element_text(size = 14, face = "bold"),
      strip.text.y = element_text(size = 14, face = "bold"),
      strip.background = element_blank(),
      legend.key = element_blank(),
      legend.key.height = unit(1, "lines"),
      legend.key.width = unit(2, "lines"),
      legend.position = "top"
    )
  
  p.tmp
}


# Plots a scatter plot with marginal histograms
# Points are connected by a line (grouping by cellID)
#
# Assumes an input of data.table with
# x, y - columns with x and y coordinates
# id - a unique point identifier (here corresponds to cellID)
# mid - a (0,1) column by which points are coloured (here corresponds to whether cells are within bounds)

myGgplotScat = function(dt.arg,
                        a.arg = 1,
                        b.arg = 0,
                        band.arg = 0.5,
                        facet.arg = NULL,
                        facet.ncol.arg = 2,
                        xlab.arg = "x label",
                        ylab.arg = "y label",
                        plotlab.arg = "plot title",
                        alpha.arg = 1,
                        group.col.arg = NULL) {
  p.scat = ggplot(dt.arg, aes(x = x, y = y))
  
  if(is.null(group.col.arg)) {
    p.scat = p.scat +
      geom_point(alpha = alpha.arg)
      #geom_path(alpha = alpha.arg)
    
  } else {
    p.scat = p.scat +
      geom_point(aes(colour = as.factor(get(group.col.arg)), group = id), alpha = alpha.arg) +
      geom_path(aes(colour = as.factor(get(group.col.arg)), group = id), alpha = alpha.arg)
  }
  
  if(band.arg < 0) 
    p.scat = p.scat +
      stat_smooth(method = function(formula, data, weights = weight) rlm(formula, data, weights=weight, method='MM'),
                 fullrange=FALSE, level = 0.95, colour = 'red')
  else {
    p.scat = p.scat + 
      geom_abline(slope = a.arg, intercept = b.arg) +
      geom_abline(
        slope = a.arg,
        intercept = b.arg * (1 + band.arg),
        linetype = 'dashed'
      ) +
      geom_abline(
        slope = a.arg,
        intercept = b.arg * (1 - band.arg),
        linetype = 'dashed'
      )
  }

  if (!is.null(facet.arg)) {
    p.scat = p.scat +
      facet_wrap(as.formula(paste("~", facet.arg)),
               ncol = facet.ncol.arg)
    
  }
  
  p.scat = p.scat +
    xlab(paste0(xlab.arg, "\n")) +
    ylab(paste0("\n", ylab.arg)) +
    ggtitle(paste0(plotlab.arg, "\n")) +
    scale_color_manual(name = "", values = rhg_cols[c(7, 3)]) +
    theme_bw(base_size = 18, base_family = "Helvetica") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.line.x = element_line(color = "black", size = 0.25),
      axis.line.y = element_line(color = "black", size = 0.25),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      strip.text.x = element_text(size = 14, face = "bold"),
      strip.text.y = element_text(size = 14, face = "bold"),
      strip.background = element_blank(),
      legend.key = element_blank(),
      legend.key.height = unit(1, "lines"),
      legend.key.width = unit(2, "lines"),
      legend.position = "none"
    )
  
  if (is.null(facet.arg))
    ggExtra::ggMarginal(p.scat, type = "histogram",  bins = 100)
  else 
    p.scat
}


# ggplots with labels on top
myGplotBox = function(dt.arg, dt.label.arg = NULL, 
                      x.arg, y.arg, xlab.arg = "", ylab.arg = "",
                      plotlab.arg = "", scalelab.arg = "",
                      facet.arg, facet.ncol.arg = 3, col.arg, ylim.arg = c(0,1)) {
  p.loc = ggplot(dt.arg, aes_string(x = x.arg, y = y.arg)) +
    geom_boxplot(aes_string(fill = col.arg), outlier.shape  = NA) +
    facet_wrap(~ get(facet.arg), ncol = facet.ncol.arg) +
    coord_cartesian(ylim = ylim.arg) +
    xlab(paste0(xlab.arg, "\n")) +
    ylab(paste0("\n", ylab.arg)) +
    ggtitle(paste0(plotlab.arg, "\n")) +
    scale_fill_manual(name = scalelab.arg, values = rhg_cols[c(7, 3)]) +
    theme_bw(base_size = 18, base_family = "Helvetica") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.line.x = element_line(color = "black", size = 0.25),
      axis.line.y = element_line(color = "black", size = 0.25),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      strip.text.x = element_text(size = 14, face = "bold"),
      strip.text.y = element_text(size = 14, face = "bold"),
      strip.background = element_blank(),
      legend.key = element_blank(),
      legend.key.height = unit(1, "lines"),
      legend.key.width = unit(2, "lines"),
      legend.position = "top"
    )
  
  
  if (!is.null(dt.label.arg))
    p.loc = p.loc + geom_text(data = dt.label.arg, aes( y = ylim.arg[2] * 0.95, label = p.value.range), size = 12)
  
  p.loc
}



## Wilcoxon test

myTestWilcox = function(dt.arg, var.indep) {
  require(dplyr)
  require(data.table)
  require(broom)
  
  dt.res.test = dt.arg %>%
    group_by(RealTime, metadata.site.stim) %>%
    do(tidy(wilcox.test(val ~ get(var.indep), data = ., exact = TRUE, correct = TRUE, conf.int = TRUE)))
  
  dt.res.test = as.data.table(dt.res.test)
  dt.res.test[, p.value.range := ifelse(p.value < 0.001, "***", 
                                        ifelse(p.value < 0.01, "**", 
                                               ifelse(p.value < 0.05, "*", 
                                                      ifelse(p.value < 0.1, ".", " "))))]
}
