####
## Analyzes objNuclei.csv in all sub-directories of output
## Extracts single-cell trajectories

# Execute in a directory where output is located


require(data.table)
require(dplyr)
require(broom)
require(xlsx)
require(RCurl)

# Source file with auxilary functions
#source('calcTraj_customFns.R')

# auxScript = getURL("https://www.dropbox.com/s/44jyl6ozdusyfu7/calcTraj_customFns.R?dl=0", ssl.verifypeer = FALSE)
# eval(parse(text = auxScript))

# s.sysname = Sys.info()["nodename"]
# if (s.sysname == "mattekPro2")
#  source(pipe(paste("wget -O -", "https://www.dropbox.com/s/44jyl6ozdusyfu7/calcTraj_customFns.R?dl=0"))) #else
  # source('~/scripts/r/calcTraj_customFns.R')


####
## Main parameters

# some constants
# these are columns created in this script

s.met.well = 'Metadata_Well' # Well column in experimental description is changed to this
s.met.sitestim = 'Metadata_SiteStim' # Column with merged site and stimulation conditions

# name of the file with paramaters for analysis and experiment
# this file should be located one folder up of cp.out from which the script is executed
s.par.plot = '../plotFormat.xlsx'

####
## Read parameters from plotFormat.xlsx file
l.par = myParRead(s.par.plot)

# Create directory for plots and output in the currenty working directory
ifelse(!dir.exists(file.path(".", l.par$dir.plot)), dir.create(file.path(".", l.par$dir.plot)), FALSE)


####
## Load experiment description

# The file with experimental description has name dependent on experiment
# look for all 'xlsx' files in a directory above current working directory
s.fname.tmp = list.files(
  path = '..',
  pattern = '*.xlsx',
  recursive = FALSE,
  full.names = TRUE
)

# take only those filenames that aren't equal to s.par.plot
s.par.exp = setdiff(s.fname.tmp, s.par.plot)

if (length(s.par.exp) > 1) {
  cat(file = stderr(), "Warning: More than one xlsx file with experiment description detected!")
  s.par.exp = s.par.exp[1]  
}


dt.exp = myExpRead(s.par.exp, inStartRow = 3)


# obtain stimulation time points from the experiment file
# this entire dt is then passed to plotting f-n
# Times provided in minutes!
# Searches for columns named:
# Stimulation_time_1...N
# stimulation_time_1..N, etc.
dt.t.stim = dt.exp[, lapply(.SD, median), .SDcols = names(dt.exp)[names(dt.exp) %like% "timulation.*time"]]

# check whether multipulse experiment
# b.multipulse flag set accordingly
if (ncol(dt.t.stim) > 1) {
  b.multipulse = TRUE
} else {
  b.multipulse = FALSE
  n.t.stim1 = dt.t.stim[[1]]
}


# Acquisition frequency from experimental description file
# Searches for a column that includes a string *requency*
# Frequency in minutes!
n.t.freq = 1
if (sum(names(dt.exp) %like% 'requency') > 0) {
  n.t.freq = as.numeric(dt.exp[1, names(dt.exp)[names(dt.exp) %like% 'requency'], with = FALSE])
}




#####
## Loading files

# search subdirectories for csv files
s.files.nuc = list.files(
  path = paste0(l.par$dir.out, '/.'),
  pattern = l.par$files.nuc,
  recursive = TRUE,
  full.names = TRUE
)


## Load main data files using custom file reading function
dt.nuc = do.call(rbind, lapply(s.files.nuc, myFread))




####
## Assignment of columns:
## s.metadata.site
## s.metadata.time
## s.trackObjectsLabel
## s.flint.nuc/cyt/corr/raw
## s.pos.x/y


s.met.site = l.par$metadata.site
s.met.time = l.par$metadata.time
s.met.trackabel = l.par$metadata.track


# extract columns with measurements of Erk-KTR, NucMem, position XY
l.met = myMetExtract(names(dt.nuc))

# select all columns with measurements;
# here from modules Intensity, AreaShape, TrackObjects
# Also, remove s.met.trackabel from this list because it's defined from plotFormat file
s.cols.meas = setdiff(names(dt.nuc)[grep('.*_Intensity_.*|.*_AreaShape_.*|.*_TrackObjects_.*|.*_Neighbors_.*', names(dt.nuc))],
                      s.met.trackabel)

# add position and ObjectNumber columns
s.cols.meas = c(s.cols.meas, l.met$pos.x, l.met$pos.y, l.met$objLabel)

# s.cols.meas stores columns relevant for further calculation
# These will be retained in the output of myTrajExtr f-n
# and saved to tCoursesSelected.csv


##### 
## Remove cells (not tracks!) without cytosol identified

dt.nuc = dt.nuc[get(l.met$flErk.mn.cyt.corr) > 0]

#####
## Extract uninterrupted timecourses
## max 1-frame break (can be many breaks)
## Trajectories that last entire experiment

dt.sel = myTrajExtr(in.dt = dt.nuc, 
                        in.max.break = 1, 
                        in.met.series = s.met.site, 
                        in.met.t = s.met.time, 
                        in.met.tracklabel = s.met.trackabel, 
                        in.aggr.cols = s.cols.meas)


# Add a column with real time
# New column is named RealTime
dt.sel = myAddRealTime(dt.sel, s.met.time, n.t.freq)


# Merge with experiment description
# Assumes that listed column exist in the xlsx file
dt.sel = merge(dt.sel, dt.exp[, c(
  'Position',
  'Well',
  'Stimulation_duration',
  'Stimulation_intensity',
  'Stimulation_treatment'
), with = FALSE],
by.x = s.met.site, by.y = 'Position')


# Change column name
setnames(dt.sel, "Well", s.met.well)


# number of cells per Site
dt.ncells.site = dt.sel[get(s.met.time) == min(dt.sel[[s.met.time]]), .N, by = c(s.met.site, 'Stimulation_duration', 'Stimulation_intensity', 'Stimulation_treatment')]
dt.ncells.well = dt.sel[get(s.met.time) == min(dt.sel[[s.met.time]]), .N, by = c(s.met.well, 'Stimulation_duration', 'Stimulation_intensity', 'Stimulation_treatment')]
dt.ncells.cond = dt.sel[get(s.met.time) == min(dt.sel[[s.met.time]]), .N, by = c('Stimulation_duration', 'Stimulation_intensity', 'Stimulation_treatment')]

if (l.par$plot.save) {
  pdf(paste0(l.par$dir.plot, "/tab_nCells_perSite.pdf"),
      height = 7,
      width = 10)
  grid.table(format(dt.ncells.site, digits = 3))
  dev.off()
  
  pdf(paste0(l.par$dir.plot, "/tab_nCells_perWell.pdf"),
      height = 7,
      width = 10)
  grid.table(format(dt.ncells.well, digits = 3))
  dev.off()

  pdf(paste0(l.par$dir.plot, "/tab_nCells_perCond.pdf"),
      height = 7,
      width = 10)
  grid.table(format(dt.ncells.well, digits = 3))
  dev.off()
}

setkeyv(dt.sel, c(s.met.site, s.met.trackabel, s.met.time))


#####
## save file with selected trajectories omitted
s.cols.omit = c(s.met.time, s.met.sitestim, 'TrackObjects_Label_uni')

if (l.par$plot.save) {
  write.csv(x = dt.sel[, setdiff(names(dt.sel), s.cols.omit), with = FALSE], 
            file = "tCoursesSelected.csv", row.names = FALSE)
}


#####
## Create and save file with stimuation sequence
## The file comes directly from experimental description.
## It's reshaped to fit plotting function
## Vertical bars under trajectories are plotted based on this dt
dt.t.stim.resh = reshape(dt.t.stim, dir = 'long', varying = list(1:ncol(dt.t.stim)), v.names = 'Stimulation_time')

if (l.par$plot.save) 
  write.csv(x = dt.t.stim.resh, file = 'stimT.csv', row.names = FALSE)




#####
## Plotting


p.out = list()

# Use facet.arg = s.met.sitestim to plot per FOV/site with facet labels 
#including FOV suber and the corresponding experimental condition
p.out$traj_ERK_cytoVnuc_noIllumCorr = myGgplotTraj(
  dt.arg = dt.sel,
  x.arg = "RealTime",
  y.arg = paste(l.met$flErk.mn.cyt.raw,  '/', l.met$flErk.mn.nuc.raw),
  group.arg = "TrackObjects_Label_uni",
  xlab.arg = "Time (min)",
  ylab.arg = "Erk-KTR: cytoplasmic vs nuclear (mean fl.int.)",
  plotlab.arg = "Raw data from uncorrected images",
  facet.arg = paste0(s.met.site, ' + Stimulation_duration + Stimulation_intensity + Stimulation_treatment'),
  dt.stim.arg = dt.t.stim.resh,
  tfreq.arg = n.t.freq,
  maxrt.arg = max(max(dt.sel[[s.met.time]]) + 60),
  xaxisbreaks.arg = l.par$plot.x.interv,
  facet.ncol.arg = l.par$plot.facets.ncol.site,
  ylim.arg = c(0, 1.2),
  stim.bar.height.arg = 0.05,
  stim.bar.width.arg = 1
)

p.out$traj_ERK_cytoVnuc_illumCorr = myGgplotTraj(
  dt.arg = dt.sel,
  x.arg = "RealTime",
  y.arg = paste(l.met$flErk.mn.cyt.corr,  '/', l.met$flErk.mn.nuc.corr),
  group.arg = "TrackObjects_Label_uni",
  xlab.arg = "Time (min)",
  ylab.arg = "Erk-KTR: cytoplasmic vs nuclear (mean fl.int.)",
  plotlab.arg = "Raw data from bg-corrected images",
  facet.arg = paste0(s.met.site, ' + Stimulation_duration + Stimulation_intensity + Stimulation_treatment'),
  facet.ncol.arg = l.par$plot.facets.ncol.site,
  dt.stim.arg = dt.t.stim.resh,
  stim.bar.height.arg = 0.05,
  stim.bar.width.arg = 1,
  tfreq.arg = n.t.freq,
  maxrt.arg = max(max(dt.sel[[s.met.time]]) + 60),
  xaxisbreaks.arg = l.par$plot.x.interv,
  ylim.arg = c(0, 1.2)
)


# You can combine different variables for plot facetting,
# e.g. to facet based on merged experimental conditions use
# facet.arg = 'Stimulation_duration + Stimulation_intensity + Stimulation_treatment'
p.out$traj_ERK_cytoVnuc_illumCorr_perCond = myGgplotTraj(
  dt.arg = dt.sel,
  x.arg = "RealTime",
  y.arg = paste(l.met$flErk.mn.cyt.corr,  '/', l.met$flErk.mn.nuc.corr),
  group.arg = "TrackObjects_Label_uni",
  xlab.arg = "Time (min)",
  ylab.arg = "Erk-KTR: cytoplasmic vs nuclear (mean fl.int.)",
  plotlab.arg = "Raw data from bg-corrected images",
  facet.arg = 'Stimulation_duration + Stimulation_intensity + Stimulation_treatment',
  facet.ncol.arg = l.par$plot.facets.ncol.cond,
  dt.stim.arg = dt.t.stim.resh,
  stim.bar.height.arg = 0.05,
  stim.bar.width.arg = 1,
  tfreq.arg = n.t.freq,
  maxrt.arg = max(max(dt.sel[[s.met.time]]) + 60),
  xaxisbreaks.arg = l.par$plot.x.interv,
  ylim.arg = c(0, 1.2)
)


p.out$traj_ERK_nuc_illumCorr = myGgplotTraj(
  dt.arg = dt.sel,
  x.arg = "RealTime",
  y.arg = l.met$flErk.mn.nuc.corr,
  group.arg = "TrackObjects_Label_uni",
  xlab.arg = "Time (min)",
  ylab.arg = "Erk-KTR:  nuclear (mean fl.int.)",
  plotlab.arg = "Raw data from bg-corrected images",
  facet.arg = paste0(s.met.site, ' + Stimulation_duration + Stimulation_intensity + Stimulation_treatment'),
  dt.stim.arg = dt.t.stim.resh,
  tfreq.arg = n.t.freq,
  maxrt.arg = max(max(dt.sel[[s.met.time]]) + 60),
  xaxisbreaks.arg = l.par$plot.x.interv,
  ylim.arg = c(0, 0.3),
  facet.ncol.arg = l.par$plot.facets.ncol.site,
  stim.bar.height.arg = 0.05,
  stim.bar.width.arg = 1
)



p.out$traj_ERK_nucInv_illumCorr = myGgplotTraj(
  dt.arg = dt.sel,
  x.arg = "RealTime",
  y.arg = paste("1 /", l.met$flErk.mn.nuc.corr),
  group.arg = "TrackObjects_Label_uni",
  xlab.arg = "Time (min)",
  ylab.arg = "Erk-KTR: 1 / nuclear (mean fl.int.)",
  plotlab.arg = "Raw data from bg-corrected images",
  facet.arg = paste0(s.met.site, ' + Stimulation_duration + Stimulation_intensity + Stimulation_treatment'),
  facet.ncol.arg = l.par$plot.facets.ncol.site,
  dt.stim.arg = dt.t.stim.resh,
  stim.bar.height.arg = 0.05,
  stim.bar.width.arg = 1,
  tfreq.arg = n.t.freq,
  maxrt.arg = max(max(dt.sel[[s.met.time]]) + 60),
  xaxisbreaks.arg = l.par$plot.x.interv,
  ylim.arg = c(0, 30)
)



## Save plots to files
if (l.par$plot.save) {
  lapply(names(p.out),
         function(x)
           ggsave(
             filename = paste0("plots/", x, ".pdf"),
             plot = p.out[[x]],
             width = l.par$plot.width,
             height = l.par$plot.width
           ))
}



