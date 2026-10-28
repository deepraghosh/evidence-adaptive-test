# EA: EDIT FIVE VALUES, SIMULATE, THEN ASSESS INDIVIDUAL GRAPHS.
# In RStudio run SETUP once. Edit one vector, select its simulation block S1--S6,
# then select either/both corresponding FIGURE blocks. One simulation produces
# BOTH the known- and estimated-covariance results; do not rerun it for the second
# plot. Selecting this ENTIRE file runs all six sweeps and displays all 12 figures.
# No PNG is saved automatically. In RStudio use Plots > Export to save a graph
# after inspecting it; use the plot-history arrows to return to an earlier graph.
#
# budget="quick": B=199,R=200 (check code/shape); "pilot": B=1999,R=2000;
# "full": B=19999,R=10000 (more precise, slower). Optional B/R override either.
# R objects include power, MCSE, Wilson intervals, paired comparisons and seeds.
# The estimated panel uses the manuscript's oracle common-bank calibration.
# Results stay in R memory. No simulation files or resumable caches are saved.
# Closing R loses these results; run a simulation once, then draw its two panels.
# Pilot correlation/dimension sweeps use 19,995 fresh experiments each. Repeated
# precision estimation can be slow, especially at large d with backend="base".
# Use budget="quick" for a first look; its power estimates are much noisier.
# Rerunning an S block repeats simulations; rerunning a FIGURE block only redraws.

# ---- SETUP: select and run this section once ----
# Keep all four .R files together in any local folder, including your Desktop.
# In RStudio open this file, then choose Session > Set Working Directory >
# To Source File Location before running SETUP. No internet is needed.
# Locate this saved script when sourced/Rscript; selected RStudio lines use getwd().
manuscript_dir <- local({
  paths <- unlist(lapply(sys.frames(), function(f) f$ofile), use.names = FALSE)
  if (length(paths)) {
    dirname(normalizePath(tail(paths, 1L), mustWork = TRUE))
  } else {
    arg <- grep("^--file=", commandArgs(), value = TRUE)
    if (length(arg)) {
      path <- sub("^--file=", "", arg[1L])
      if (!file.exists(path)) path <- gsub("~+~", " ", path, fixed = TRUE)
      dirname(normalizePath(path, mustWork = TRUE))
    } else getwd()
  }
})
if(!file.exists(file.path(manuscript_dir,"EA_manuscript_simulations.R")))
  stop("Set the working directory to the folder containing the four EA R files.")
engine_file <- normalizePath(file.path(manuscript_dir,"EA_manuscript_simulations.R"))
message("Loading simulation engine: ",engine_file)
source(engine_file)
.ea_runtime$source_file <- engine_file
if(!all(c("line_width","point_size") %in% names(formals(plot_manuscript_curves))))
  stop("An older engine was loaded from ",engine_file,". Replace it with the updated R file.")
if(!"workers" %in% names(formals(run_setting)))
  stop("An older engine without parallel replicates was loaded. Replace the engine and rerun SETUP.")
custom_root <- NULL                     # MEMORY ONLY: no result folders/files
budget <- "quick"                       # "quick", "pilot" or "full"
workers <- 4L                          # 4 workers + coordinating R = at most 5
# Use workers <- 1L for serial execution. All six S blocks use the same limit.
# Each setting splits its replicate scoring across workers; no nested pools.
backend <- "base"                       # standard R; no compiler files
seed_base <- 2026101500
B <- NULL                               # leave NULL for budget defaults
R <- NULL

# PLOT APPEARANCE: edit/run these lines, then run a FIGURE block only.
# This changes ALL method curves. No S1--S6 simulation needs to be rerun.
line_width <- 4                         # try 6 for even thicker lines
point_size <- 1.5                       # marker size for ALL methods; try 1.6

# EDIT THESE SIX VECTORS: exactly five distinct, finite values per vector.
main_correlations <- c(.5,.6,.7,.8,.9)  # -1 < rho < 1; zero/negative allowed
main_dimensions <- c(20,50,100,150,200) # integers >=4
main_strengths <- c(2.5,3,3.5,4,4.5)    # nonnegative; zero is permitted
original_correlations <- c(.1,.3,.5,.7,.9) # first five-block needs -.25 < rho <1
original_dimensions <- c(50,100,150,200,250) # integers >=10; off-grid is allowed
original_strengths <- c(1,1.5,2,2.5,3)   # nonnegative; zero is permitted

# FIXED PARAMETERS used while correlation or strength varies.
fixed_d <- 50
fixed_n <- 249                          # Wishart degrees, never a sample-size label
main_fixed_rho <- .95
main_fixed_strength <- 3.5
original_fixed_rho <- .9
original_fixed_strength <- 2
original_dimension_strength <- 3       # mean=(a,.5,...,.5,0,...), nine .5 entries
original_dimension_rho <- .7            # first ten-coordinate block; -.111...<rho<1

# DIMENSION SWEEPS derive n=5*d-1 (main), n=3*d-1 (original).
# Optional override: one common n, or five integer degrees matching the vector.
main_dimension_n <- NULL
original_dimension_n <- NULL
# In the original dimension model, first block size=10 (default rho=.7); background blocks
# have rho=.1,size=10, with a SHORTER LAST block if d is not a multiple of ten.

# Select one FIGURE block to display its graph in RStudio's Plots pane.
# Save it manually with Plots > Export when you are satisfied with it.
draw_custom_figure <- function(number,result) {
  plot_manuscript_figure(number,results=result,line_width=line_width,point_size=point_size)
}

# ---- S1: SIMULATE MAIN CORRELATIONS once; supplies FIGURES 1 and 4 ----
main_rho_result <- run_custom_sweep("rho",main_correlations,study="main",
  d=fixed_d,n=fixed_n,amplitude=main_fixed_strength,budget=budget,B=B,R=R,
  workers=workers,seed=seed_base,outdir=custom_root,backend=backend)

# ---- S2: SIMULATE MAIN DIMENSIONS once; supplies FIGURES 2 and 5 ----
main_d_result <- run_custom_sweep("d",main_dimensions,study="main",
  rho=main_fixed_rho,amplitude=main_fixed_strength,dimension_n=main_dimension_n,
  budget=budget,B=B,R=R,workers=workers,seed=seed_base+100,
  outdir=custom_root,backend=backend)

# ---- S3: SIMULATE MAIN STRENGTHS once; supplies FIGURES 3 and 6 ----
main_amplitude_result <- run_custom_sweep("amplitude",main_strengths,study="main",
  d=fixed_d,n=fixed_n,rho=main_fixed_rho,budget=budget,B=B,R=R,
  workers=workers,seed=seed_base+200,outdir=custom_root,backend=backend)

# ---- S4: SIMULATE ORIGINAL CORRELATIONS once; supplies FIGURES 7 and 10 ----
original_rho_result <- run_custom_sweep("rho",original_correlations,study="original",
  d=fixed_d,n=fixed_n,amplitude=original_fixed_strength,budget=budget,B=B,R=R,
  workers=workers,seed=seed_base+300,outdir=custom_root,backend=backend)

# ---- S5: SIMULATE ORIGINAL DIMENSIONS once; supplies FIGURES 8 and 11 ----
original_d_result <- run_custom_sweep("d",original_dimensions,study="original",
  rho=original_dimension_rho,amplitude=original_dimension_strength,dimension_n=original_dimension_n,
  budget=budget,B=B,R=R,workers=workers,seed=seed_base+400,
  outdir=custom_root,backend=backend)

# ---- S6: SIMULATE ORIGINAL STRENGTHS once; supplies FIGURES 9 and 12 ----
original_amplitude_result <- run_custom_sweep("amplitude",original_strengths,study="original",
  d=fixed_d,n=fixed_n,rho=original_fixed_rho,budget=budget,B=B,R=R,
  workers=workers,seed=seed_base+500,outdir=custom_root,backend=backend)

# ---- FIGURE 1: main known correlation (requires S1) ----
draw_custom_figure(1,main_rho_result)

# ---- FIGURE 2: main known dimension (requires S2) ----
draw_custom_figure(2,main_d_result)

# ---- FIGURE 3: main known amplitude (requires S3) ----
draw_custom_figure(3,main_amplitude_result)

# ---- FIGURE 4: main estimated correlation (requires S1) ----
draw_custom_figure(4,main_rho_result)

# ---- FIGURE 5: main estimated dimension (requires S2) ----
draw_custom_figure(5,main_d_result)

# ---- FIGURE 6: main estimated amplitude (requires S3) ----
draw_custom_figure(6,main_amplitude_result)

# ---- FIGURE 7: original known correlation (requires S4) ----
draw_custom_figure(7,original_rho_result)

# ---- FIGURE 8: original known dimension (requires S5) ----
draw_custom_figure(8,original_d_result)

# ---- FIGURE 9: original known amplitude (requires S6) ----
draw_custom_figure(9,original_amplitude_result)

# ---- FIGURE 10: original estimated correlation (requires S4) ----
draw_custom_figure(10,original_rho_result)

# ---- FIGURE 11: original estimated dimension (requires S5) ----
draw_custom_figure(11,original_d_result)

# ---- FIGURE 12: original estimated amplitude (requires S6) ----
draw_custom_figure(12,original_amplitude_result)

# Direct display without saving: plot_manuscript_figure(4,results=main_rho_result).
# While R stays open, redraw a graph from its in-memory result without simulating:
# plot_manuscript_figure(1,main_rho_result)

# ---- END: reset subsequent simulations to serial execution ----
# Every S block already stops its local workers automatically on completion/error.
# Run this line after your chosen block(s); set workers <- 4L to parallelize again.
workers <- 1L
