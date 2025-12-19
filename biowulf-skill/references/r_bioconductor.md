# R/Bioconductor Reference

## Loading R

```bash
module avail R
# R/4.3, R/4.4, R/4.5, etc.

module load R/4.5  # Load default/latest
```

## Interactive R Session

```bash
# Always use sinteractive (not on login node!)
sinteractive --mem=16g --cpus-per-task=4 --gres=lscratch:10
module load R
R
```

## Installing Personal Packages

### Default Library Location
Packages install to `/data/$USER/R/rhel8/X.X/` (where X.X is R version)

### Create Library Directory
```bash
mkdir -p /data/$USER/R/rhel8/4.5/
```

### Install Packages
```r
# CRAN packages
install.packages("tidyverse")

# Bioconductor packages
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install("DESeq2")

# From GitHub
devtools::install_github("user/package")
```

### Using pacman
```r
library(pacman)
p_install(ggplot2, dplyr, tidyr)
```

## Batch Job Script

```bash
#!/bin/bash
#SBATCH --mem=16g
#SBATCH --cpus-per-task=4
#SBATCH --time=4:00:00
#SBATCH --gres=lscratch:10

module load R/4.5

# Set TMPDIR to local scratch
export TMPDIR=/lscratch/$SLURM_JOB_ID

Rscript analysis.R
# or
R --no-echo --no-restore --no-save < analysis.R > output.txt
```

## Parallel R with 'parallel' Package

### Detecting CPUs Correctly
```r
# CORRECT: Use allocated CPUs
library(parallelly)
ncpus <- availableCores()

# OR manually:
ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "2"))

# WRONG: Detects all node CPUs
# ncpus <- parallel::detectCores()  # Don't use!
```

### mclapply Example
```r
library(parallel)
ncpus <- parallelly::availableCores()

# Parallel lapply
results <- mclapply(1:100, function(x) {
    sqrt(x)
}, mc.cores = ncpus)
```

### foreach with doParallel
```r
library(foreach)
library(doParallel)
library(doMC)

registerDoMC(cores = parallelly::availableCores())

results <- foreach(i = 1:100, .combine = c) %dopar% {
    sqrt(i)
}
```

## BiocParallel

### Set Workers Correctly
```r
library(BiocParallel)

# IMPORTANT: BiocParallel doesn't know about Slurm!
# Set workers explicitly:
register(MulticoreParam(workers = parallelly::availableCores()), default = TRUE)

# Or in your script/session:
options(MulticoreParam = quote(MulticoreParam(workers = parallelly::availableCores())))
```

## Swarm of R Jobs

### Swarmfile
```bash
Rscript /data/$USER/scripts/analysis.R sample1
Rscript /data/$USER/scripts/analysis.R sample2
Rscript /data/$USER/scripts/analysis.R sample3
```

### Submit
```bash
swarm -f r_jobs.swarm -g 8 -t 2 --module R/4.5
```

### With Command Line Arguments
```r
# In R script
args <- commandArgs(trailingOnly = TRUE)
sample_name <- args[1]
cat("Processing:", sample_name, "\n")
```

## Implicit Multithreading

R uses MKL for matrix operations. Control threads with:
```bash
# R module sets OMP_NUM_THREADS=1 by default
# Set higher if your code benefits from parallelized math:
export OMP_NUM_THREADS=4
```

Or in R:
```r
# For dist() function
.Internal(setMaxNumMathThreads(4))
.Internal(setNumMathThreads(4))
```

## AnnotationHub/ExperimentHub

These need proxy settings:
```r
# Set proxy
Sys.setenv(ANNOTATION_HUB_PROXY = Sys.getenv("http_proxy"))
Sys.setenv(EXPERIMENT_HUB_PROXY = Sys.getenv("http_proxy"))

# Or use options
setAnnotationHubOption("PROXY", Sys.getenv("http_proxy"))
setExperimentHubOption("PROXY", Sys.getenv("http_proxy"))
```

## Using renv for Reproducible Environments

```r
# Initialize renv in project
renv::init()

# Install packages
install.packages("tidyverse")

# Save state
renv::snapshot()

# Restore on another system
renv::restore()
```

## RStudio on Biowulf

Use Open OnDemand: https://hpcondemand.nih.gov → Interactive Apps → RStudio Server

Or via `svis` for visualization.

## Common Issues

### Package Installation Fails
```r
# Check library path
.libPaths()
# Should include /data/$USER/R/rhel8/X.X/

# Create if missing
dir.create("/data/$USER/R/rhel8/4.5", recursive = TRUE)
```

### rlang Version Conflicts
```bash
# Remove local rlang if causing issues
rm -rf /data/$USER/R/rhel8/4.5/rlang
```

### Migrate Packages to New R Version
```r
# Get packages from old version
old_pkgs <- installed.packages(lib.loc = "/data/$USER/R/rhel8/4.4")[,"Package"]

# Find missing in new version
new_pkgs <- installed.packages(lib.loc = "/data/$USER/R/rhel8/4.5/")[,"Package"]
to_install <- setdiff(old_pkgs, new_pkgs)

# Install
BiocManager::install(to_install)
```

### Graphics in Batch Jobs
```r
# Use non-interactive backend
options(bitmapType = 'cairo')

# For ggplot2
ggsave("plot.png", width = 8, height = 6, dpi = 300)
```

## Shiny Apps

```r
# In R script with tunnel
port <- as.integer(Sys.getenv("PORT1"))
shinyApp(ui, server, options = list(
    port = port, 
    launch.browser = FALSE, 
    host = "127.0.0.1"
))
```

```bash
sinteractive --tunnel --mem=8g
module load R
Rscript shiny_app.R
# Access at http://localhost:PORT1 after creating tunnel
```

## Example: Complete Batch Workflow

```bash
#!/bin/bash
#SBATCH --mem=32g
#SBATCH --cpus-per-task=8
#SBATCH --time=8:00:00
#SBATCH --gres=lscratch:50

module load R/4.5
export TMPDIR=/lscratch/$SLURM_JOB_ID

# Run analysis
Rscript - <<'EOF'
library(parallel)
library(DESeq2)

ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "2"))
cat("Using", ncpus, "CPUs\n")

# Your analysis here
# ...

cat("Done!\n")
EOF
```
