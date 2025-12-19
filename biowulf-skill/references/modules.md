# Environment Modules Reference

Environment modules dynamically configure software environments on Biowulf.

## Basic Commands

| Command | Description |
|---------|-------------|
| `module avail` | List all available modules |
| `module avail appname` | List versions of specific app |
| `module load appname` | Load default version |
| `module load appname/version` | Load specific version |
| `module unload appname` | Unload a module |
| `module list` | Show loaded modules |
| `module purge` | Unload all modules |
| `module switch old new` | Switch versions |
| `module display appname` | Show what module does |
| `module spider keyword` | Search for modules |

## Finding Modules

### List All Available
```bash
module avail

# Show only default versions
module -d avail
```

### Search for Specific Software
```bash
# Case-insensitive search
module spider bowtie

# Regex search (useful for R, Python)
module -r avail '^R$'
module -r avail '^python'
```

### Check Available Versions
```bash
module avail samtools
# samtools/1.17  samtools/1.18  samtools/1.19(default)
```

## Loading Modules

### Load Default Version
```bash
module load samtools
```

### Load Specific Version
```bash
module load samtools/1.17
```

### Load Multiple Modules
```bash
module load bowtie2 samtools bcftools
```

### Check What's Loaded
```bash
module list
# Currently Loaded Modules:
#   1) samtools/1.19   2) bowtie2/2.5.1
```

## Unloading Modules

### Unload Specific Module
```bash
module unload samtools
```

### Unload All Modules
```bash
module purge
```

## Switching Versions

### Using switch
```bash
module load python/3.10
module switch python python/3.12
```

### Or just load new version
```bash
module load python/3.10
module load python/3.12  # Automatically unloads 3.10
```

## Examining Modules

### See What a Module Does
```bash
module display samtools
# Shows: PATH changes, environment variables set, etc.
```

### See Module Location
```bash
module show samtools
```

## Using Modules in Scripts

### Batch Script
```bash
#!/bin/bash
#SBATCH --mem=8g

module load samtools bowtie2
samtools --version
bowtie2 --version
```

### Swarm
```bash
# Using --module option
swarm -f commands.swarm --module samtools,bcftools

# Or in swarmfile
module load samtools && samtools index file1.bam
module load samtools && samtools index file2.bam
```

### Quiet Loading (suppress messages)
```bash
module -q load samtools  # No output messages
```

## Common Module Combinations

### Python with CUDA (for deep learning)
```bash
module load python/3.12 cuDNN/8.9.2/CUDA-12 CUDA/12.1
```

### Bioinformatics Pipeline
```bash
module load bowtie2 samtools bcftools bedtools
```

### R with Dependencies
```bash
module load R/4.5
```

## Personal Modulefiles

### Create Personal Module Directory
```bash
mkdir -p ~/modulefiles
```

### Add to Module Path
```bash
module use --append ~/modulefiles
# Or prepend (your modules take priority):
module use --prepend ~/modulefiles
```

### Create a Modulefile
```lua
-- File: ~/modulefiles/myapp/1.0.lua
whatis("My custom application")
prepend_path("PATH", "/data/$USER/apps/myapp/bin")
prepend_path("LD_LIBRARY_PATH", "/data/$USER/apps/myapp/lib")
setenv("MYAPP_HOME", "/data/$USER/apps/myapp")
```

## Shared Group Modules

### Setup
```bash
# Create shared modulefiles directory
mkdir /data/MyGroup/modulefiles

# Each user creates symlink
ln -s /data/MyGroup/modulefiles ~/modulefiles/shared

# Load personal modules
module load use.own
module avail  # Now shows shared modules
```

## Shell Compatibility

### Bash (default)
```bash
module load appname  # Works directly
```

### tcsh Scripts
```tcsh
#!/bin/tcsh
source /etc/profile.d/modules.csh
module load appname
```

## Troubleshooting

### Module Command Not Found
Check if `~/.bashrc` was modified incorrectly:
```bash
cp /etc/skel/.bashrc ~/.bashrc
# Then reconnect
```

### Application Not Found After Module Load
```bash
# Verify module is loaded
module list

# Check PATH
echo $PATH | tr ':' '\n' | grep appname

# See what module sets
module display appname
```

### Version Conflicts
```bash
# Unload all, reload what you need
module purge
module load app1 app2
```

### Python Module Conflicts
Some Python modules include tools like samtools. Loading order matters:
```bash
# Load specific tools AFTER Python if needed
module load python/3.12
module load samtools  # This samtools will be in PATH first
```

## Tips

1. **Check versions** - `module avail appname` before loading
2. **Be specific** - Use `appname/version` for reproducibility
3. **In scripts** - Always specify version for reproducibility
4. **Purge first** - `module purge` before loading for clean environment
5. **Use swarm --module** - Cleaner than module load in each command
