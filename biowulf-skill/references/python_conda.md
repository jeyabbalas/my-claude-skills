# Python & Conda Reference

## System Python Modules

### Available Versions
```bash
module avail python
# python/3.8, python/3.9, python/3.10, python/3.11, python/3.12
```

### Load Python
```bash
module load python/3.12  # Latest recommended
```

Includes: numpy, scipy, pandas, scikit-learn, matplotlib, pytorch, tensorflow

## Creating Personal Conda Environments

### Recommended: Use mamba_install Wrapper

```bash
# Start interactive session (required!)
sinteractive --mem=20g --gres=lscratch:20

# Install miniforge (first time only)
module load mamba_install
mamba_install
# Installs to /data/$USER/conda by default
# Creates activation script at ~/bin/myconda

# Activate conda
source myconda

# Create environment
mamba create -n myenv python=3.10 numpy pandas
conda activate myenv
```

### Manual Installation

```bash
sinteractive --mem=20g --gres=lscratch:20
cd /data/$USER
export TMPDIR=/lscratch/$SLURM_JOB_ID

# Download miniforge
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh

# Install to /data (NOT home!)
bash Miniforge3-Linux-x86_64.sh -p /data/$USER/conda -b
rm Miniforge3-Linux-x86_64.sh

# Activate
source /data/$USER/conda/etc/profile.d/conda.sh
conda activate base
```

## Managing Environments

### Create Environment
```bash
# With specific packages
mamba create -n myproject python=3.10 numpy scipy pandas matplotlib

# From environment file
mamba env create -f environment.yml

# Clone existing
mamba create -n newenv --clone oldenv
```

### Activate/Deactivate
```bash
conda activate myenv
conda deactivate
```

### Install Packages
```bash
# Using mamba (faster)
mamba install packagename
mamba install bioconda::samtools

# Using conda
conda install packagename

# Using pip (in activated env)
pip install packagename
```

### List Environments
```bash
conda info --envs
```

### Remove Environment
```bash
conda env remove -n myenv
```

## Important Configuration

### Channel Setup (for bioinformatics)
```bash
conda activate myenv
conda config --env --add channels bioconda
conda config --env --add channels conda-forge
conda config --env --set channel_priority strict
```

### Prevent Auto-Activation
**Do NOT** let conda modify your `.bashrc`. Remove or comment out:
```bash
# >>> conda initialize >>>
...
# <<< conda initialize <<<
```

Instead, activate manually with `source myconda` or `source /data/$USER/conda/etc/profile.d/conda.sh`

## Using in Batch Jobs

### sbatch Script
```bash
#!/bin/bash
#SBATCH --mem=16g
#SBATCH --cpus-per-task=4

source /data/$USER/conda/etc/profile.d/conda.sh
conda activate myenv
python myscript.py
```

### Swarm
```bash
# In swarmfile
source /data/$USER/conda/etc/profile.d/conda.sh && conda activate myenv && python script.py arg1
source /data/$USER/conda/etc/profile.d/conda.sh && conda activate myenv && python script.py arg2
```

Or activate before submitting swarm:
```bash
source myconda
conda activate myenv
swarm -f commands.swarm
```

## Common Issues

### "403 error" for conda channels
Add to `/data/$USER/conda/.condarc`:
```yaml
channels:
  - conda-forge
  - bioconda
defaults: []
channel_priority: strict
```

### Home directory full
Install conda to `/data/$USER/conda`, not home directory

### VNC/graphics problems
Don't activate conda in `.bashrc` - can break desktop sessions

### OMP_NUM_THREADS
System Python sets `OMP_NUM_THREADS=1` to prevent overloading.
Set manually if needed:
```bash
export OMP_NUM_THREADS=8
```

### matplotlib display issues
For batch jobs without display:
```python
import matplotlib
matplotlib.use('agg')  # Before importing pyplot
import matplotlib.pyplot as plt
```

Or set environment:
```bash
export MPLBACKEND=agg
```

## Multiprocessing

### Detecting CPUs Correctly
```python
import os
# CORRECT: Use allocated CPUs
ncpus = int(os.environ.get('SLURM_CPUS_PER_TASK', '2'))

# WRONG: Detects all node CPUs
# import multiprocessing
# ncpus = multiprocessing.cpu_count()  # Don't use this!
```

### With multiprocessing Pool
```python
import os
import multiprocessing
import signal

def init_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)

ncpus = int(os.environ.get('SLURM_CPUS_PER_TASK', '2'))
pool = multiprocessing.Pool(ncpus, init_worker)
results = pool.map(my_function, data)
pool.close()
pool.join()
```

## Environment Files

### Export Current Environment
```bash
conda env export > environment.yml
# Or without build strings (more portable)
conda env export --no-builds > environment.yml
```

### Create from File
```bash
mamba env create -f environment.yml
```

### Example environment.yml
```yaml
name: myproject
channels:
  - conda-forge
  - bioconda
dependencies:
  - python=3.10
  - numpy
  - pandas
  - scikit-learn
  - pip
  - pip:
    - some-pip-package
```
