---
name: biowulf
description: Comprehensive guide for NIH Biowulf HPC cluster. Use when users need help with job submission (sbatch, swarm), interactive sessions (sinteractive), GPU allocation, Python/conda environments, PyTorch deep learning, R/Bioconductor, Jupyter notebooks, data transfer (Globus), SSH tunneling, environment modules, or any Biowulf/HPC-related questions. Covers connecting, storage, partitions, walltimes, and troubleshooting.
---

# NIH Biowulf HPC Skill

Biowulf is the NIH High Performance Computing (HPC) cluster with ~93,000 cores and ~1,000 GPUs. This skill provides guidance for using Biowulf effectively.

## Quick Reference

### Connecting
```bash
ssh biowulf.nih.gov  # Use NIH username/password
```

### Basic Job Submission
```bash
# Single batch job
sbatch jobscript.sh

# Swarm of jobs (many similar commands)
swarm -f commands.swarm -g 4 -t 2 --module appname

# Interactive session
sinteractive --mem=10g --cpus-per-task=4
sinteractive --gres=gpu:v100x:1 --partition=gpu  # GPU session
```

### Storage Locations
| Location | Purpose | Quota |
|----------|---------|-------|
| `/home/$USER` | Config files, small scripts | 16 GB |
| `/data/$USER` | Data, conda envs, analysis | 100+ TB shared |
| `/lscratch/$SLURM_JOB_ID` | Fast local scratch (allocate with `--gres=lscratch:N`) | Per-job |

### Key Environment Variables
```bash
$SLURM_CPUS_PER_TASK  # CPUs allocated
$SLURM_JOB_ID         # Job ID
$SLURM_MEM_PER_NODE   # Memory allocated
```

## Reference Files

Read the appropriate reference file based on the user's question:

| Topic | Reference File | When to Read |
|-------|----------------|--------------|
| Hardware specs, nodes, GPUs | `references/hardware.md` | Questions about cluster capacity, node types, GPU types |
| Job submission (sbatch, partitions) | `references/job_submission.md` | Questions about submitting jobs, partitions, walltimes, dependencies |
| Interactive sessions & tunneling | `references/interactive_jobs.md` | sinteractive, GPU allocation, SSH tunnels, visualization |
| Swarm (parallel jobs) | `references/swarm.md` | Running many similar commands, job arrays |
| Python & conda environments | `references/python_conda.md` | Python setup, creating conda envs, pip, virtual environments |
| PyTorch deep learning | `references/deep_learning.md` | GPU training, PyTorch examples, CUDA modules |
| R/Bioconductor | `references/r_bioconductor.md` | R jobs, parallelization, package installation |
| Jupyter notebooks | `references/jupyter.md` | Running Jupyter on compute nodes, kernels |
| Data transfer (Globus) | `references/globus.md` | Moving data to/from Biowulf, Globus setup |
| Environment modules | `references/modules.md` | Loading software, module commands |
| FAQ & troubleshooting | `references/faq.md` | Common problems, job failures, debugging |

## Common Workflows

### Submit a batch job with specific resources
```bash
sbatch --cpus-per-task=8 --mem=32g --time=4:00:00 script.sh
```

### Run GPU job (e.g., PyTorch)
```bash
sinteractive --gres=gpu:a100:1 --partition=gpu --mem=64g --cpus-per-task=8
module load python/3.12 cuDNN/8.9.2/CUDA-12 CUDA/12.1
python train.py
```

### Create conda environment
```bash
sinteractive --mem=20g --gres=lscratch:20
module load mamba_install
mamba_install  # First time only
source myconda
mamba create -n myenv python=3.10 numpy pytorch
```

### Run Jupyter notebook
Use Open OnDemand: https://hpcondemand.nih.gov
Or via tunneling:
```bash
sinteractive --tunnel --mem=10g
module load jupyter
jupyter lab --ip localhost --port $PORT1 --no-browser
# Then create SSH tunnel from local machine
```

## Important Notes

- **No computation on login node**: Always use `sinteractive` or `sbatch`
- **Default allocation**: 2 CPUs, 4 GB memory if not specified
- **Default walltime**: Check with `batchlim` command
- **Citation**: "This work utilized the computational resources of the NIH HPC Biowulf cluster (https://hpc.nih.gov)"
- **Help**: staff@hpc.nih.gov

## Key Commands

| Command | Purpose |
|---------|---------|
| `module load appname` | Load software |
| `module avail appname` | List available versions |
| `freen` | Show free nodes/resources |
| `sjobs` | Show your running jobs |
| `jobhist JOBID` | Job history/resource usage |
| `checkquota` | Check disk quotas |
| `batchlim` | Show partition limits |
