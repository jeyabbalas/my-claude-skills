# Job Submission Reference

## sbatch Command

### Basic Syntax
```bash
sbatch [options] jobscript.sh
```

### Essential Options

| Option | Description | Example |
|--------|-------------|---------|
| `--cpus-per-task=N` | CPUs per task | `--cpus-per-task=8` |
| `--mem=Ng` | Memory (note the 'g') | `--mem=32g` |
| `--time=HH:MM:SS` | Walltime limit | `--time=24:00:00` |
| `--partition=name` | Partition | `--partition=gpu` |
| `--gres=resource` | Generic resources | `--gres=gpu:a100:1` |
| `--exclusive` | Exclusive node access | `--exclusive` |
| `--output=file` | Stdout file | `--output=job_%j.out` |
| `--error=file` | Stderr file | `--error=job_%j.err` |
| `--job-name=name` | Job name | `--job-name=myanalysis` |

### Default Allocations
- Without options: 2 CPUs, 4 GB memory
- Memory scales: 2 GB per CPU requested

### Example Job Scripts

**Single-threaded job:**
```bash
#!/bin/bash
#SBATCH --mem=8g
#SBATCH --time=4:00:00

module load samtools
samtools sort input.bam -o sorted.bam
```

**Multi-threaded job:**
```bash
#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=32g
#SBATCH --time=8:00:00

module load bowtie2
bowtie2 -p $SLURM_CPUS_PER_TASK -x index -1 r1.fq -2 r2.fq -S out.sam
```

**GPU job:**
```bash
#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64g
#SBATCH --time=24:00:00

module load python/3.12 cuDNN/8.9.2/CUDA-12 CUDA/12.1
python train_model.py
```

## Partitions

| Partition | Purpose | Notes |
|-----------|---------|-------|
| `norm` | Default, single-node jobs | Most jobs go here |
| `multinode` | Multi-node parallel jobs | No single-node jobs allowed |
| `largemem` | High memory (>350GB) | Must request ≥350GB |
| `quick` | Jobs <4 hours | Higher priority |
| `gpu` | GPU jobs | Must specify GPU type |
| `visual` | Visualization | Hardware-accelerated graphics |
| `unlimited` | Jobs >10 days | Small partition, use sparingly |

### Using Multiple Partitions
```bash
sbatch --partition=norm,quick script.sh  # Max 2 partitions
```

## Walltimes

Check limits with `batchlim`:
```bash
batchlim
```

Set walltime:
```bash
sbatch --time=24:00:00 script.sh      # 24 hours
sbatch --time=5-00:00:00 script.sh    # 5 days
```

## Job Dependencies

```bash
# Run job2 after job1 completes (any exit status)
sbatch job1.sh  # Returns 12345
sbatch --dependency=afterany:12345 job2.sh

# Run only if previous succeeded
sbatch --dependency=afterok:12345 job2.sh

# Run only if previous failed
sbatch --dependency=afternotok:12345 job2.sh

# Multiple dependencies
sbatch --dependency=afterany:12345,12346 job3.sh
```

## Local Scratch Disk

Allocate with `--gres=lscratch:N` (N = GB):
```bash
#!/bin/bash
#SBATCH --gres=lscratch:100

# Use local scratch
cd /lscratch/$SLURM_JOB_ID
cp /data/$USER/input.bam .
# ... process data ...
cp output.bam /data/$USER/
```

**Important**: Data in `/lscratch/$SLURM_JOB_ID` is deleted when job ends!

## Job Management

```bash
# Check job status
squeue -u $USER
sjobs

# Cancel job
scancel JOBID
scancel -u $USER  # Cancel all your jobs

# Job history
jobhist JOBID
sacct -j JOBID

# Modify pending job
scontrol update JobId=JOBID TimeLimit=48:00:00
```

## Environment Variables in Jobs

| Variable | Description |
|----------|-------------|
| `$SLURM_JOB_ID` | Job ID |
| `$SLURM_CPUS_PER_TASK` | CPUs allocated |
| `$SLURM_MEM_PER_NODE` | Memory allocated |
| `$SLURM_ARRAY_TASK_ID` | Array job index |
| `$SLURM_ARRAY_JOB_ID` | Array job parent ID |

## Email Notifications

```bash
sbatch --mail-type=BEGIN,END,FAIL --mail-user=you@nih.gov script.sh

# Available types: BEGIN, END, FAIL, REQUEUE, ALL
# TIME_LIMIT_50, TIME_LIMIT_80, TIME_LIMIT_90
```

## Job Directives in Script

Place after `#!/bin/bash`:
```bash
#!/bin/bash
#SBATCH --job-name=myanalysis
#SBATCH --cpus-per-task=8
#SBATCH --mem=32g
#SBATCH --time=8:00:00
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
```
