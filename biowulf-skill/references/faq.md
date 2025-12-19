# Frequently Asked Questions

## How do I acknowledge Biowulf?

In publications:
> *This work utilized the computational resources of the NIH HPC Biowulf cluster (https://hpc.nih.gov).*

## Why is my job pending?

Check reason with `squeue -u $USER` or `sjobs`:

| Reason | Meaning | Solution |
|--------|---------|----------|
| QOSMaxCpusPerUserLimit | Hit CPU limit | Wait or reduce job size |
| QOSJobLimit | Hit job count limit | Wait for jobs to finish |
| Resources | No available resources | Wait; check `freen` |
| Priority | Lower priority than others | Wait; use `--partition=quick` for short jobs |
| Dependency | Waiting on another job | Normal behavior |

### Speed Up Scheduling
- Request two partitions: `--partition=norm,quick`
- Request accurate resources (don't over-request)
- Bundle swarms with `-b` for short jobs
- Use `-p 2` for single-threaded swarms

## Why did my job get killed?

### Out of Memory
```bash
# Check with jobhist
jobhist JOBID
# Look at MemUsed vs MemReq

# Check dashboard
https://hpcnihapps.cit.nih.gov/auth/dashboard/
```

**Solution**: Resubmit with more memory (`--mem=Ng`)

### Time Limit
Job exceeded walltime.

**Solution**: Increase `--time` or optimize code

### Node Failure
Check `.e` file for `slurmstepd` errors.

**Solution**: Resubmit; consider `--requeue`

## Why is my home directory full?

Home directory quota: 16 GB (cannot be increased)

### Find What's Using Space
```bash
dust $HOME
```

### Common Causes
1. **Conda in home** - Move to `/data/$USER/conda`
2. **pip cache** - `rm -rf ~/.cache/pip`
3. **Singularity cache** - Move with `export SINGULARITY_CACHEDIR=/data/$USER/.singularity`
4. **Data files** - Move to `/data/$USER/`

## Can't find application/command

### Did you load the module?
```bash
module avail appname
module load appname
which appname
```

### Are you on Helix vs Biowulf?
Helix has limited software. SSH to biowulf.nih.gov for full access.

## 'module load' fails

### Check .bashrc for errors
```bash
cp /etc/skel/.bashrc ~/.bashrc
# Reconnect to Biowulf
```

### Module doesn't exist
```bash
module spider keyword  # Search
module avail appname   # Check exact name
```

## Graphics applications won't run

### Check environment
```bash
# Remove customizations temporarily
mv ~/.bashrc ~/.bashrc.bak
cp /etc/skel/.bashrc ~/.bashrc
module purge
# Reconnect and test
```

### Conda interfering?
Deactivate conda environments before graphics apps.

### Not enough memory?
Graphics apps may need more memory - increase `--mem`.

### Need GPU acceleration?
Use visual partition:
```bash
svis
```

## How do I resubmit a failed swarm?

### Get swarm info
```bash
jobhist JOBID
# Shows Swarm Command used
```

### Resubmit whole swarm
```bash
cd /path/to/submission/dir
swarm -f original.swarm [same options]
```

### Resubmit specific failed subjobs
```bash
# Find failed subjob scripts
ls /spin1/swarm/$USER/JOBID/

# Extract failed commands
cat /spin1/swarm/$USER/JOBID/cmd.3 /spin1/swarm/$USER/JOBID/cmd.5 > retry.swarm

# Resubmit with more resources if needed
swarm -f retry.swarm -g 16  # More memory
```

## How do I create a swarmfile programmatically?

### Bash loop
```bash
for f in /data/$USER/input/*.bam; do
    base=$(basename $f .bam)
    echo "samtools sort $f -o /data/$USER/output/${base}.sorted.bam"
done > commands.swarm
```

### Using find
```bash
find /data/$USER/input -name "*.fq" -exec echo "gzip {}" \; > compress.swarm
```

## Can I recover /lscratch data if job fails?

No automatic recovery. Add to script:
```bash
#!/bin/bash
#SBATCH --gres=lscratch:100

cd /lscratch/$SLURM_JOB_ID
# ... run commands ...

# If command fails, save lscratch
if [ $? -ne 0 ]; then
    tar -czf /data/$USER/rescue_${SLURM_JOB_ID}.tar.gz /lscratch/$SLURM_JOB_ID
fi
```

## Best practices for HPC

1. **Don't flood with short jobs** - Bundle if <15 min each
2. **Use /lscratch** - For temp files, set `TMPDIR=/lscratch/$SLURM_JOB_ID`
3. **< 1000 files per directory** - Use subdirectories
4. **Don't bloat .bashrc** - Keep startup files simple
5. **Load minimal modules** - Only what you need
6. **Test before scaling** - Run one job before submitting 1000

## Mac popup: "connect to server"

System Preferences → Users & Groups → Login Items
- Remove any HPC-related items
- Clear "Recent Servers" in Finder → Go → Connect to Server

## Getting help

### Email
staff@hpc.nih.gov

### Include in your email
- Are you on Biowulf or Helix?
- Job ID(s) of failed jobs
- Working directory (`pwd`)
- Modules loaded (`module list`)
- Exact command run
- Full error message (not screenshot)
- Path to error files

### Dashboard
https://hpcnihapps.cit.nih.gov/auth/dashboard/
- Job history
- Resource usage
- Error messages

## Common Error Messages

| Error | Likely Cause | Solution |
|-------|--------------|----------|
| `slurmstepd: Exceeded job memory limit` | Out of memory | Increase `--mem` |
| `CANCELLED DUE TO TIME LIMIT` | Job took too long | Increase `--time` |
| `command not found` | Module not loaded | `module load appname` |
| `Permission denied` | Wrong permissions | Check file permissions |
| `No space left on device` | Quota exceeded | Check `checkquota` |
| `cannot allocate memory` | System memory issue | Contact staff |
