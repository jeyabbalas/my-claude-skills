# Swarm Reference

Swarm submits many similar commands as a Slurm job array.

## Basic Usage

### Create a Swarmfile
Each line is one command:
```bash
# file: commands.swarm
program input1.txt > output1.txt
program input2.txt > output2.txt
program input3.txt > output3.txt
```

### Submit
```bash
swarm -f commands.swarm
```

## Key Options

| Option | Description | Example |
|--------|-------------|---------|
| `-f FILE` | Swarmfile | `-f commands.swarm` |
| `-g N` | GB memory per command | `-g 8` |
| `-t N` | Threads per command | `-t 4` |
| `-b N` | Bundle N commands per subjob | `-b 10` |
| `-p N` | Pack N single-threaded commands per subjob | `-p 2` |
| `--time HH:MM:SS` | Walltime per command | `--time 2:00:00` |
| `--module name` | Load module(s) | `--module samtools,python` |
| `--partition name` | Partition | `--partition quick` |
| `--gres resource` | Generic resources | `--gres=lscratch:50` |
| `--logdir dir` | Output directory | `--logdir /data/$USER/logs` |
| `--merge-output` | Combine stdout/stderr | |
| `--devel` | Show what would run (don't submit) | |

## Common Patterns

### Multi-threaded Commands
```bash
# Each command uses 8 threads
swarm -f commands.swarm -g 16 -t 8 --module bowtie2
```

Swarmfile uses `$SLURM_CPUS_PER_TASK`:
```bash
bowtie2 -p $SLURM_CPUS_PER_TASK -x index -1 r1.fq -2 r2.fq -S out1.sam
bowtie2 -p $SLURM_CPUS_PER_TASK -x index -1 r3.fq -2 r4.fq -S out2.sam
```

### Using Local Scratch
```bash
swarm -f commands.swarm -g 8 --gres=lscratch:100 --module samtools
```

Swarmfile:
```bash
cd /lscratch/$SLURM_JOB_ID; cp /data/$USER/in1.bam .; samtools sort in1.bam -o s1.bam; cp s1.bam /data/$USER/
cd /lscratch/$SLURM_JOB_ID; cp /data/$USER/in2.bam .; samtools sort in2.bam -o s2.bam; cp s2.bam /data/$USER/
```

### Bundling Large Swarms
For >1000 commands or very short commands (<15 min):
```bash
# Bundle 50 commands per subjob (run sequentially)
swarm -f commands.swarm -b 50 -g 4
```

### Packing Single-Threaded Commands
Use both hyperthreads (2 CPUs = 1 core):
```bash
swarm -f commands.swarm -p 2 -g 2  # 2 commands run in parallel per subjob
```

### GPU Swarm
```bash
swarm -f gpu_commands.swarm --partition=gpu --gres=gpu:v100x:1 -g 32 -t 8
```

### Quick Partition (Jobs <4 hours)
```bash
swarm -f commands.swarm --time 1:00:00 --partition=quick -g 4
```

## Swarmfile Directives

Options can be in the swarmfile:
```bash
#SWARM -g 8 -t 4 --time 4:00:00
#SWARM --module samtools
command1 input1.txt
command2 input2.txt
```

## Output Files

Default naming: `swarm_JOBID_SUBJOBID.{o,e}`
```
swarm_12345_0.o  swarm_12345_0.e
swarm_12345_1.o  swarm_12345_1.e
```

Custom naming:
```bash
swarm -f cmd.swarm --job-name myswarm --logdir /data/$USER/logs
# Creates: myswarm_12345_0.o, etc.
```

## Monitoring

```bash
# Check swarm status
sjobs
squeue -u $USER

# Detailed job info
jobhist JOBID

# Cancel swarm
scancel JOBID
```

## Generating Swarmfiles

### Bash loop
```bash
for f in /data/$USER/input/*.fq; do
    base=$(basename $f .fq)
    echo "program $f > /data/$USER/output/${base}.out"
done > commands.swarm
```

### Using find
```bash
find /data/$USER/input -name "*.bam" -exec echo "samtools index {}" \; > index.swarm
```

## Environment Variables in Swarm

| Variable | Description |
|----------|-------------|
| `$SLURM_CPUS_PER_TASK` | CPUs allocated |
| `$SLURM_JOB_ID` | Job ID |
| `$SLURM_ARRAY_JOB_ID` | Array parent job ID |
| `$SLURM_ARRAY_TASK_ID` | Subjob index (0, 1, 2...) |
| `$SWARM_PROC_ID` | Process ID when using `-p` |

## Tips

1. **Test first**: Use `--devel` to see what would run
2. **Memory**: Check one job's memory usage, then set `-g` accordingly
3. **Bundling**: Auto-bundles if >1000 commands
4. **Walltime**: Multiplied by bundle factor automatically
5. **Comments**: Lines starting with `#` are ignored
6. **Bash required**: Use bash syntax (not csh) unless `--usecsh`

## Example: Complete Workflow

```bash
# 1. Create swarmfile
cat > align.swarm << 'EOF'
#SWARM -t 8 -g 16 --time 4:00:00
#SWARM --module bowtie2,samtools
bowtie2 -p $SLURM_CPUS_PER_TASK -x /fdb/igenomes/hg38 -1 s1_R1.fq -2 s1_R2.fq | samtools sort -o s1.bam
bowtie2 -p $SLURM_CPUS_PER_TASK -x /fdb/igenomes/hg38 -1 s2_R1.fq -2 s2_R2.fq | samtools sort -o s2.bam
bowtie2 -p $SLURM_CPUS_PER_TASK -x /fdb/igenomes/hg38 -1 s3_R1.fq -2 s3_R2.fq | samtools sort -o s3.bam
EOF

# 2. Preview
swarm -f align.swarm --devel

# 3. Submit
swarm -f align.swarm

# 4. Monitor
sjobs
```
