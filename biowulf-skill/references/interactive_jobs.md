# Interactive Jobs & Tunneling Reference

## sinteractive Command

### Basic Usage
```bash
sinteractive                    # Default: 2 CPUs, 1.5 GB memory
sinteractive --cpus-per-task=4  # 4 CPUs
sinteractive --mem=10g          # 10 GB memory
sinteractive --time=8:00:00     # 8 hour session
```

### Interactive GPU Sessions
```bash
# A100 GPU
sinteractive --partition=gpu --gres=gpu:a100:1 --mem=64g --cpus-per-task=8

# V100x GPU (32GB VRAM)
sinteractive --partition=gpu --gres=gpu:v100x:1 --mem=64g --cpus-per-task=8

# P100 GPU
sinteractive --partition=gpu --gres=gpu:p100:1 --mem=32g --cpus-per-task=8

# Multiple GPUs
sinteractive --partition=gpu --gres=gpu:a100:4 --mem=128g --cpus-per-task=32
```

### With Local Scratch
```bash
sinteractive --gres=lscratch:50 --mem=16g
cd /lscratch/$SLURM_JOB_ID
```

### Interactive Session Limits
- Max concurrent: 2 sessions
- Max walltime: 36 hours
- Check with `batchlim`

## SSH Tunneling

### Purpose
Connect to applications running on compute nodes (Jupyter, TensorBoard, etc.)

### Setup with sinteractive
```bash
sinteractive --tunnel  # Creates one tunnel
sinteractive -TT       # Creates two tunnels

# Output shows:
# Created 1 generic SSH tunnel(s)...
# Please create a SSH tunnel from your workstation:
#     ssh -L 45000:localhost:45000 biowulf.nih.gov
```

### Creating the Tunnel

**From Mac/Linux terminal:**
```bash
ssh -L 45000:localhost:45000 biowulf.nih.gov
# Keep this connection open
```

**From Windows PowerShell:**
```powershell
ssh -L 45000:localhost:45000 user@biowulf.nih.gov
```

**From Windows PuTTY:**
1. Session → Host: biowulf.nih.gov
2. SSH → Tunnels → Source port: 45000, Destination: localhost:45000
3. Click Add, then Open

### Using $PORT1 Variable
```bash
# After sinteractive --tunnel
echo $PORT1  # Shows the assigned port number

# Start application on that port
jupyter lab --ip localhost --port $PORT1 --no-browser
```

### Reconnecting Tunnels
```bash
# On biowulf, get tunnel command for existing sessions
reconnect_tunnels
```

## Visualization Jobs (svis)

For hardware-accelerated graphics:
```bash
svis  # Allocates visual partition node

# Follow instructions to set up VNC tunnel
# Then connect with VNC client
```

## Maintaining Sessions

### Using tmux (Recommended)
```bash
# On login node
module load tmux
tmux

# Start sinteractive inside tmux
sinteractive --tunnel

# Detach: Ctrl+b, then d
# Reconnect later:
tmux attach
```

### Using screen
```bash
screen
sinteractive
# Detach: Ctrl+a, then d
# Reconnect: screen -r
```

## GPU Monitoring

```bash
# On GPU node
nvidia-smi           # GPU status
nvidia-smi -l 1      # Update every second
watch nvidia-smi     # Continuous monitoring
```

## Example: Interactive PyTorch Session

```bash
# 1. Start GPU session with tunnel
sinteractive --partition=gpu --gres=gpu:a100:1,lscratch:50 \
    --mem=64g --cpus-per-task=8 --tunnel

# 2. Load modules
module load python/3.12 cuDNN/8.9.2/CUDA-12 CUDA/12.1

# 3. Test GPU
python -c "import torch; print(torch.cuda.is_available())"

# 4. Optional: Start TensorBoard
tensorboard --logdir=./logs --port $PORT1 --bind_all

# 5. In separate terminal, create tunnel
ssh -L 45000:localhost:45000 biowulf.nih.gov

# 6. Access TensorBoard at http://localhost:45000
```

## Tips

- Always use `sinteractive` for computation (not on login node)
- Request realistic resources to reduce wait time
- Use `--exclusive` sparingly (wastes resources)
- Set TMPDIR for programs that use temp files:
  ```bash
  export TMPDIR=/lscratch/$SLURM_JOB_ID
  ```
