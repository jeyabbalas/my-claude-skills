# Jupyter Notebooks Reference

## Recommended: Open OnDemand

The easiest way to use Jupyter on Biowulf:

1. Go to https://hpcondemand.nih.gov
2. Click "Interactive Apps" → "Jupyter"
3. Configure resources (CPUs, memory, time, partition)
4. Click "Launch"
5. When ready, click "Connect to Jupyter"

## Manual Setup via SSH Tunnel

### Step 1: Start Interactive Session with Tunnel
```bash
# Use tmux to preserve session if disconnected
module load tmux
tmux

# Start interactive session with tunnel
sinteractive --gres=lscratch:5 --mem=10g --tunnel

# Output shows something like:
# ssh -L 33327:localhost:33327 biowulf.nih.gov
```

### Step 2: Start Jupyter
```bash
module load jupyter

# Start Jupyter Lab (recommended)
jupyter lab --ip localhost --port $PORT1 --no-browser

# Or Jupyter Notebook (classic)
jupyter notebook --ip localhost --port $PORT1 --no-browser

# Note the URL with token, e.g.:
# http://localhost:33327/lab?token=xxxxx
```

### Step 3: Create Tunnel from Local Machine

**Mac/Linux:**
```bash
ssh -L 33327:localhost:33327 biowulf.nih.gov
# Keep this connection open
```

**Windows PowerShell:**
```powershell
ssh -L 33327:localhost:33327 user@biowulf.nih.gov
```

### Step 4: Open Browser
Navigate to the URL shown by Jupyter (including the token)

## GPU Jupyter Session

```bash
sinteractive --partition=gpu --gres=gpu:a100:1,lscratch:50 \
    --mem=64g --cpus-per-task=8 --tunnel

module load jupyter
jupyter lab --ip localhost --port $PORT1 --no-browser
```

## Available Kernels

```bash
jupyter kernelspec list

# Common kernels:
# py3.10, py3.11, py3.12  - Python with scientific stack
# ir43, ir44              - R kernels
# bash                    - Bash kernel
```

## Using Custom Conda Environment as Kernel

### Install ipykernel in Your Environment
```bash
conda activate myenv
conda install ipykernel
# or: pip install ipykernel
```

### Register Kernel
```bash
conda activate myenv
python -m ipykernel install --user --name myenv --display-name "Python (myenv)"
```

### Verify
```bash
jupyter kernelspec list
# Should show ~/.local/share/jupyter/kernels/myenv
```

### Remove Kernel
```bash
jupyter kernelspec uninstall myenv
```

## Using R Magic in Python Notebooks

```bash
# Load rpy2 before starting Jupyter
module load jupyter rpy2
jupyter lab --ip localhost --port $PORT1 --no-browser
```

In notebook:
```python
%load_ext rpy2.ipython

%%R
x <- c(1, 2, 3)
print(x)
```

## Common Issues

### "jupyter-lab not found"
Don't load a Python module alongside Jupyter:
```bash
# WRONG:
module load python jupyter

# CORRECT:
module load jupyter
```

### PermissionError: /run/user/xxxx
```bash
unset XDG_RUNTIME_DIR
```

### Missing packages in kernel
Use `py3.10`, `py3.11`, or `py3.12` kernels (not just "Python 3") - these have the scientific stack.

### Export to PDF fails
```bash
module load tex
```
Then export from Jupyter.

## Jupyter Configuration

### Check Running Servers
```bash
jupyter server list
```

### Custom Config
```bash
jupyter lab --generate-config
# Creates ~/.jupyter/jupyter_lab_config.py
```

## Tips

1. **Use Open OnDemand** - much simpler than manual tunneling
2. **Use tmux/screen** - prevents losing session if network drops
3. **Request lscratch** - for faster local data access
4. **Don't run on login node** - always use sinteractive

## Example: Complete GPU Jupyter Setup

```bash
# 1. Start session in tmux
ssh biowulf.nih.gov
module load tmux
tmux

# 2. Get GPU node with tunnel
sinteractive --partition=gpu --gres=gpu:v100x:1,lscratch:100 \
    --mem=64g --cpus-per-task=8 --tunnel --time=8:00:00

# 3. Note the port and tunnel command shown

# 4. Load modules and start Jupyter
module load jupyter
jupyter lab --ip localhost --port $PORT1 --no-browser

# 5. In SEPARATE terminal on your local machine:
ssh -L 45000:localhost:45000 biowulf.nih.gov
# (use the actual port from step 3)

# 6. Open browser to URL shown by Jupyter

# 7. In Jupyter, verify GPU:
# import torch
# print(torch.cuda.is_available())
# print(torch.cuda.get_device_name(0))
```

## Reconnecting

If your local tunnel drops:
```bash
# On biowulf (in another terminal)
reconnect_tunnels
# Shows command to recreate tunnel
```

If your whole session drops:
```bash
# Reconnect to tmux
ssh biowulf.nih.gov
tmux attach
# Your Jupyter should still be running
```
