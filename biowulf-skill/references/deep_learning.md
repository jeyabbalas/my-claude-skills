# Deep Learning (PyTorch) Reference

## GPU Types Available

| GPU | VRAM | Nodes | SLURM gres |
|-----|------|-------|------------|
| A100 | 80 GB | 76 | `gpu:a100:N` |
| V100-SXM2 | 32 GB | 56 | `gpu:v100x:N` |
| V100 | 16 GB | 8 | `gpu:v100:N` |
| P100 | 16 GB | 48 | `gpu:p100:N` |
| K80 | 24 GB | 72 | `gpu:k80:N` |

**Note**: Python 3.11+ modules are NOT compatible with K80 GPUs (compute capability <5.0). Use Python 3.10 or earlier for K80.

## Loading PyTorch

### Interactive Session
```bash
# Request A100 GPU
sinteractive --partition=gpu --gres=gpu:a100:1 --mem=64g --cpus-per-task=8

# Load modules
module load python/3.12 cuDNN/8.9.2/CUDA-12 CUDA/12.1

# Verify GPU
python -c "import torch; print(torch.cuda.is_available())"
python -c "import torch; print(torch.cuda.get_device_name(0))"
```

### For K80 GPUs (older)
```bash
sinteractive --partition=gpu --gres=gpu:k80:1 --mem=32g --cpus-per-task=8
module load python/3.10  # 3.10 or earlier for K80
```

## Batch Job Script

```bash
#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64g
#SBATCH --time=24:00:00
#SBATCH --gres=lscratch:100

# Load modules
module load python/3.12 cuDNN/8.9.2/CUDA-12 CUDA/12.1

# Use local scratch for data
export TMPDIR=/lscratch/$SLURM_JOB_ID
cd $TMPDIR

# Copy data to local scratch (faster I/O)
cp -r /data/$USER/dataset ./

# Run training
python train.py --data ./dataset --epochs 100

# Copy results back
cp -r ./checkpoints /data/$USER/results/
```

## PyTorch Example: MNIST

```python
import torch
import torch.nn as nn
import torch.optim as optim
from torchvision import datasets, transforms

# Check GPU
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Using device: {device}")

# Data loading
transform = transforms.Compose([
    transforms.ToTensor(),
    transforms.Normalize((0.1307,), (0.3081,))
])

train_dataset = datasets.MNIST('./data', train=True, download=True, transform=transform)
train_loader = torch.utils.data.DataLoader(train_dataset, batch_size=64, shuffle=True)

# Simple model
class Net(nn.Module):
    def __init__(self):
        super(Net, self).__init__()
        self.fc1 = nn.Linear(784, 128)
        self.fc2 = nn.Linear(128, 10)
    
    def forward(self, x):
        x = x.view(-1, 784)
        x = torch.relu(self.fc1(x))
        return self.fc2(x)

model = Net().to(device)
optimizer = optim.Adam(model.parameters())
criterion = nn.CrossEntropyLoss()

# Training loop
for epoch in range(5):
    for batch_idx, (data, target) in enumerate(train_loader):
        data, target = data.to(device), target.to(device)
        optimizer.zero_grad()
        output = model(data)
        loss = criterion(output, target)
        loss.backward()
        optimizer.step()
    print(f"Epoch {epoch+1}, Loss: {loss.item():.4f}")
```

## Multi-GPU Training

### Single Node, Multiple GPUs
```bash
#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:4
#SBATCH --cpus-per-task=32
#SBATCH --mem=256g

module load python/3.12 cuDNN/8.9.2/CUDA-12 CUDA/12.1
python -m torch.distributed.launch --nproc_per_node=4 train.py
```

### DataParallel (Simple)
```python
import torch.nn as nn

model = Net()
if torch.cuda.device_count() > 1:
    model = nn.DataParallel(model)
model.to(device)
```

### DistributedDataParallel (Recommended)
```python
import torch.distributed as dist
from torch.nn.parallel import DistributedDataParallel

dist.init_process_group(backend='nccl')
model = DistributedDataParallel(model)
```

## GPU Memory Management

```python
# Check GPU memory
print(torch.cuda.memory_allocated() / 1e9, "GB allocated")
print(torch.cuda.memory_reserved() / 1e9, "GB reserved")

# Clear cache
torch.cuda.empty_cache()

# Gradient checkpointing (reduce memory)
from torch.utils.checkpoint import checkpoint
```

## Mixed Precision Training

```python
from torch.cuda.amp import autocast, GradScaler

scaler = GradScaler()

for data, target in train_loader:
    optimizer.zero_grad()
    
    with autocast():
        output = model(data)
        loss = criterion(output, target)
    
    scaler.scale(loss).backward()
    scaler.step(optimizer)
    scaler.update()
```

## TensorBoard

```bash
# In training script
from torch.utils.tensorboard import SummaryWriter
writer = SummaryWriter('./logs')
writer.add_scalar('Loss/train', loss, epoch)

# Launch TensorBoard (with tunnel)
sinteractive --partition=gpu --gres=gpu:a100:1 --tunnel
module load python/3.12
tensorboard --logdir=./logs --port $PORT1 --bind_all

# Access via http://localhost:PORT1 after creating tunnel
```

## Monitoring GPU Usage

```bash
# On compute node
nvidia-smi              # Current status
nvidia-smi -l 1         # Update every second
watch -n 1 nvidia-smi   # Continuous monitoring

# GPU utilization in script
nvidia-smi --query-gpu=utilization.gpu,memory.used --format=csv -l 1
```

## Common Issues

### CUDA Out of Memory
- Reduce batch size
- Use gradient checkpointing
- Use mixed precision (fp16)
- Clear cache: `torch.cuda.empty_cache()`

### Module Compatibility
```bash
# Always load in this order:
module load python/3.12 cuDNN/8.9.2/CUDA-12 CUDA/12.1
```

### DataLoader Workers
```python
# Set num_workers based on allocated CPUs
import os
num_workers = int(os.environ.get('SLURM_CPUS_PER_TASK', 4))
loader = DataLoader(dataset, num_workers=num_workers, pin_memory=True)
```

### Reproducibility
```python
import torch
import numpy as np
import random

seed = 42
torch.manual_seed(seed)
torch.cuda.manual_seed_all(seed)
np.random.seed(seed)
random.seed(seed)
torch.backends.cudnn.deterministic = True
```

## Swarm of GPU Jobs

```bash
# gpu_train.swarm
python train.py --lr 0.001 --batch 32
python train.py --lr 0.01 --batch 32
python train.py --lr 0.001 --batch 64
python train.py --lr 0.01 --batch 64
```

```bash
swarm -f gpu_train.swarm --partition=gpu --gres=gpu:v100x:1 -g 32 -t 8 \
    --module python/3.12,cuDNN/8.9.2/CUDA-12,CUDA/12.1
```
