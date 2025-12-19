# Biowulf Hardware Reference

## Cluster Overview
- **Total**: ~2,911 nodes, ~93,292 cores, ~1,040 GPUs
- **Storage**: 35+ Petabytes across multiple systems

## CPU Nodes

| Nodes | Cores/Node | CPUs/Node | Memory | Network | SLURM Features |
|-------|------------|-----------|--------|---------|----------------|
| 64 | 96 (AMD Epyc 9454) | 192 | 768 GB | HDR200 IB | e9454,cpu192,g768 |
| 72 | 64 (AMD Epyc 7543) | 128 | 512 GB | HDR200 IB | e7543,cpu128,g512 |
| 243 | 36 (Intel Gold 6140) | 72 | 384 GB | HDR100 IB | x6140,cpu72,g384 |
| 1152 | 28 (Intel E5-2680v4) | 56 | 256 GB | FDR IB | x2680,cpu56,g256 |
| 1080 | 28 (Intel E5-2695v3) | 56 | 256 GB | FDR IB | x2695,cpu56,g256 |

## GPU Nodes

| Nodes | Processor | GPUs | GPU Memory | Node Memory | Features |
|-------|-----------|------|------------|-------------|----------|
| 76 | AMD Epyc 7543p (32 cores) | 4x A100 | 80 GB each | 256 GB | gpua100,ibhdr200 |
| 56 | Intel Gold 6140 (36 cores) | 4x V100-SXM2 | 32 GB each | 384 GB | gpuv100x,ibhdr |
| 8 | Intel E5-2680v4 (28 cores) | 4x V100 | 16 GB each | 128 GB | gpuv100 |
| 48 | Intel E5-2680v4 (28 cores) | 4x P100 | 16 GB each | 128 GB | gpup100 |
| 72 | Intel E5-2680v4 (28 cores) | 4x K80 | 24 GB each | 256 GB | gpuk80 |

## Large Memory Nodes

| Nodes | Cores | CPUs | Memory | Features |
|-------|-------|------|--------|----------|
| 16 | 96 (AMD Epyc 9454) | 192 | 3 TB | g3072,ibhdr200 |
| 4 | 72 (Intel E7-8860v4) | 144 | 3 TB | g3072 |
| 20 | 72 (Intel E7-8860v4) | 144 | 1.5 TB | g1536 |

## GPU Allocation Examples

```bash
# Request A100 GPU
sbatch --partition=gpu --gres=gpu:a100:1 script.sh

# Request V100x GPU (32GB VRAM)
sbatch --partition=gpu --gres=gpu:v100x:2 script.sh

# Request P100 GPU
sbatch --partition=gpu --gres=gpu:p100:1 script.sh

# Interactive GPU session
sinteractive --partition=gpu --gres=gpu:a100:1 --mem=64g --cpus-per-task=8
```

## Checking Available Resources

```bash
# Show free nodes and GPUs
freen

# Show GPU availability specifically
freen | grep -E 'Partition|----|gpu'

# Output example:
# Partition      FreeNds      FreeCPUs      FreeGPUs  Cores  CPUs  GPUs    Mem
# gpu (a100)      5 / 76     1200 / 4864    20 / 304     32    64     4   249g
# gpu (v100x)     2 / 53     2042 / 3816    17 / 212     36    72     4   373g
```

## Storage Systems

| System | Type | Capacity | Best For |
|--------|------|----------|----------|
| VAST (NVMe Flash) | NFS over RDMA | 35 PB | High-performance I/O |
| DDN SFA18K | GPFS | 7.6 PB | Large sequential I/O |
| DDN SFA12KX | GPFS | 21.6 PB | General storage |
| NetApp | NFS | 450 TB | Home directories |

## Network
- Nodes connected via 56-200 Gb/s Infiniband
- GPU nodes have 1:1 non-blocking connection to core
- Direct 100 Gb/s Ethernet to NIHnet
