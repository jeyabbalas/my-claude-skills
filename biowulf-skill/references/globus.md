# Globus Data Transfer Reference

Globus provides reliable, high-performance file transfers to/from Biowulf.

## NIH HPC Globus Collections

| Collection Name | UUID | Purpose |
|-----------------|------|---------|
| **NIH HPC Data Transfer (Biowulf)** | e2620047-6d04-11e5-ba46-22000b92c6ec | Main collection for /home, /data, shared directories |
| NIH HPC Internet2 - AWS S3 | c24547a8-ef53-4b86-bcf6-d050c55d00f4 | Transfer to/from AWS S3 |
| NIH HPC Google Cloud Collection | 46312f97-8565-456d-a1ea-cb3e28e49caa | Transfer to/from Google Cloud |
| NIH HPC Google Drive Collection | 22629017-758c-469c-8f75-51eaebcf0417 | Transfer to/from Google Drive |
| NIH HPC OneDrive Collection | ba595c6f-8822-4905-8ce9-6e072bb49ce4 | Transfer to/from NIH OneDrive |

## Getting Started

### 1. Log into Globus
1. Go to https://www.globus.org
2. Click "Log In"
3. Select "National Institutes of Health"
4. Authenticate with PIV card or NIH credentials

### 2. Install Globus Connect Personal (for your desktop)
Download from: https://www.globus.org/globus-connect-personal

**Important**: Install while OFF the VPN

#### Installation Notes
- **Windows**: May need admin rights; change install directory if not admin
- **Mac/Linux**: Follow standard installation
- **Don't select "High Assurance"** during setup (NIH doesn't have this subscription)

### 3. Configure Local Access
By default, Globus Connect Personal only sees your home directory.
To access network drives or other folders:
- **Windows**: https://docs.globus.org/how-to/globus-connect-personal-windows/#configuration
- **Mac**: https://docs.globus.org/how-to/globus-connect-personal-mac/#configuration
- **Linux**: https://docs.globus.org/faq/globus-connect-endpoints/#how_do_i_configure_accessible_directories_on_globus_connect_personal_for_linux

## Transferring Files

### Web Interface
1. Log into https://app.globus.org
2. In File Manager, set:
   - **Left panel**: Source collection (e.g., your desktop)
   - **Right panel**: Destination (e.g., "NIH HPC Data Transfer (Biowulf)")
3. Navigate to folders, select files
4. Click "Start" to transfer

### Transfer to Biowulf
1. Source: Your Globus Connect Personal endpoint
2. Destination: "NIH HPC Data Transfer (Biowulf)"
3. Navigate to `/data/$USER/` on Biowulf side
4. Select files/folders → Start

### Transfer from Biowulf
1. Source: "NIH HPC Data Transfer (Biowulf)"
2. Destination: Your Globus Connect Personal endpoint
3. Select files from Biowulf → Start

## Command Line Interface

### Install CLI
```bash
pip install globus-cli
globus login  # Authenticate via browser
```

### List Endpoints
```bash
globus endpoint search "NIH HPC"
```

### Transfer Files
```bash
# Get endpoint IDs
BIOWULF=e2620047-6d04-11e5-ba46-22000b92c6ec
LOCAL=your-endpoint-id

# Single file
globus transfer $LOCAL:/path/to/file.txt $BIOWULF:/data/$USER/file.txt

# Directory (recursive)
globus transfer --recursive $LOCAL:/path/to/dir $BIOWULF:/data/$USER/dir

# Check transfer status
globus task list
globus task show TASK_ID
```

## Sharing Data

### Create a Shared Endpoint
1. In Globus web interface, navigate to folder to share
2. Click "Share" in right panel
3. Give the share a name
4. Add users/groups with their email or Globus identity
5. Set permissions (read-only or read-write)

### Access Shared Data
- Recipient logs into Globus
- Finds shared endpoint in "Shared With You"
- Transfers files normally

## Cloud Transfers

### AWS S3
Collection: "NIH HPC Internet2 - AWS S3"
Requires: S3 credentials configured in Globus

### Google Cloud Storage
Collection: "NIH HPC Google Cloud Collection"
Requires: Google Cloud credentials

### Google Drive
Collection: "NIH HPC Google Drive Collection"
Authenticate with Google account when prompted

### NIH OneDrive
Collection: "NIH HPC OneDrive Collection"
Authenticate with NIH Microsoft 365 account

## Scheduled/Recurring Transfers

Set up through Globus web interface:
1. Start a transfer
2. Before clicking Start, click "Transfer & Timer Options"
3. Configure schedule (daily, weekly, etc.)
4. Set start date/time

## Best Practices

1. **Use for large transfers** - Globus handles retries, checksums
2. **Use sync option** - Only transfers changed files
3. **Check quotas** - `checkquota` on Biowulf before large transfers
4. **No PII/PHI** - Unless special arrangements made
5. **Verify transfers** - Check file counts and sizes after transfer

## Troubleshooting

### "Browser login did not complete"
Log off VPN and retry installation

### Can't see network drive
Configure accessible paths in Globus Connect Personal settings

### Transfer failed
- Check disk quota on destination
- Verify permissions
- Check Globus task details for specific error

### Slow transfers
- Use "NIH HPC Data Transfer (Biowulf)" for transfers within NIH
- Use Internet2 endpoint only for transfers to other Internet2 sites

## Security Notes

- Transfers are encrypted in transit
- Authentication via NIH SSO
- Don't transfer PII/PHI without proper authorization
- See: https://hpc.nih.gov/policies/index.html#PII
