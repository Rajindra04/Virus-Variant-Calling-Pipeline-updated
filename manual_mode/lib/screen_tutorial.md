# GNU Screen Mini-Tutorial

Several pipeline steps take 10-30 minutes per sample. If your SSH connection drops, the job dies. **Screen** keeps your session alive on the server.

## Essential Commands

| Action | Command |
|--------|---------|
| Start a new session | `screen -S pipeline` |
| Detach (leave running) | `Ctrl+A`, then `D` |
| Reattach | `screen -r pipeline` |
| List sessions | `screen -ls` |
| Kill current session | `Ctrl+A`, then `K`, then `Y` |

## Workflow Example

```bash
# Start a screen session
screen -S denv2_mapping

# Run a long command
bash 04_mapping.sh SRR35818859

# Detach: press Ctrl+A, then D
# You'll see: [detached from 12345.denv2_mapping]

# Safe to close your laptop / terminal now

# Later, reconnect:
ssh htcf
screen -r denv2_mapping
# You're back where you left off
```

## Running Multiple Samples in Parallel

```bash
# Start separate sessions for each sample
screen -S sample1
bash 04_mapping.sh SRR35818859
# Ctrl+A, D to detach

screen -S sample2
bash 04_mapping.sh SRR35818860
# Ctrl+A, D to detach

# Check on them
screen -ls
# There are screens on:
#   12345.sample1   (Detached)
#   12346.sample2   (Detached)

# Reattach to check progress
screen -r sample1
```

## If Screen Says "Already Attached"

```bash
# Force detach the other connection, then reattach
screen -d -r pipeline
```

## Alternative: tmux

If `screen` is not available, `tmux` works similarly:

| Screen | tmux |
|--------|------|
| `screen -S name` | `tmux new -s name` |
| `Ctrl+A, D` | `Ctrl+B, D` |
| `screen -r name` | `tmux attach -t name` |
| `screen -ls` | `tmux ls` |
