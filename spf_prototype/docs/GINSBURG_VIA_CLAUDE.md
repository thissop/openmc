# Running on Ginsburg through Claude — SSH multiplexing (no password ever seen by Claude)

Standard instructions for any Claude Code session that needs to drive the Columbia **Ginsburg** cluster
(SLURM). The model: **you authenticate once (password + Duo) in your own terminal; Claude reuses that
already-authenticated connection over a shared local socket and never sees, handles, or transmits any
credential.** Project-agnostic — the paths/host below are examples; swap in your own.

## TL;DR
1. **You** run `ssh ginsburg` in a terminal once, complete the Duo push, and leave it (the connection
   persists ~8 h in the background even if you close the window).
2. **Claude** then runs `ssh ginsburg "<cmd>"` in its Bash tool. These reuse your live master
   connection — no new login, no password, no Duo. Claude cannot see the secret; it only borrows the
   authenticated channel.
3. Quick things (`squeue`, `ls`, `cat`, `sbatch`) run on the login node; **anything heavy runs as a
   SLURM batch job** (`sbatch`), never directly on the login node.

## Why this is safe (the security property)
- Authentication (password + Duo) happens **only in your terminal**, establishing an SSH
  **ControlMaster** — a persistent master connection with a local Unix socket at
  `~/.ssh/cm-<user>@<host>:<port>`.
- Claude's `ssh ginsburg …` commands run with **`ControlMaster auto`**, so they **multiplex over your
  existing socket** instead of authenticating. No credential is sent, prompted, or stored anywhere
  Claude can read. The socket is an *authenticated channel*, not a credential store — there is no
  password in it to extract.
- Claude always uses **`ssh -o BatchMode=yes`**, which disables all interactive prompts: if the master
  is down, Claude's ssh **fails immediately** ("Permission denied") rather than ever seeing a password
  prompt. So Claude structurally cannot capture your password.
- Net: **you hold the 2FA; Claude holds only a borrowed, already-open pipe.** Revoke it any time by
  closing the master (`ssh -O exit ginsburg`).

## One-time setup (in `~/.ssh/config`)
```sshconfig
Host ginsburg
  HostName ginsburg.rcs.columbia.edu
  User <your-uni>
  ControlMaster auto
  ControlPath ~/.ssh/cm-%r@%h:%p
  ControlPersist 8h            # master lives 8 h after the last use
  ServerAliveInterval 60       # keepalive so it doesn't idle out
```

## Each session
- **You:** `ssh ginsburg` → approve Duo → you can close the terminal (ControlPersist keeps it alive).
- **Check it's up** (Claude or you): `ssh -O check ginsburg` → `Master running (pid=…)`.
- If it dropped (Claude sees `Permission denied (publickey…)`), just `ssh ginsburg` + Duo again. The
  socket does expire (8 h, or on network changes / sleep); re-auth is the only manual step, and it is
  yours by design.

## How Claude should use it (patterns learned in practice)
- **Quick queries** (login node, seconds only): `ssh ginsburg "squeue -u <user>"`, `cat`, `ls`, `sbatch`.
- **Compute** → **always `sbatch`** a job script; never run solvers/builds/transport on the login node.
  Poll with a background Bash loop (`squeue -h -j <id> | wc -l`) that notifies on completion.
- **File transfer:** `rsync` over the multiplexed connection can flake ("unexpected end of file");
  **tar-over-ssh is robust**: `( cd src && tar czf - files ) | ssh ginsburg "cd dst && tar xzf -"`.
  On macOS add `COPYFILE_DISABLE=1` before `tar` to skip AppleDouble xattr noise.
- **Env gotchas on this cluster:** activate the conda env *and* `export PATH="$CONDA_PREFIX/bin:$PATH"`
  (conda activate alone may leave base Python first); run inline scripts from a neutral `cwd` (e.g.
  `/tmp`) so a repo dir literally named `openmc/` doesn't shadow the installed package.
- **Do NOT** ask for, cache, or type the password/Duo. If a command returns `Permission denied`, report
  it and ask the user to re-`ssh ginsburg`; do not attempt `ssh -MNf` or any interactive auth.

## What to tell a fresh Claude session
> "Ginsburg is reachable via `ssh ginsburg` using an SSH ControlMaster I (the user) authenticate. Use
> `ssh -o BatchMode=yes ginsburg "<cmd>"`; quick queries on the login node, all compute via `sbatch`;
> tar-over-ssh for files. You never handle my password/Duo — if you get `Permission denied`, tell me and
> I'll re-auth. See `spf_prototype/docs/GINSBURG_VIA_CLAUDE.md`."
