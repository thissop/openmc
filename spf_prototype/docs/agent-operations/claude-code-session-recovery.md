# Claude Code session monitoring and context-aware recovery

## Purpose

This is a local macOS operating procedure for finding independent Claude Code sessions, reading enough of each session's recent state to understand what stopped, and safely sending a continuation prompt after usage becomes available again.

This is not a Claude-native cross-session messaging feature. It is host-level automation: inspect local processes and transcripts, identify the correct macOS Terminal tab, and type into that tab as though a user had done so.

The procedure is intentionally generic. Never preserve a PID, TTY, reset time, or transcript filename as permanent truth; resolve all of them again immediately before acting.

## Preconditions

- The Claude Code sessions are running locally in macOS Terminal.
- The operator has permission to inspect those processes and repositories.
- Terminal automation is permitted. The visible paste/Return fallback also requires macOS Accessibility permission.
- Claude's local transcripts are available under ~/.claude/projects.
- A scheduler such as a Codex heartbeat is available if recovery must happen after the present conversation ends.

If any of these assumptions is false, stop and report the limitation rather than guessing.

## Safety rules

1. Match a session by repository working directory and TTY, not by window order, tab title, or an old PID alone.
2. Read the newest state before sending anything. A session may have recovered, completed, compacted, or changed tasks since the monitor was scheduled.
3. Resume only an explicit usage- or session-limit block. Do not confuse the following with a limit:
   - active thinking or tool execution;
   - normal context compaction;
   - a completed summary followed by an ordinary prompt;
   - waiting for user approval or credentials;
   - a healthy long-running job or monitor.
4. Never broadcast a prompt to every Terminal tab.
5. Never start a second writer in a repository that is already being modified by an active session.
6. A successful automation call is not proof that Claude received the message. Verify new Terminal or transcript activity.
7. Keep prompts scoped to the work already authorized by the user.
8. Avoid copying secrets, tokens, private keys, or unrelated transcript content into restart messages or logs.

## 1. Discover live Claude processes

Start with a read-only inventory:

~~~sh
ps -axo pid=,ppid=,tty=,command= \
  | rg -i '(^|[ /])claude([[:space:]]|$)|Claude Code' \
  | rg -v 'rg -i|Codex'
~~~

For each plausible interactive Claude PID, resolve its working directory:

~~~sh
task_pid=12345
ps -p "$task_pid" -o pid=,ppid=,tty=,stat=,command=
lsof -a -p "$task_pid" -d cwd -Fn
~~~

The lsof result is normally shaped like:

~~~text
p12345
fcwd
n/Users/name/Documents/GitHub/example-repo
~~~

Record a temporary mapping:

~~~text
repository root -> current PID -> /dev/ttysNNN
~~~

Treat daemon and background-PTY processes separately from interactive Terminal sessions. A daemon's origin metadata can be a clue, but it is not proof that an interactive tab is still present.

## 2. Confirm that Claude is the foreground program

Before typing into a TTY, inspect all processes attached to it:

~~~sh
tty_name=ttys001
ps -t "$tty_name" -o pid=,ppid=,pgid=,tpgid=,stat=,command=
~~~

The foreground process group has PGID equal to TPGID. Confirm that the foreground program is the intended Claude TUI. If a shell, editor, SSH password prompt, or some other program is in front, do not inject a Claude prompt.

## 3. Locate the freshest transcript

Claude commonly stores project transcripts below a path-encoded project directory:

~~~text
~/.claude/projects/-Users-name-Documents-GitHub-example-repo/*.jsonl
~~~

List the newest candidates rather than assuming the transcript UUID:

~~~sh
project_dir="$HOME/.claude/projects/-Users-name-Documents-GitHub-example-repo"

find "$project_dir" -maxdepth 1 -type f -name '*.jsonl' \
  -exec stat -f '%m %Sm %N' -t '%Y-%m-%d %H:%M:%S' {} \; \
  | sort -nr \
  | head
~~~

Forked or resumed sessions can create a newer JSONL than the UUID shown on the Claude command line. Compare modification times and recent message ancestry. Use the file containing the current TUI's newest messages.

The JSONL format is an implementation detail and may change. Inspect defensively. A useful current-style extraction is:

~~~sh
transcript=/absolute/path/to/current-session.jsonl

tail -n 300 "$transcript" \
  | jq -r '
      select(.type == "user" or .type == "assistant")
      | .timestamp as $timestamp
      | .message.role as $role
      | (
          .message.content
          | if type == "array" then
              map(select(.type == "text") | .text) | join(" ")
            elif type == "string" then
              .
            else
              ""
            end
        ) as $text
      | select($text != "")
      | "\($timestamp)\t\($role)\t\($text)"
    '
~~~

Search for explicit blockers while retaining enough surrounding context to understand the unfinished work:

~~~sh
rg -n -i \
  "hit your (session|usage|monthly spend) limit|resets|usage-credits|rate limit" \
  "$transcript" \
  | tail -20
~~~

Also inspect the latest user request, assistant plan, todo state, background-task notifications, and repository status. Do not reduce a complex session to only the limit message.

## 4. Classify the current state

Use evidence from both the Terminal tab and the JSONL tail.

| State | Evidence | Action |
|---|---|---|
| Limit-blocked | Latest response explicitly reports a usage/session limit; Claude is idle at that result | Eligible for a context-aware restart after the stated reset |
| Active | Transcript is growing, tools are running, or the TUI shows thinking/work | Leave it alone |
| Normally idle | Claude completed a response and awaits ordinary user input | Do not call this frozen |
| Compacting | The session is summarizing or restoring context without an explicit limit | Let compaction finish |
| Waiting on external input | SSH, MFA, approval, credentials, or a user decision is required | Report the real blocker; do not spam Continue |
| Completed | The requested work is finished or a final handoff was delivered | Do not restart |
| Ambiguous | Terminal and transcript disagree, or the foreground process cannot be proven | Do not inject; investigate or report |

A restart decision should be conservative. False negatives cost a later manual nudge; false positives can interrupt work or create competing edits.

## 5. Build the continuation prompt

Construct the prompt from the target session's freshest evidence, not from an old generic summary.

Include:

- the immediate task that was underway before the limit;
- the last unanswered user request;
- relevant todo or background-task state;
- any repository, branch, job, workflow, or run identifier needed to resume;
- a request to inspect current state before changing anything;
- an instruction to continue within the user's existing scope;
- a reminder not to duplicate work already completed while the session was paused.

Template:

~~~text
Your usage is available again. Continue the existing task from the point where
the session limit stopped you. First reread your newest transcript/todo state
and inspect the current repository and running-job state so you do not repeat
completed work.

Immediate stopping point:
- [fresh, concrete stopping point]

Outstanding user request:
- [last unanswered request]

Relevant live state:
- [branch, job, workflow, run ID, monitor, or file]

Resume the authorized work, preserve existing changes, verify that work actually
continues, and report any new blocker honestly.
~~~

Keep session-specific facts in the message for that session only.

## 6. Submit to one exact Terminal tab

The proven macOS Terminal route is:

1. Find the tab whose tty exactly matches the refreshed TTY.
2. Insert the prompt with Terminal's do script command.
3. Send an empty do script to the same tab; in the observed Claude TUI this supplies the Return needed to submit.
4. Confirm that the text was submitted rather than merely left in the input box.

Example:

~~~sh
target_tty=/dev/ttys001
prompt_text='Your usage is available again. Reread your newest state and continue the authorized work from the exact stopping point.'

osascript - "$target_tty" "$prompt_text" <<'APPLESCRIPT'
on run argv
    set targetTTY to item 1 of argv
    set promptText to item 2 of argv

    tell application "Terminal"
        repeat with terminalWindow in windows
            repeat with terminalTab in tabs of terminalWindow
                try
                    if (tty of terminalTab) is targetTTY then
                        do script promptText in terminalTab
                        delay 0.2
                        do script "" in terminalTab
                        return "submitted to " & targetTTY
                    end if
                end try
            end repeat
        end repeat
    end tell

    error "No Terminal tab matched " & targetTTY
end run
APPLESCRIPT
~~~

Important: Terminal do script is not a Claude API. If Claude is not the foreground TUI, the text may be interpreted by a shell or another program. That is why the foreground-process and visible-history checks are mandatory.

If the text appears in Claude's input box but does not submit, use a visible paste/Return fallback only when Accessibility permission is enabled. Preserve and restore the user's clipboard if that fallback uses it.

## 7. Verify recovery

Capture the pre-send transcript modification time:

~~~sh
before_mtime=$(stat -f '%m' "$transcript")
~~~

After submission, wait briefly and inspect:

~~~sh
after_mtime=$(stat -f '%m' "$transcript")
tail -n 30 "$transcript"
ps -p "$task_pid" -o pid=,tty=,stat=,command=
~~~

Recovery is verified when at least one strong signal appears:

- new assistant/tool records are appended after the restart prompt;
- the Terminal tab visibly enters thinking or tool execution;
- expected background jobs or workflow records begin updating.

If Claude immediately returns another limit message, record that the reset has not propagated. Do not send repeated prompts in a tight loop.

## 8. Scheduling a reset-time recovery

For a future reset, use a durable scheduler such as a Codex heartbeat rather than leaving a shell sleep running.

The scheduled instructions must be self-contained because PIDs, TTYs, transcripts, and work state can change before the wake-up. Include:

- exact repository roots, with current PIDs/TTYs only as hints;
- the explicit reset time shown by Claude and its timezone;
- a small post-reset margin;
- the full discover -> inspect -> classify -> prompt -> submit -> verify sequence;
- context hints for each session, clearly labeled as potentially stale;
- a prohibition on touching sessions that already resumed;
- a requirement to delete or disable the one-time heartbeat after a terminal check.

When the scheduler wakes, it must re-read the latest state before acting. Never let the scheduled prompt blindly send a prewritten Continue message.

## 9. Common failure modes

### AppleScript reports success but Claude does not respond

The first do script may have placed text in the input box without submitting it. Inspect the tab and send the explicit empty do script/Return step. Then verify JSONL growth.

### No Terminal tab matches the old TTY

The session restarted, moved, or closed. Re-run process and cwd discovery. Do not reuse the old mapping.

### Several Claude processes share a repository

Use TTY, foreground process group, command line, transcript recency, and visible context together. If the intended writer remains ambiguous, do not inject.

### The newest transcript is not the resumed UUID

A fork or resume can create a new transcript. Use timestamps and message continuity rather than the original UUID alone.

### A background-task failure contains a limit message

Determine whether the parent session itself is blocked. A failed subagent notification does not automatically mean the main TUI needs another prompt.

### The reset time passed but usage is still blocked

Allow a short propagation margin and check again once. Repeated automatic submission can duplicate prompts when service state catches up.

### The session needs credentials or approval

A continuation prompt cannot solve MFA, SSH authentication, or a missing user decision. Report the actual blocker.

## 10. Completion checklist

- [ ] Live Claude processes were rediscovered.
- [ ] Each target was matched by exact cwd and refreshed TTY.
- [ ] Claude was confirmed as the foreground program.
- [ ] The current transcript was selected by recency and continuity.
- [ ] The latest state was classified correctly.
- [ ] Only explicit limit-blocked sessions were touched.
- [ ] Each restart prompt used fresh, session-specific context.
- [ ] Submission targeted one exact Terminal tab.
- [ ] Recovery was verified from real activity.
- [ ] The user received a concise result report.
- [ ] Any one-time scheduler was removed or disabled.

## Provenance and maintenance

This playbook records a locally demonstrated workaround used across several research repositories. It should be updated when Claude's process layout, transcript format, Terminal behavior, or the scheduler changes.

Do not describe the method as built-in Claude cross-session communication. The durable principle is broader than any one command: resolve identity from the host, reconstruct fresh context, intervene conservatively, and verify the outcome.
