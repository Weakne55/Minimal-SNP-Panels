# Recommended GitHub configuration

## 1. Repository visibility

The current repository should preferably remain **Private** until:
- ownership of the code is clear;
- publication/release is approved;
- no unpublished metadata or institute-specific information is accidentally exposed.

Do not publish the 12 in-house polar bear samples or raw data in GitHub.

## 2. GitHub Project fields

Keep GitHub Issues as the work items and use Project fields for planning.

### Status
- Backlog
- Ready
- In progress
- Blocked
- Review
- Done

### Phase
- Literature
- Data acquisition
- WGS processing
- Variant QC
- Population structure
- SNP selection
- Validation
- Benchmarking
- Reporting

### Species
- Beluga
- Polar bear
- Both / generic

### Priority
- P0 — critical path
- P1 — high
- P2 — normal
- P3 — optional

### Target date
Date field.

### Iteration
Weekly iteration field. Use one-week iterations to match work reporting.

## 3. Recommended Project views

### Board — Execution
Layout: Board  
Group by: Status  
Use for daily work.

### Table — Master plan
Layout: Table  
Columns: Title, Status, Phase, Species, Priority, Target date, Iteration  
Group by: Phase.

### Roadmap — Deadline
Layout: Roadmap  
Date field: Target date  
Group by: Phase.

### Table — This week
Layout: Table  
Filter to current iteration and open items.

### Table — Literature
Filter to Phase = Literature.

## 4. Built-in Project workflows

Enable:
- item added to project → Status = Backlog or Ready;
- issue/PR closed → Status = Done;
- merged PR → Status = Done;
- auto-add issues from the `Minimal-SNP-Panels` repository.

Optional later:
- auto-archive old completed items.

## 5. Repository labels

Use labels only for information that should also be useful outside GitHub Projects.

Recommended:
- `pipeline`
- `analysis`
- `literature`
- `documentation`
- `bug`
- `blocked`
- `decision-needed`

Do not duplicate Species/Priority/Phase as labels if they already exist as Project fields.

## 6. Milestone

Create one repository milestone:

**v0.1 — project delivery**  
Due: **2026-12-25**

Attach all mandatory delivery issues to this milestone.

## 7. Main branch protection / ruleset

For a one-person project, keep this lightweight:
- block force pushes to `main`;
- block deletion of `main`;
- optionally require passing CI checks once CI exists.

Do not require PR approval from another person unless there is actually a reviewer.

## 8. Security and dependencies

Useful once dependency files exist:
- Dependency graph;
- Dependabot alerts;
- secret scanning where available.

Never commit tokens, passwords, SSH keys or credentials. Keep local secrets in `.env` or cluster-side configuration that is ignored by Git.

## 9. Actions / CI

Add CI only after the first Bash/Python modules exist.

Recommended checks:
- `shellcheck` for Bash;
- `ruff` for Python;
- `pytest` for tests;
- a small test dataset / toy VCF for smoke tests.

Do not run real WGS processing in GitHub Actions.

## 10. License

Do **not** add an open-source license automatically if ownership of the work/code has not been agreed with the employer/institute.

A public repository without an explicit license is visible but does not grant normal reuse rights. Decide licensing before public release.
