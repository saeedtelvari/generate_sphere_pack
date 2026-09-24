# Issue tracker: Local Markdown

Issues and wayfinder maps for this repo live as markdown files in `.scratch/`.

## Conventions

- One effort per directory: `.scratch/<effort>/` (e.g. `.scratch/porous-media-generator/`)
- The map is `.scratch/<effort>/map.md`
- Tickets are individual markdown files at `.scratch/<effort>/issues/NN-<slug>.md`, numbered from `01`
- Triage/claim state is recorded as a `Status:` line (`open`, `claimed`, `resolved`)
- Ticket type is recorded as a `Type:` line (`research`, `prototype`, `grilling`, `task`)
- Dependencies are recorded as a `Blocked by: NN, NN` line
- Comments and resolutions append to the bottom of the file

## Wayfinding Operations

- **Map**: `.scratch/<effort>/map.md`
- **Child ticket**: `.scratch/<effort>/issues/NN-<slug>.md`
- **Blocking**: A ticket is unblocked when every ticket listed in `Blocked by:` has `Status: resolved`
- **Frontier**: Scan `.scratch/<effort>/issues/` for open, unblocked, unclaimed tickets; lowest number first
- **Claim**: Set `Status: claimed`
- **Resolve**: Append resolution under `## Answer`, set `Status: resolved`, add context pointer to `map.md`
