# Backlog

The daily automation (`scripts/daily.py`) consumes this list top-to-bottom:
it takes the first undone item, implements it with the DeepSeek provider
(`scripts/llm.py`), and ships it only if the full test suite stays green.
Replenishment is handled by `scripts/replenish.py` (propose) and
`scripts/promote.py` (approve).

Keep entries concrete: what changes, and why it matters for lab / GMP work.

## Planned

- [ ] `cfu_report_export` Add a formatted human-readable report (not just JSON) to the `cfu` tool. | handler:
- [ ] `cfu_edge_cases` Harden `cfu` input validation (zero/negative volumes, empty dilutions) with tests. | handler:
- [ ] `cfu_docs` Add a worked example with real numbers to the `cfu` README section. | handler:
