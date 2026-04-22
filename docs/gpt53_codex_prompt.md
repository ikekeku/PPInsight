# GPT-5.3 Codex Handoff Prompt

Use this prompt as-is with GPT-5.3 Codex when you want it to perform the same class of work in this repo.

---

You are working in the `PPInsight` repo.

Goals:

1. Inspect the current code/docs first. Do not assume previous changes are present.
2. Make the fetch-to-pairs handoff user-friendly:
   - Fetched structures should be saved and reusable with user-friendly accession-based names.
   - Users should be able to put those same identifiers directly into `proteinA` / `proteinB`.
   - Preserve backward compatibility with legacy `.ent` naming if older files still exist.
   - Add or update tests for the naming/path-resolution behavior.
3. Do a repo-wide sweep for stale wording and workflow drift:
   - Prefer umbrella CLI wording: `ppinsight <subcommand>`.
   - Make docs/help/tutorials match the real workflow, especially `fetch`, `batch`, `collect`, and `compare`.
   - Remove stale references to retired or unsupported compare plots.
   - Keep wording accessible for users with little or no CLI familiarity.
   - Prefer user-facing names over internal/Biopython-style filenames in generic docs, but keep concrete bundled example identifiers where they are the truthful names of shipped repo assets.
   - In tutorials, distinguish clearly between verified shipped-example commands and reference commands that require a real aligned benchmark or extra runtime dependencies.
4. Treat maintainer-context docs as approval-gated:
   - Do not edit `docs/key_decisions_log.md` or `docs/non_negotiables.md` unless the maintainer explicitly asks for it in the current task.
   - If approval is given, read the current contents first.
   - Update them additively; do not replace, rewrite wholesale, or delete existing content.
   - Keep additions dense and bullet-oriented.
   - In `docs/key_decisions_log.md`, use the structure: `Decision made`, `Alternatives considered`, `Why rejected`, `Tradeoffs accepted`.
   - In `docs/non_negotiables.md`, keep additions aligned with the existing sections rather than inventing a new structure.
5. Validate in the real project environment:
   - Use the `ppinsight` conda environment, not system Python.
   - If Ruff is missing, install it in that environment.
   - Run Ruff on the whole repo.
   - Run targeted `python -m pytest` on the files/areas you changed.
   - If you touch runtime code, do not stop until the touched area is lint-clean and the targeted tests pass.
6. Report back with:
   - What changed.
   - Which files were edited.
   - What validation commands were run and whether they passed.
   - Any remaining uncertainties that require maintainer check-in.

Guardrails:

- Do not do unrelated refactors.
- Do not restructure, slim down, or delete repo content unless explicitly asked.
- Do not delete, replace, or edit the maintainer docs without explicit approval for that task.
- Keep changes minimal, truthful, and consistent with the existing repo style.
- Fix root-cause UX/documentation mismatches when possible instead of layering more explanation onto a confusing workflow.
- After behavior changes, sweep nearby docs/help/tutorial text in the same pass.
- Preserve the product contract that supported workflow interfaces ship with PPInsight, while exceptionally heavy, separately licensed, or post-hoc runtimes may remain install-when-needed with lazy imports and explicit install guidance.
- Use `python -m pytest`, not bare `pytest`.

Validation commands to include if needed:

```bash
conda run -p /Users/okik/miniconda3/envs/ppinsight python -m pip install ruff
conda run -p /Users/okik/miniconda3/envs/ppinsight ruff check .
```

If the environment path changes, use the active `ppinsight` environment equivalent instead.