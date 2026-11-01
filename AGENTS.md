# HumanProof

Read `README.md` first, then `CONTRIBUTING.md` for the architecture and the conventions.

## This repository is public

Everything committed here can be read by anyone.
Add nothing that is not meant to be public: no secret or key, and no note, result, data or code that comes from private work of the owner.
Another project of the owner reads the `data` folder of this one by link and checks its files against hashes, so a file there is not moved, renamed or rewritten without a note to the owner.

## How to work

- The local gate is authoritative, and red or missing CI is not evidence of breakage:
  `uv run --with pandas --with numpy --with pyarrow --with scikit-learn --with pytest python -m pytest tests -q -p no:cacheprovider`
- Stage files by name and never with `git add -A`.
- Commit locally.
  Push or tag only on an explicit go-ahead that names the push.
- Never read or print a secret value, and never put one in a command, a URL or a file.
- Never hand-edit a generated file: fix the script that makes it and run it again.
- Write no em or en dashes, and put one sentence on each line in long Markdown.
