# Contributing to PyNetworKIN

This project is a modular, research-grade pipeline for kinase–substrate prediction.
Contributions should maintain correctness, reproducibility, and clarity.

---

## Development Setup

We use `uv` for dependency and environment management.

### Clone and install (dev mode)

```bash
git clone https://github.com/bibymaths/pynetworkin
cd pynetworkin
uv sync --group dev
````

### Run CLI locally

```bash
uv run pynetworkin --help
```

### Run tests

```bash
uv run pytest
```

---

## Project Structure

Core code lives in:

```
src/pynetworkin/
```

Key modules:

* `networkin.py` → main pipeline
* `motif_scoring.py` → NetPhoREST wrapper
* `graph_scoring.py` → STRING-based scoring
* `recovery.py` → false-negative recovery
* `cli.py` → CLI entrypoint

Tests are in:

```
tests/
```

---

## Code Style

* Python ≥ 3.10
* Follow **PEP8**, but prioritize clarity over strict formatting
* Use type hints where meaningful
* Prefer small, composable functions over monolithic logic
* Avoid hidden side effects in pipeline components

Logging:

* Use the existing logging utilities (`logger.py`)
* Do not introduce ad-hoc print statements

---

## Branching Strategy

* `main` → stable, releasable code
* feature branches → `feature/<short-description>`
* bug fixes → `fix/<short-description>`

Examples:

```
feature/add-yeast-support
fix/string-cache-bug
```

---

## Pull Request Guidelines

Before submitting:

* [ ] Code runs (`uv run pynetworkin predict ...`)
* [ ] Tests pass (`uv run pytest`)
* [ ] New functionality includes tests
* [ ] No breaking changes without justification
* [ ] Documentation updated if behavior changes

PRs should be:

* Focused (one feature or fix)
* Technically justified
* Minimal in scope

---

## Data & Reproducibility

* Do not commit large datasets
* Use existing caching mechanisms
* Respect external data source constraints (OmniPath, STRING)

---

## Questions / Discussions

Open an issue for:

* API design changes
* Algorithmic modifications
* Data integration decisions

````

---

# 2. `CODE_OF_CONDUCT.md`

Keep this standard, not over-engineered.

```markdown
# Code of Conduct

## Scope

This project is an open, research-oriented software effort.
All participants are expected to maintain a professional and respectful environment.

---

## Expected Behavior

- Be direct, but respectful
- Focus discussions on technical merit
- Provide constructive feedback
- Respect differing levels of experience

---

## Unacceptable Behavior

- Personal attacks or harassment
- Dismissive or non-constructive criticism
- Discriminatory language or behavior

---

## Enforcement

Maintainers reserve the right to:

- Moderate discussions
- Reject contributions that do not meet standards
- Remove participants who violate this code

---

## Reporting

If issues arise, open an issue or contact maintainers directly.
