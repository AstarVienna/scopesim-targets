# Testing the MyST notebooks

`docs/conftest.py` collects the MyST text-notebooks in `docs/` and executes them
as pytest items. Plain markdown pages are ignored.

```bash
pytest docs --notebooks                      # all notebooks
pytest docs/quickstart.md --notebooks        # one notebook
pytest docs --notebooks -n auto --dist loadgroup   # parallel
```

## Front matter

### Detection (MyST-NB keys, also used by Sphinx)

A `.md` file is treated as a notebook if **either** of these is present:

```yaml
file_format: mystnb
```

```yaml
jupytext:
  text_representation:
    format_name: myst
```

Everything else in `docs/` is skipped, so no naming convention or ignore list is
needed.

### `kernelspec`

```yaml
kernelspec:
  name: python3
  display_name: Python 3 (ipykernel)   # optional
  language: python                     # optional
```

Only `name` is required. `display_name` and `language` are filled in
automatically on the in-memory copy when missing — MyST-NB tolerates a bare
`name`, but nbformat's schema does not.

### `myst_test` — test-only, invisible to Sphinx

```yaml
myst_test:
  needs: 01_setup.md        # or a list
  setup: false
```

**`needs`** — notebooks this one must run after, as paths relative to the
declaring file. Used when a chain hands off through the filesystem (one notebook
writes a file a later one reads). Effects:

- collected notebooks are topologically sorted;
- each connected chain is pinned to a single xdist worker via `xdist_group`, so
  `-n auto --dist loadgroup` is safe. Unrelated notebooks still run in parallel.

Error behaviour: a target that doesn't exist on disk is a collection error; a
cycle is a collection error; a target that exists but isn't part of the current
run (e.g. running one file directly) skips the dependent with a reason rather
than passing off a stale file.

Note this constraint exists in the Sphinx build too — MyST-NB gives every
document its own kernel, so a split chain can only hand off through the
filesystem there as well. Declaring `needs` documents an ordering that RTD
currently satisfies by accident of filename order.

**`setup: false`** — opt this notebook out of the injected preamble (below).

### `mystnb`

```yaml
mystnb:
  execution_mode: off
```

Skips execution. MyST-NB's own key, so it applies to the Sphinx build as well.
Notebooks with no code cells are skipped automatically.

## Cell tags

Executed through nbclient, so the usual tags work and behave as they do under
Sphinx:

- `raises-exception` — the cell is expected to fail
- `skip-execution` — the cell is not run

## Injected setup

If `docs/_notebook_setup.py` exists, its contents run in each notebook's kernel
before that notebook's own cells. This is for test-only environment wiring that
must not appear in the published docs — linking the IRDB, for example:

```python
from scopesim import link_irdb

link_irdb()
```

It is executed in the kernel, not inserted as a cell, so reported cell numbers
still match the `{code-cell}` blocks in the source and the notebook object is
unchanged. If it fails, the notebook fails immediately with
`notebook setup failed` and its traceback.

Override the path with `--notebook-setup=PATH` (useful when reusing this plugin
in another repo). Absent file, no injection. Keep it to environment wiring: the
more it stubs or patches, the less a passing test says about the page a reader
actually gets.

## Flags

| Flag                    | Effect                                                                        |
| ----------------------- | ----------------------------------------------------------------------------- |
| `--notebooks`           | Enable collection of notebooks. Without it, nothing in `docs/` is collected.  |
| `--notebook-setup=PATH` | Use `PATH` instead of `docs/_notebook_setup.py`.                              |
| `--notebook-keep-going` | Run a notebook's remaining cells after a cell fails and report every failure. |

## Defaults

- Kernel cwd is the notebook's own directory, matching the Sphinx build.
- 600 s per notebook, 120 s for kernel startup (`TIMEOUT` / `STARTUP_TIMEOUT` in
  `conftest.py`).
- Coverage from the kernel subprocess is collected via pytest-cov; run with
  `--cov=scopesim_targets`.
