# Changelog

All notable changes to pickagm are listed here. Version numbers follow
<https://semver.org/>; while the version starts with `0.`, any release may
change the public interface.

## 0.1.0

First public release.

- The code in `src/pickagm/` is identical to commit
  [`9175fe6`](https://github.com/nc-hsu/pickagm/commit/9175fe6c825feb42087765537bc4b377a1ef5f43)
  (tag `record-selection-run`), the commit that the record selection for
  the accompanying paper ran on. Check it with:

  ```console
  git diff --exit-code record-selection-run v0.1.0 -- src/pickagm
  ```

- Packaging for PyPI: declared runtime dependencies, the bundled ESM
  flatfile and correlation coefficients installed as package data, an SPDX
  licence expression (`GPL-3.0-or-later`).
- Documentation on Read the Docs: <https://pickagm.readthedocs.io>.
- A smoke test that checks an installed wheel (`tests/test_installed_package.py`).
