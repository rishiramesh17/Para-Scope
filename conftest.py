"""Make `src/` importable without depending on the editable install's .pth file.

This repo lives under ~/Desktop, which is iCloud-synced, and iCloud intermittently
sets the macOS UF_HIDDEN flag on files inside .venv. Python 3.13's
site.addpackage() skips hidden .pth files, so the editable install
(__editable__.para_scope-*.pth) stops being processed at random and every
`import para_scope` fails with ModuleNotFoundError -- including the two .pth
files setuptools itself relies on.

Clearing the flag works until iCloud re-applies it, so the durable fix is not to
need the .pth at all. A src-layout conftest that prepends src/ is standard
practice independent of this bug; here it also makes the suite reproducible on a
machine whose venv is being modified underneath it.

The real fix is to move this repo off the iCloud-synced Desktop -- see
docs/ROADMAP.md, "Environment".
"""

import sys
from pathlib import Path

_SRC = Path(__file__).parent / "src"
if _SRC.is_dir() and str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))
