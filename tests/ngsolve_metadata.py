"""scikit-build-core dynamic-metadata provider for the ngsolve wheel.

Computes the PEP 440 package version from the git history, reusing the logic in
``tests/utils.py`` (the single source of truth, also used by the CI ``--check-pip``
and ``--wait-pip`` helpers).

Referenced from ``pyproject.toml``::

    [tool.scikit-build.metadata.version]
    provider = "ngsolve_metadata"
    provider-path = "tests"

scikit-build-core inserts ``provider-path`` (``tests``) onto ``sys.path`` and
imports this module, so ``import utils`` resolves to ``tests/utils.py``.
"""

from __future__ import annotations

import os
import sys

# scikit-build-core only keeps provider-path on sys.path while importing this
# module, so resolve and import ``utils`` (tests/utils.py) now, at import time.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils  # noqa: E402

__all__ = ["dynamic_metadata", "get_requires_for_dynamic_metadata"]


def __dir__() -> list[str]:
    return __all__


# BLAS/LAPACK is shipped as a separate wheel on Linux/Windows; macOS uses the
# system Accelerate framework.
_OPENBLAS_REQUIREMENT = "ngsolve-openblas==0.3.33; platform_system != 'Darwin'"


def _netgen_requirement() -> str:
    """Pin netgen-mesher to the exact version ngsolve is built against.

    The compiled ngsolve libraries link against netgen's libngcore etc., so the
    installed netgen-mesher must be ABI-compatible. The build environment has the
    matching netgen-mesher installed (see the CI build scripts / cibuildwheel
    before-build step), so read its version from there -- mirroring what the old
    setup.py did with ``netgen.config.NETGEN_VERSION_PYTHON``.
    """
    try:
        import netgen.config

        return f"netgen-mesher=={netgen.config.NETGEN_VERSION_PYTHON}"
    except Exception:
        # netgen not importable (e.g. plain sdist build): leave it unpinned.
        return "netgen-mesher"


def dynamic_metadata(field: str, settings: dict | None = None):
    if settings:
        msg = "No inline configuration is supported"
        raise ValueError(msg)

    if field == "version":
        return utils.get_version(os.getcwd())
    if field == "dependencies":
        return [
            _netgen_requirement(),
            "numpy",
            _OPENBLAS_REQUIREMENT,
        ]

    msg = f"Unsupported field: {field}"
    raise ValueError(msg)


def get_requires_for_dynamic_metadata(_settings: dict | None = None) -> list[str]:
    # get_version() only needs git and the standard library; netgen-mesher is
    # provided by the (non-isolated) build environment, not declared here.
    return []
