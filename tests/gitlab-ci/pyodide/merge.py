import hashlib
import json
import os

pyodide_version = os.environ["PYODIDE_VERSION"]

repo_file = "pyodide-lock.json"
ori = json.load(open(repo_file))
pkg = ori["packages"]

for name in pkg:
    p = pkg[name]
    if not p["file_name"].startswith("https://") and not os.path.exists(
        os.path.join("pyodide", p["file_name"])
    ):
        p["file_name"] = (
            f"https://cdn.jsdelivr.net/pyodide/v{pyodide_version}/full/"
            + p["file_name"]
        )
    pkg[name] = p

pkg["webgui-jupyter-widgets"]["depends"].append("numpy")


def getHash(filename):
    with open(filename, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


pkg.update(
    {
        "widgetsnbextension": {
            "version": "4.0.999",
            "file_name": "../../extensions/@jupyterlite/pyodide-kernel-extension/static/pypi/./widgetsnbextension-4.0.999-py3-none-any.whl",
            "install_dir": "site",
            "sha256": "07fb49113ea4f6c943684b67b73a4f4209f4194ef97658d5285adffe51acfe8a",
            "imports": [],
            "depends": [],
        },
        "netgen": {
            "name": "netgen",
            "version": "6.2.2406",
            "file_name": "netgen.zip",
            "install_dir": "stdlib",
            "sha256": getHash("netgen.zip"),
            "package_type": "cpython_module",
            "imports": ["netgen"],
            "depends": ["pyngcore", "webgui-jupyter-widgets"],
            "unvendored_tests": False,
            "shared_library": True,
        },
        "ngsolve": {
            "name": "ngsolve",
            "version": "6.2.2406",
            "file_name": "ngsolve.zip",
            "install_dir": "stdlib",
            "sha256": getHash("ngsolve.zip"),
            "package_type": "cpython_module",
            "imports": ["ngsolve"],
            "depends": ["openblas", "numpy", "netgen"],
            "unvendored_tests": False,
            "shared_library": True,
        },
        "pyngcore": {
            "name": "pyngcore",
            "version": "6.2.2406",
            "file_name": "pyngcore.zip",
            "install_dir": "stdlib",
            "sha256": getHash("pyngcore.zip"),
            "package_type": "cpython_module",
            "imports": ["pyngcore"],
            "depends": [],
            "unvendored_tests": False,
            "shared_library": True,
        },
    }
)

json.dump(ori, open(repo_file, "w"), indent=2)
