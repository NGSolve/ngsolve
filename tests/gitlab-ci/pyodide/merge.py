import hashlib
import json
import os

pyodide_version = os.environ['PYODIDE_VERSION']

repo_file = 'pyodide-lock.json'
ori = json.load(open(repo_file))
pkg = ori['packages']

for name in pkg:
    p = pkg[name]
    if not p['file_name'].startswith('https://') and not os.path.exists(os.path.join('pyodide', p['file_name'])):
        p['file_name'] = f"https://cdn.jsdelivr.net/pyodide/v{pyodide_version}/full/"+p['file_name']
    pkg[name] = p

pkg['webgui-jupyter-widgets']['depends'].append('numpy')

def getHash(filename):
    with open(filename,"rb") as f:
        return hashlib.sha256(f.read()).hexdigest()

pkg.update({
    "ipykernel": {
      "name": "ipykernel",
      "version": "6.9.2",
      "file_name": "../../extensions/@jupyterlite/pyodide-kernel-extension/static/pypi/./ipykernel-6.9.2-py3-none-any.whl",
      "install_dir": "site",
      "sha256": "9525bc401a9257575e3564669404a0680276bb961538a3242bc191e35d7e2c75",
      "imports": [],
      "depends": []
    },
    "piplite": {
      "name": "piplite",
      "version": "0.4.2",
      "file_name": "../../extensions/@jupyterlite/pyodide-kernel-extension/static/pypi/piplite-0.4.2-py3-none-any.whl",
      "install_dir": "site",
      "sha256": "e508c140070243c175173ae64f8dda2d150276ea8fbdfd4f0af58c884fa0b30c",
      "imports": [],
      "depends": []
    },
    "pyodide-kernel": {
      "name": "pyodide-kernel",
      "version": "0.4.2",
      "file_name": "../../extensions/@jupyterlite/pyodide-kernel-extension/static/pypi/./pyodide_kernel-0.4.2-py3-none-any.whl",
      "install_dir": "site",
      "sha256": "6dd6f1d5124a77fa973b503cff9d8b0154c193022d2c5bd1d294fe2c1b39515c",
      "imports": [],
      "depends": []
    },
    "widgetsnbextension": {
      "name": "widgetsnbextension",
      "version": "4.0.999",
      "file_name": "../../extensions/@jupyterlite/pyodide-kernel-extension/static/pypi/./widgetsnbextension-4.0.999-py3-none-any.whl",
      "install_dir": "site",
      "sha256": "9d88fa0e2e496413454a4bc215d2af9f0c33612ce6bb459bd1a61f4eec59358c",
      "imports": [],
      "depends": []
    },
    "netgen": {
      "name": "netgen",
      "version": "6.2.2406",
      "file_name": "netgen.zip",
      "install_dir": "stdlib",
      "sha256": getHash("netgen.zip"),
      "package_type": "cpython_module",
      "imports": [
        "netgen"
      ],
      "depends": [
        "pyngcore",
        "webgui-jupyter-widgets"
      ],
      "unvendored_tests": False,
      "shared_library": True
    },
    "ngsolve": {
      "name": "ngsolve",
      "version": "6.2.2406",
      "file_name": "ngsolve.zip",
      "install_dir": "stdlib",
      "sha256": getHash("ngsolve.zip"),
      "package_type": "cpython_module",
      "imports": [
        "ngsolve"
      ],
      "depends": [
        "openblas",
        "numpy",
        "netgen"
      ],
      "unvendored_tests": False,
      "shared_library": True
    },
    "pyngcore": {
      "name": "pyngcore",
      "version": "6.2.2406",
      "file_name": "pyngcore.zip",
      "install_dir": "stdlib",
      "sha256": getHash("pyngcore.zip"),
      "package_type": "cpython_module",
      "imports": [
        "pyngcore"
      ],
      "depends": [],
      "unvendored_tests": False,
      "shared_library": True
    }
})

json.dump(ori, open(repo_file, "w"), indent=2)
