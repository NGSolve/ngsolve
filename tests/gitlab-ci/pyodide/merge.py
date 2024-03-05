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
      "sha256": "83d4bd34f31029710294f53e6f11ac06f48f6777098a5b8cb0ae9aa46fb38764",
      "imports": [],
      "depends": []
    },
    "piplite": {
      "name": "piplite",
      "version": "0.3.0a0",
      "file_name": "../../extensions/@jupyterlite/pyodide-kernel-extension/static/pypi/piplite-0.3.0a0-py3-none-any.whl",
      "install_dir": "site",
      "sha256": "20734a5170829f6b8b56179d601f10e36a6289f00a85eb0ded9c7a37d085f982",
      "imports": [],
      "depends": []
    },
    "pyodide-kernel": {
      "name": "pyodide-kernel",
      "version": "0.3.0a0",
      "file_name": "../../extensions/@jupyterlite/pyodide-kernel-extension/static/pypi/./pyodide_kernel-0.3.0a0-py3-none-any.whl",
      "install_dir": "site",
      "sha256": "9d9a78d0305b359cc494b7343257e12049c5753819f3fe7c82bc9f4e2424de69",
      "imports": [],
      "depends": []
    },
    "widgetsnbextension": {
      "name": "widgetsnbextension",
      "version": "4.0.10",
      "file_name": "../../extensions/@jupyterlite/pyodide-kernel-extension/static/pypi/./widgetsnbextension-4.0.10-py3-none-any.whl",
      "install_dir": "site",
      "sha256": "28b5e8b2a9fa39d187b88e4f5dea10441b0cd3a5afc4a3b8d452e9c256f22c46",
      "imports": [],
      "depends": []
    },
    "netgen": {
      "name": "netgen",
      "version": "6.2.2401",
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
      "version": "6.2.2401",
      "file_name": "ngsolve.zip",
      "install_dir": "stdlib",
      "sha256": getHash("ngsolve.zip"),
      "package_type": "cpython_module",
      "imports": [
        "ngsolve"
      ],
      "depends": [
        "numpy",
        "netgen"
      ],
      "unvendored_tests": False,
      "shared_library": True
    },
    "pyngcore": {
      "name": "pyngcore",
      "version": "0.0.1",
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
