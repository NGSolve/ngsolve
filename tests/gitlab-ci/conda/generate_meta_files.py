import subprocess
import os

# CI_COMMIT_TAG is set when a git tag is built -> a new release
# otherwise build a nightly package
if "CI_COMMIT_TAG" in os.environ:
    package_name_suffix = ''
else:
    package_name_suffix = '-nightly'

def parseVersion(s):
    s = s[1:]
    tok = s.split('-')
    v = tok[0].strip()
    if len(tok)==1:
        return v,v+".0"
    else:
        v += "."+tok[1].strip()
        return v,v


# get current version descriptions and generate meta.yaml with these versions
subprocess.run(["git", "remote", "update"], cwd="../../../external_dependencies/netgen")
v_netgen,v_netgen_dep = parseVersion(subprocess.run(["git", "describe", "--tags"], capture_output=True, cwd="../../../external_dependencies/netgen").stdout.decode())
print("Netgen version", v_netgen)

meta = open("netgen/meta_template.yaml","r").read() \
        .replace("{{version}}", v_netgen)     \
        .replace("{{suffix}}", package_name_suffix)
open("netgen/meta.yaml","w").write(meta)

v_ngs, dummy = parseVersion(subprocess.run(["git", "describe", "--tags"], capture_output=True, cwd="../../..").stdout.decode())
print("NGSolve version", v_ngs)

meta = open("ngsolve/meta_template.yaml","r").read()   \
        .replace("{{version}}", v_ngs)       \
        .replace("{{netgen_version}}", v_netgen_dep) \
        .replace("{{suffix}}", package_name_suffix)
open("ngsolve/meta.yaml","w").write(meta)
