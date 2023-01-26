import json

policy_file = "/opt/_internal/pipx/venvs/auditwheel/lib/python3.10/site-packages/auditwheel/policy/manylinux-policy.json"
data = json.load(open(policy_file))
# todo: read libs from netgen.libs directory and remove hashes from names, then add just libmkl_rt.so.2
additional_libs = [
        "libGLU.so.1",
        "libGLX.so.0",
        "libGLdispatch.so.0",
        "libOpenGL.so.0",
        "libXmu.so.6",
        "libbz2.so.1.0.6",
        "libcsg.so",
        "libcsgvis.so",
        "libfontconfig.so.1.11.1",
        "libfreetype.so.6.14.0",
        "libgeom2d.so",
        "libgeom2dvis.so",
        "libgui.so",
        "libinterface.so",
        "libmesh.so",
        "libmkl_rt.so.2",
        "libngcore.so",
        "libnggui.so",
        "libnglib.so",
        "libocc.so",
        "libpng15.so.15.13.0",
        "libstl.so",
        "libstlvis.so",
        "libtcl8.so",
        "libtk8.so",
        "libuuid.so.1.3.0",
        "libvisual.so",
        "libz.so.1.2.7",
        ]

for entry in data:
    if 'manylinux' in entry['name']:
        entry['lib_whitelist'] += additional_libs
        
json.dump(data, open(policy_file, 'w'))
