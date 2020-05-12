from glob import glob
import sys, os
from base64 import b64encode
p = sys.argv[1]+'/'
code = open(p+"webgui.py","r").read()

shader_codes = {}
for f in glob(p+"shader/*"):
    name = os.path.basename(f)
    source = b64encode(open(f,"r").read().encode("ascii")).decode("ascii")
    shader_codes[name] = source

jscode = open(p+"ngsolve_webgui.js","r").read()
jscode = jscode.replace("const shaders = {}", "const shaders = "+str(shader_codes))
code = code.replace('render_js_code = ""', 'render_js_code = """'+ jscode + '"""')

open("webgui.py","w").write(code)
