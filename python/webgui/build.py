from glob import glob
import sys, os
p = sys.argv[1]+'/'
code = open(p+"webgui.py","r").read()

shader_codes = {}
for f in glob(p+"shader/*"):
    name = os.path.basename(f)
    shader_codes[name] = open(f,"r").read()

code = code.replace("shader_codes = {}", "shader_codes = " + str(shader_codes))
code = code.replace('render_js_code = ""', 'render_js_code = """'+ open(p+"ngsolve_webgui.js","r").read() + '"""')

open("webgui.py","w").write(code)
