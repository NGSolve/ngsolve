from glob import glob
import sys, os
from base64 import b64encode
import json
sdir = sys.argv[1]+'/'
bdir = sys.argv[2]+'/'
code = open(sdir+"webgui_template.py","r").read()

jscode = open(bdir+"standalone.js","r").read()
code = code.replace('render_js_code = ""', 'render_js_code = r"""'+ jscode + '"""')

widgets_version = '^'+json.load(open(sdir+'/js/package.json'))['version']
code = code.replace('widgets_version = ""', 'widgets_version = "'+ widgets_version + '"')

open("webgui.py","w").write(code)
