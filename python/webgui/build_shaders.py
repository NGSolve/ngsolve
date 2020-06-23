from glob import glob
import sys, os
from base64 import b64encode
import json
p = sys.argv[1]+'/'

shader_codes = {}
for f in glob(p+"shader/*"):
    name = os.path.basename(f)
    source = b64encode(open(f,"r").read().encode("ascii")).decode("ascii")
    shader_codes[name] = source

source = "export const shaders = " + json.dumps(shader_codes)
open(p+"/js/src/shaders.ts","w").write(source)
