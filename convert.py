import glob
import subprocess

procs = []

# convert all ipynb files to html
for f in glob.iglob('**/*.ipynb', recursive=True):
    procs.append(subprocess.Popen(['jupyter','nbconvert','--to','html',f]))

# wait for processes to finish
for p in procs:
    p.wait()

# replace links to .ipynb files by .html
for f in glob.iglob('**/*.html', recursive=True):
    subprocess.run(['sed','-i','s/\.ipynb/\.html/g',f])
