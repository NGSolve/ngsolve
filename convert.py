import glob
import subprocess
import os.path

procs = []

dest = '../docs/i-tutorials'

def ProcessNotebooks(args):
    for f in glob.iglob('**/*.ipynb', recursive=True):
        procs.append(subprocess.Popen([*args,f]))

    # wait for processes to finish
    for p in procs:
        p.wait()

# clear all output in ipynb files
ProcessNotebooks(['jupyter','nbconvert','--to','notebook', '--inplace', '--ClearOutputPreprocessor.enabled=True'])

# convert all ipynb files to html
ProcessNotebooks(['jupyter','nbconvert','--to','html'])

# replace links to .ipynb files by .html
for f in glob.iglob('**/*.html', recursive=True):
    if not f.endswith('before_we_start.html'):
        subprocess.run(['sed','-i','s/\.ipynb/\.html/g',f])

# create dir to copy html files to
subprocess.run(['mkdir',dest])

# zip all ipynb files
nbs = [f for f in glob.iglob('**/*.ipynb', recursive=True)]
subprocess.run(['zip', dest+'/i-tutorials.zip',*nbs])

# copy all html files and resources (e.g. png) to destination
all_files = [f for f in glob.iglob('**', recursive=True)]
all_files = filter(lambda f: not os.path.isdir(f), all_files)
all_files = filter(lambda f: not f.endswith('.ipynb'), all_files)
all_files = filter(lambda f: not f.endswith('.py'), all_files)
for f in all_files:
    subprocess.run(['cp',f,'--parents',dest])
