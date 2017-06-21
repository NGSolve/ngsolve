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

# create dir to copy html files to
subprocess.run(['mkdir','../html_i-tutorials'])

# copy all html files to dir for upload
htmls = [f for f in glob.iglob('**/*.html', recursive=True)]
print(htmls)
for f in htmls:
    subprocess.run(['cp',f,'--parents','../html_i-tutorials'])

# zip all ipynb files
nbs = [f for f in glob.iglob('**/*.ipynb', recursive=True)]
subprocess.run(['zip','../html_i-tutorials/i-tutorials.zip',*nbs])
