import sys
from subprocess import check_output

cwd = '.'
if len(sys.argv)>1:
    cwd = sys.argv[1]

git_version = check_output(['git', 'describe', '--tags'], cwd=cwd).decode('utf-8').strip()

version = git_version[1:].split('-')
if len(version)>2:
    version = version[:2]
if len(version)>1:
    version = '.post'.join(version) + '.dev'
else:
    version = version[0]

print(version)

