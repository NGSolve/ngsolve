import sys
import os
from subprocess import check_output

cwd = '.'
if len(sys.argv)>1:
    cwd = sys.argv[1]

git_version = check_output(['git', 'describe', '--tags'], cwd=cwd).decode('utf-8').strip()

def is_dev_build():
    if 'NG_NO_DEV_PIP_VERSION' in os.environ:
        return False
    if 'CI_COMMIT_REF_NAME' in os.environ and os.environ['CI_COMMIT_REF_NAME'] == 'release':
        return False
    return True

version = git_version[1:].split('-')
if len(version)>2:
    version = version[:2]
if len(version)>1:
    version = '.post'.join(version)
    if is_dev_build():
        version += '.dev0'
else:
    version = version[0]

print(version)

