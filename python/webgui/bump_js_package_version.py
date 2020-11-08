#! /usr/bin/python3

# 1) set new version in package.json
# 2) build ngsolve js package -> generates .tgz with new package
# 3) release widgets package with 'npm publish'
# 4) update SHA1 value in ngsolve_jupyter_widgets.sha1sum
# 5) commit changes in ngsolve_jupyter_widgets and package.json (done manually)

import sys, json, subprocess
package_json = json.load(open("js/package.json"))
old_version = package_json['version']

if len(sys.argv) < 2:
    print("Usage: {} new_version".format(sys.argv[0]))
    print("Current version: "+old_version)
    exit(1)

new_version = sys.argv[1]

def onError():
    package_json['version'] = old_version
    json.dump(package_json, open("js/package.json", "w"), indent=2)
    exit(1)

def runCommand(cmd):
    proc = subprocess.run(cmd, capture_output=True, cwd="js")
    if proc.returncode!=0:
        print('ERROR:')
        print(proc.stdout.decode())
        print(proc.stderr.decode())
        onError()

    return proc.stdout.decode()

answer = ''
while answer not in ['y','n']:
    answer = input('bumping version from {} to {}, are you sure? [y/n] '.format(old_version, new_version))

if new_version <= old_version:
    print("error: cannot decrease version number")
    exit(1)

if answer == 'y':
    package_json['version'] = new_version
    json.dump(package_json, open("js/package.json", "w"), indent=2)
    print("new version set to {}".format(new_version))
    print("building js package...", end=" ", flush=True)
    runCommand(['npm', 'run', 'build'])
    runCommand(['npm', 'pack'])
    print("done")
    print("publishing js package...", end=" ", flush=True)
    runCommand(['npm', 'publish'])
    print("done")
    print("update SHA1 sum in CMakeLists.txt...")
    sha1sum = runCommand(["sha1sum", "ngsolve_jupyter_widgets-{}.tgz".format(new_version)])
    open("js/ngsolve_jupyter_widgets.sha1sum","w").write(sha1sum.split(' ')[0])
    print("git add package.json ngsolve_jupyter_widgets.sha1sum")
    runCommand(["git", "add", "package.json", "ngsolve_jupyter_widgets.sha1sum"])
    print("Done! Run git commit to commit the version bump")

