import re
from subprocess import run, PIPE, STDOUT

deps = snakemake.params[0]

OKBLUE = '\033[94m'
ENDC = '\033[0m'
FAIL = '\033[91m'

def dep_check(dep):
        try:
                output = run([ dep["command"], dep["version_arg"] ], stdout = PIPE, stderr = STDOUT)
                lines = output.stdout.decode('utf-8').splitlines()
                version = lines.pop(0)
                if lines and re.search('version', lines[0], re.IGNORECASE): version = lines.pop(0)
                print(OKBLUE + "Using %s: %s" % (dep["command"], version) + ENDC, file = sys.stderr)
                return True
        except FileNotFoundError as e:
                print(FAIL + "%s not found on path, download from %s" % (dep["name"], dep["url"]) + ENDC, file = sys.stderr)
                return False

ok = True
for dep in deps:
    ok *= dep_check(dep)

if not ok:
    exit("Some of the dependencies missing")
else:
    with open('data/dependencies-cli.txt', 'w') as fp:
        fp.write('OK')
