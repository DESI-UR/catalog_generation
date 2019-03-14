
import os
import re

file_path = os.path.abspath(os.path.expandvars(os.path.expanduser("./py/paramock/_version.py")))
print(file_path)

def get_version(out_type=None):
    found_version = False
    with open(file_path, "r") as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        version_info = re.match("__version__.*=.*([0-9]+)\.([0-9]+)\.([0-9]+).*", line)
        major, minor, patch = version_info.groups()
        found_version = True
    if found_version is False:
        raise RuntimeError("Version information could not be accessed.")
    if out_type == 'string':
        return "{}.{}.{}".format(major, minor, patch)
    return int(major), int(minor), int(patch)

def update_version():
    major, minor, patch = get_version()
    print(type(major))
    with open(file_path, "w+") as fid:
        fid.write("__version__ = '{}.{}.{}'".format(major, minor, patch+1))

if __name__=="__main__":
    update_version()
