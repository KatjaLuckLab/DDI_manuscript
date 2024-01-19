
import os

def find_file(filename):
    for dirpath, dirs, files in os.walk('/home/milenadjokic/NDD_project/'):
        for fname in dirs + files:
            if fname == filename:
                filepath = os.path.join(dirpath, fname)
    return filepath
