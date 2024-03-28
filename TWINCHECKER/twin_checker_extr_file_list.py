#!/usr/bin/env python3

# python
import sys
import json5

if __name__ == "__main__":
    with open(sys.argv[1], "r") as fp:
        json = json5.load(fp)
        files = json['files']
        for group_name in json['specific']:
            group = json['specific'][group_name]
            files += group['files']

        print(' '.join(set(files)))
