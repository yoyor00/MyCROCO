#!/usr/bin/env python3

from enum import Enum

class AccessMode(Enum):
    UNDEFINED = 0
    READ = 1
    WRITE = 2
    READ_WRITE = 3
