
from __future__ import annotations
import os
import re
from typing import List
import sys

from SCons.Node.FS import FS


class Preprocessor:
    def __init__(self, node: FS.File):
        self.node = node

    def gen_deps_exists(self, paths: List[FS.Dir], inc_files=None):
        """
        Generate found include dependencies which are also in `paths`

        :param self: Description
        :param path: Description
        """

        if hasattr(self.node.attributes, "cache_preprocess_inc_files"):
            return self.node.attributes.cache_preprocess_inc_files

        if not inc_files:
            # Note that this is internally cached in `re` module
            inc_regex = r'^\s*#\s*include\s*["\']([^"\']+)["\']'
            regexp = re.compile(inc_regex, re.M)

            if not hasattr(self.node.attributes, "cache_content"):
                self.node.attributes.contents = self.node.get_text_contents()
            contents = self.node.attributes.contents

            inc_files = regexp.findall(
                contents,
                re.MULTILINE,
            )

        # Search for matching file in paths
        existing_inc_deps = []

        for inc in inc_files:
            found = False
            for dirname in paths:
                file = str(dirname) + os.sep + inc
                if os.path.exists(file):
                    if file.startswith("./"):
                        file = file[2:]
                    existing_inc_deps.append(file)
                    found = True
                    break

            if not found:
                # Ignore dependencies which we did not find, but make the
                # omission visible to the user.
                print(
                    "WARNING: Preprocessor scanner - File "
                    f"'{inc}' not found in paths "
                    f"'{[str(_) for _ in paths]}'",
                    file=sys.stderr,
                )

        self.node.attributes.cache_preprocess_inc_files = existing_inc_deps
        return existing_inc_deps
