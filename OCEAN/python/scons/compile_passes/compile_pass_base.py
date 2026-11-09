
from __future__ import annotations
from typing import List
from SCons.Node.FS import File
from abc import ABC, abstractmethod


class CompilePassBase(ABC):
    @abstractmethod
    def apply(self, pipe: List[File]) -> List[File]:
        """Apply the compile pass to a pipeline."""
        raise NotImplementedError
