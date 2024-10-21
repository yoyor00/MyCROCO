##########################################################
#  CROCO PSYCLONE scripts, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
'''
Implement a simple graphviz helper to generate graphs.
'''

##########################################################
# python
import tempfile
import subprocess
# internal
from .types import AccessMode

##########################################################
class RenderGraphNode:
    def __init__(self, name: str):
        self._name = name

    @property
    def name(self):
        return self._name

##########################################################
class RenderGraph:
        def __init__(self, digraph: bool = False):
            self.next_id = 0
            self.digraph = digraph
            if digraph:
                self.lines = [
                    'digraph CROCO',
                    '{',
                ]
            else:
                self.lines = [
                    'graph CROCO',
                    '{',
                ]

        def add_schedule(self, parent_name: RenderGraphNode=None) -> RenderGraphNode:
            id = f"node_{self.next_id}"
            self.next_id += 1
            self.lines.append(f"{id} [label=\"SCHEDULE\", shape=box]")
            self.add_link(parent_name, RenderGraphNode(id))
            return RenderGraphNode(id)

        def add_call(self, name: str, parent_name: RenderGraphNode=None) -> RenderGraphNode:
            id = f"node_{self.next_id}"
            self.next_id += 1
            self.lines.append(f"{id} [label=\"CALL {name}\", shape=diamond]")
            self.add_link(parent_name, RenderGraphNode(id))
            return RenderGraphNode(id)

        def add_loop(self, variable, parent_name: RenderGraphNode=None) -> RenderGraphNode:
            id = f"node_{self.next_id}"
            self.next_id += 1
            self.lines.append(f"{id} [label=\"LOOP {variable}\", shape=doublecircle]")
            self.add_link(parent_name, RenderGraphNode(id))
            return RenderGraphNode(id)

        def add_kernel(self, name, parent_name: RenderGraphNode=None) -> RenderGraphNode:
            id = f"node_{self.next_id}"
            self.next_id += 1
            self.lines.append(f"{id} [label=\"{name}\", shape=box]")
            self.add_link(parent_name, RenderGraphNode(id))
            return RenderGraphNode(id)

        def add_if(self, parent_name: RenderGraphNode=None) -> RenderGraphNode:
            id = f"node_{self.next_id}"
            self.next_id += 1
            self.lines.append(f"{id} [label=\"IF\", shape=triangle]")
            self.add_link(parent_name, RenderGraphNode(id))
            return RenderGraphNode(id)

        def add_node(self, title: str, parent_name: RenderGraphNode=None, mode:AccessMode=AccessMode.UNDEFINED) -> RenderGraphNode:
            id = f"node_{self.next_id}"
            self.next_id += 1
            if mode == AccessMode.UNDEFINED:
                self.lines.append(f"{id} [label=\"{title}\"]")
            elif mode == AccessMode.READ or mode == 'R':
                self.lines.append(f"{id} [label=\"{title}\", color=\"green\"]")
            elif mode == AccessMode.WRITE or mode == 'W':
                self.lines.append(f"{id} [label=\"{title}\", color=\"red\"]")
            elif mode == AccessMode.READ_WRITE or mode == 'RW':
                self.lines.append(f"{id} [label=\"{title}\", color=\"blue\"]")
            else:
                raise Exception(f"Access mode not supported : {mode}")
            self.add_link(parent_name, RenderGraphNode(id))
            return RenderGraphNode(id)

        def add_link(self, source:RenderGraphNode, dest:RenderGraphNode, style=''):
            if style != '':
                style = f"[{style}]"
            if source != None:
                if self.digraph:
                    self.lines.append(f"{source.name} -> {dest.name} {style}")
                else:
                    self.lines.append(f"{source.name} -- {dest.name} {style}")

        def open_subgraph(self, name:str):
            id = f"clusterstep1{self.next_id}"
            self.next_id += 1
            self.lines.append(f"subgraph {id} {{")
            self.lines.append(f"label=\"{name}\"")
            self.lines.append(f"color=black")

        def close_subgraph(self):
            self.lines.append(" }")

        def render(self):
            return "\n".join(self.lines)+"\n}\n"

        def render_as_image(self, filename: str, format: str):
            with tempfile.NamedTemporaryFile() as fp:
                # gen .dot file
                fp.write(self.render().encode())
                fp.flush()

                # gen requested file with dot command
                subprocess.run(f"dot -T{format} -o{filename} {fp.name}", shell=True, check=True)
