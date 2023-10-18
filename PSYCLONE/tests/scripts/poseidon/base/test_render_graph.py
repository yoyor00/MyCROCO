##########################################################
#  CROCO PSYCLONE scripts, under CeCILL-C
#  From Sébastien Valat (INRIA) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
# python
import os
# internal
from scripts.poseidon.base.render_graph import RenderGraph
from scripts.poseidon.base.types import AccessMode

##########################################################
GRAPH_1='''graph CROCO
{
node_0 [label="A"]
node_1 [label="B"]
node_0 -- node_1 
}
'''

##########################################################
GRAPH_2='''graph CROCO
{
node_0 [label="ROOT"]
node_1 [label="CALL func_name", shape=diamond]
node_0 -- node_1 
node_2 [label="IF", shape=triangle]
node_1 -- node_2 
node_3 [label="LOOP i", shape=doublecircle]
node_2 -- node_3 
node_4 [label="SCHEDULE", shape=box]
node_3 -- node_4 
}
'''

##########################################################
GRAPH_3='''graph CROCO
{
node_0 [label="A"]
node_1 [label="A", color="red"]
node_2 [label="A", color="green"]
}
'''

##########################################################
def test_basic_node_link():
    graph = RenderGraph()
    node_a = graph.add_node("A")
    node_b = graph.add_node("B")
    graph.add_link(node_a, node_b)
    assert graph.render() == GRAPH_1

##########################################################
def test_various_nodes():
    graph = RenderGraph()
    node_root = graph.add_node("ROOT")
    node_call = graph.add_call("func_name", parent_name=node_root)
    node_if = graph.add_if(parent_name=node_call)
    node_loop = graph.add_loop("i", parent_name=node_if)
    graph.add_schedule(parent_name=node_loop)
    assert graph.render() == GRAPH_2

##########################################################
def test_gen_svg():
    graph = RenderGraph()
    node_a = graph.add_node("A")
    graph.add_node("B", parent_name=node_a)
    filename = "/tmp/test_poseidon_render_graph_gen_svg.svg"
    graph.render_as_image(filename, "svg")
    os.unlink(filename)

##########################################################
def test_mode():
    graph = RenderGraph()
    graph.add_node("A", mode=AccessMode.UNDEFINED)
    graph.add_node("A", mode=AccessMode.WRITE)
    graph.add_node("A", mode=AccessMode.READ)
    assert graph.render() == GRAPH_3
