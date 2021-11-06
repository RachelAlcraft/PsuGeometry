'''
from pyvis.network import Network
net = Network(notebook=True)
net.add_node(1, label='Alex')
net.add_node(2, label='Cathy')
#net.from_nx(G)
net.show('nodes.html')
'''
from pyvis import network as net
import networkx as nx
g = net.Network('800px', '1000px',notebook=True)
g.add_node(0, label='Birkbeck',color='Crimson',value=200, shape='triangle')
g.add_node(1, label='UCL',color='Crimson',value=200, shape='triangle')
g.show('nodes.html')