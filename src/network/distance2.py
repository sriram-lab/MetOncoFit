"""
Compute topological and medium epicenter distances by transforming a COBRA model into a network graph.

@author: Scott Campit
"""
import re

import networkx as nx
import cobra
import matplotlib.pyplot as plt

recon1Graph = nx.DiGraph()

modelFile = '/home/scampit/Data/CBM/MetabolicModels/RECON1/model_human_duarte.mat'
model = cobra.io.load_matlab_model(modelFile)

with open('/home/scampit/software/MetOncoFit/srv/cofactors.txt') as f:
    cofactors = f.readlines()
cofactors = [x.strip() for x in cofactors]

# Decompose network by iterating through reactions
for rxn in model.reactions:
    rxnFormula = re.sub(r'\s*\d+\.\d+\s*', '', rxn.reaction)
    # Get reactions, productions, and network directionality
    if "-->" in rxnFormula:
        reactants, products = rxnFormula.split('-->')
        reactionType = 'IRREVERSIBLE'
    elif "<=>" in rxnFormula:
        reactants, products = rxnFormula.split('<=>')
        reactionType = 'REVERSIBLE'
    # Construct directed network by products
    for prod in products.split('+'):
        for react in reactants.split('+'):
            if (prod.strip() or react.strip()) not in cofactors:
                if reactionType is 'IRREVERSIBLE':
                    recon1Graph.add_edge(u_of_edge=react.strip(), v_of_edge=prod.strip(), label=rxn.id)
                elif reactionType is 'REVERSIBLE':
                    recon1Graph.add_edge(u_of_edge=react.strip(), v_of_edge=prod.strip(), label=rxn.id)
                    recon1Graph.add_edge(u_of_edge=prod.strip(), v_of_edge=react.strip(), label=rxn.id)

# Sanity check 1: Draw out RECON1 as a graph --> hairball, as you would expect
#nx.draw(recon1, node_size=5)
#plt.show()

# Get list of metabolites to query -> biomass and medium components
with open('/home/scampit/software/MetOncoFit/srv/biomass.txt') as f:
    biomass = f.readlines()
biomass = [x.strip() for x in biomass]
with open('/home/scampit/software/MetOncoFit/srv/media.txt') as f:
    media = f.readlines()
media = [x.strip() for x in media]

biomass_rxn = list()
medium_rxn = list()
biomass_dist = list()
medium_dist = list()
# Get the shortest distance from the medium or biomass component to any given reaction.
for node in recon1Graph:
    for b in biomass:
        if (node.strip() or b.strip()) not in cofactors:
            try:
                # Distance from any given reaction to a biomass component
                biomass_dist.append(len(nx.shortest_path(G=recon1Graph, source=node.strip(), target=b.strip()))  )
                biomass_rxn.append(recon1Graph.out_edges(node.strip()))

            except: # If doesn't exist, set to 0
                biomass_rxn.append(recon1Graph.out_edges(node.strip()))
                biomass_dist.append(0)

    for m in media:
        if (node.strip() or m.strip()) not in cofactors:
            try:
                # Distance from any given medium component to a reaction
                medium_dist.append(len(nx.shortest_path(G=recon1Graph, source=m.strip(), target=node.strip())))
                medium_rxn.append(recon1Graph.in_edges(node.strip()))
            except: # If doesn't exist, set to 0
                medium_rxn.append(recon1Graph.in_edges(node.strip()))
                medium_dist.append(0)
print(medium_rxn)


