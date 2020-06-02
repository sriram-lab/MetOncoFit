"""
Distance.py
Calculate topological and medium epicenter distances

TODO:
  * Create a script that re-computes the topological and medium epicenter scores
  * Compute other distance metrics other than shortest distance

@author: Scott Campit
"""

# Load RECON1
import cobra
modelFile = '/home/scampit/Data/CBM/MetabolicModels/RECON1/model_human_duarte.mat'
model = cobra.io.load_matlab_model(modelFile)

# Get genes, reactions, metabolite names, and biomass components from the metabolic model
reactionNames = []
tmp = []
rxnGeneObj = []
for rxn in model.reactions:
    reactionNames.append(rxn.id)
#    tmp.append(rxn.gene_reaction_rule)
#print(tmp)

metaboliteNames = []
for met in model.metabolites:
    metaboliteNames.append(met.id)

geneNames = []
for gene in model.genes:
    geneNames.append(gene.id)
#print(geneNames)


# Get list to query
#with open('/home/scampit/software/MetOncoFit/srv/biomass.txt') as f:
#    biomass = f.readlines()
#biomass = [x.strip() for x in biomass]

#with open('/home/scampit/software/MetOncoFit/srv/media.txt') as f:
#    media = f.readlines()
#media = [x.strip() for x in media]
#to_query = biomass + media


# Get model stoichiometric matrix as an array
stoich_matrix = cobra.util.create_stoichiometric_matrix(model)

import networkx as nx
import matplotlib.pyplot as plt

# Construct vertices for reactants and products
import numpy as np
product_rows, product_cols = np.where(stoich_matrix == 1)
product_edges = zip(product_rows.tolist(), product_cols.tolist())

reactant_rows, reactant_cols = np.where(stoich_matrix == -1)
reactant_edges = zip(reactant_rows.tolist(), reactant_cols.tolist())

# Construct a directed graph for reactants and products
prodGraph = nx.DiGraph(products=metaboliteNames, reactions=reactionNames)
prodGraph.add_edges_from(product_edges)

reactGraph = nx.DiGraph()
reactGraph.add_edges_from(reactant_edges)

# Visualize the directed graphs
#nx.draw(prodGraph, node_size=5)
#plt.show()
#nx.draw(reactGraph, node_size=5)
#plt.show()

# Compute shortest distances for biomass components
#shortestBiomassPath = []
#for i in biomassComponents:
#    shortestBiomassPath[i] = nx.shortest_path_length(prodGraph, source=i, target=j)
print([p for p in nx.all_shortest_paths(prodGraph, source=1, target=2)])

#dijkstra_length = nx.shortest_path(gr)

