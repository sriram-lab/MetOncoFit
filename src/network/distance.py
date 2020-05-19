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
for rxn in model.reactions:
    reactionNames.append(rxn.id)

metaboliteNames = []
for met in model.metabolites:
    metaboliteNames.append(met.id)

geneNames = []
for gene in model.genes:
    geneNames.append(gene.id)
print(geneNames)

biomassComponents = [] # <-- QUERY
biomassReactants = model.reactions.get_by_id('biomass_objective').reactants
for biocmp in biomassReactants:
    biomassComponents.append(biocmp.id)

# Get model stoichiometric matrix as an array
#stoich_matrix = cobra.util.create_stoichiometric_matrix(model)

#import networkx as nx
#import matplotlib.pyplot as plt

# Construct vertices for reactants and products
#import numpy as np
#product_rows, product_cols = np.where(stoich_matrix == 1)
#product_edges = zip(product_rows.tolist(), product_cols.tolist())

#reactant_rows, reactant_cols = np.where(stoich_matrix == -1)
#reactant_edges = zip(reactant_rows.tolist(), reactant_cols.tolist())

# Construct a directed graph for reactants and products
#prodGraph = nx.DiGraph(products=metaboliteNames, reactions=reactionNames)
#prodGraph.add_edges_from(product_edges)

#reactGraph = nx.DiGraph()
#reactGraph.add_edges_from(reactant_edges)

# Visualize the directed graphs
#nx.draw(prodGraph, node_size=5)
#plt.show()
#nx.draw(reactGraph, node_size=5)
#plt.show()

# Compute shortest distances for biomass components
#shortestBiomassPath = []
#for i in biomassComponents:
#    shortestBiomassPath[i] = nx.shortest_path_length(prodGraph, source=i, target=j)


#dijkstra_length = nx.shortest_path(gr)

