"""
Compute topological and medium epicenter distances by transforming a COBRA model into a network graph.

@author: Scott Campit
"""

# Built-in libraries
import re

# Downloaded libraries
import networkx as nx
import pandas as pd
import cobra

# Visualization libraries
import matplotlib.pyplot as plt

recon1Graph = nx.DiGraph()

modelFile = '/home/scampit/Data/CBM/MetabolicModels/RECON1/model_human_duarte.mat'
model = cobra.io.load_matlab_model(modelFile)

with open('/home/scampit/software/MetOncoFit/srv/cofactors.txt') as f:
    cofactors = f.readlines()
cofactors = [x.strip() for x in cofactors]

# Construct network by iterating through reactions
print("Constructing graph network from COBRA model")
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
biomass_component = list()
medium_component = list()
# Get the shortest distance from the medium or biomass component to any given reaction.
print("Computing medium and biomass epicenter scores")
for node, _, metadata in recon1Graph.edges(data=True):
    for b in biomass:
        if (node.strip() or b.strip()) not in cofactors:
            try:
                # Distance from any given reaction to a biomass component
                biomass_dist.append(len(nx.shortest_path(G=recon1Graph, source=node.strip(), target=b.strip()))  )
                biomass_rxn.append(metadata['label'])
                biomass_component.append(b.strip())

            except: # If doesn't exist, set to 0
                biomass_rxn.append(metadata['label'])
                biomass_dist.append(0)
                biomass_component.append(b.strip())

    for m in media:
        if (node.strip() or m.strip()) not in cofactors:
            try:
                # Distance from any given medium component to a reaction
                medium_dist.append(len(nx.shortest_path(G=recon1Graph, source=m.strip(), target=node.strip())))
                medium_rxn.append(metadata['label'])
                medium_component.append(m.strip())

            except: # If doesn't exist, set to 0
                medium_rxn.append(metadata['label'])
                medium_dist.append(0)
                medium_component.append(m.strip())

# Construct dataframes for the distance metrics
print("Mapping scores from reactions to genes")
biomass_df = pd.DataFrame({"Distance": biomass_dist,
                           "Reaction": biomass_rxn,
                           'Biomass': biomass_component})
medium_df = pd.DataFrame({"Distance": medium_dist,
                          "Reaction": medium_rxn,
                          'Medium': medium_component})

# Get the minimum shortest distance that is not 0 based on reaction-medium pairs --> bottleneck
medium_df = medium_df.groupby(['Reaction', 'Medium']).apply(lambda x: x[x > 0].min(axis=1))
medium_df = medium_df.reset_index().fillna(0)

biomass_df = biomass_df.groupby(['Reaction', 'Biomass']).apply(lambda x: x[x > 0].min(axis=1))
biomass_df = biomass_df.reset_index().fillna(0)

# Map back to genes
geneRxnMap = pd.read_csv('/home/scampit/software/MetOncoFit/srv/RECON1ReactionGeneMap.txt')
medium_df = pd.merge(medium_df, geneRxnMap,
                     how='inner',
                     left_on='Reaction', right_on='Reaction')
biomass_df = pd.merge(biomass_df, geneRxnMap,
                     how='inner',
                     left_on='Reaction', right_on='Reaction')
print(biomass_df)

# Pivot so that genes are rows and medium components are columns
print("Construct final array and save to file")
#biomass_df = biomass_df.pivot(index='Gene', columns='Biomass', values=0)
#medium_df = medium_df.pivot(index='Gene', columns='Medium', values=0)
biomass_df.to_csv('/home/scampit/software/MetOncoFit/srv/biomass_distances.csv')
medium_df.to_csv('/home/scampit/software/MetOncoFit/srv/medium_distances.csv')



