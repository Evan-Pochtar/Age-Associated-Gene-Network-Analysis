import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
from bokeh.models.callbacks import CustomJS
from bokeh.layouts import row
from bokeh.io import output_file, show
from bokeh.plotting import figure, from_networkx
from bokeh.models import MultiLine, LabelSet, ColumnDataSource, Circle, Select
from bokeh import events
from bokeh.models import HoverTool
import os

# Global Vars, reading in data
# Initial Gene expression data, not used in final product
# gene_data.csv is too big for github, can be found at "https://doi.org/10.1111/acel.13280"
# Use txtToCSV.py to convert txt from article to csv
patient_data = pd.read_csv('patients.csv')
gene_data = pd.read_csv('gene_data.csv', index_col = 0)

# Geneset data from R program
geneset_data = pd.read_csv('geneset_data.csv', index_col = 0)
# Geneset PCC data from R program
geneset_pcc_data = pd.read_csv('healthy_tissue_independent_pcc_matrix.csv', index_col=0)

# Disease data, calcualted from R program
# diseases layout is: Disease : Geneset : regulation : geneset : regulation etc
significant_diseases = pd.read_csv('significant_matrix.csv', index_col = "GENESET")
diseases = significant_diseases.to_dict()

# JSD data used to create JSD network, given from R program
meanJSD_CM = pd.read_csv('mean_jsd_combined_matrix.csv', index_col = 0)
pvalJSD_CM = pd.read_csv('p_value_jsd_combined_matrix.csv', index_col = 0)

# Geneset data from R program
# Uses HOA labels of genes to convert into final HOA data
genePCCMatrix = pd.read_csv('gene_pcc_matrix.csv', index_col = 0)
hoa_component = pd.read_csv('hoa_component.csv', index_col = 0)
gene_normalized_hoa_component = pd.read_csv('gene_normalized_hoa_component.csv', index_col = 0)
adjusted_hoa_component = pd.read_csv('adjusted_hoa_component.csv', index_col = 0)

# Input: A geneset name and a threshold
# Output: Prints out correlations higher than the given threshold from Geneset PCC data
def describeGeneset(geneset, threshold=0.8):
    # Get the correlations for the given geneset
    correlations = geneset_pcc_data.loc[geneset]
    # Filter out the genesets that are not connected
    connected_genesets = correlations[correlations >= threshold]
    connected_genesets.sort_values(ascending=False, inplace=True)
    connected_genesets.drop(geneset, inplace=True)
    # Print Connected genes
    print(f"Genesets connected to {geneset} and their correlations:")
    for connected_geneset, correlation in connected_genesets.items():
        print(f"{connected_geneset}: {correlation}")

# Input: geneList is the list of gene names
# Output: List of all total pairs, as well as correlation between gene expression data
def getGeneCorrelation(geneList):
    def compute_correlation(pair):
        return gene_data.loc[pair[0]].corr(gene_data.loc[pair[1]], method='pearson')
    pairs = [(gene1, gene2) for i, gene1 in enumerate(geneList) for j, gene2 in enumerate(geneList) if i < j]
    correlations = Parallel(n_jobs=-1)(delayed(compute_correlation)(pair) for pair in pairs)
    return pairs, correlations

# Input: geneList is the list of geneset names
# Output: List of all total pairs, as well as correlation between coefficent estimate data
def getGeneCorrelation_JCC(geneList):
    def compute_correlation_JCC(pair):
        x = geneset_data.loc[pair[0]]['coefficient_estimate']
        y = geneset_data.loc[pair[1]]['coefficient_estimate']
        correlation = (x * y) / (((abs(x) + abs(y)) / 2) ** 2)
        return correlation
    pairs = [(gene1, gene2) for i, gene1 in enumerate(geneList) for j, gene2 in enumerate(geneList) if i < j]
    correlations = Parallel(n_jobs=-1)(delayed(compute_correlation_JCC)(pair) for pair in pairs)    
    return pairs, correlations

# Input: geneList is the list of geneset names
# Output: List of all total pairs, as well as correlations from geneset PCC data
def getGenesetPCC(geneList):
    def getPCCfromData(pair):
        return geneset_pcc_data.loc[pair[0]][pair[1]]
    pairs = [(gene1, gene2) for i, gene1 in enumerate(geneList) for j, gene2 in enumerate(geneList) if i < j]
    correlations = Parallel(n_jobs=-1)(delayed(getPCCfromData)(pair) for pair in pairs)
    return pairs, correlations

# Input: geneList is the list of gene names
# Output: List of all total pairs, as well as correlation from gene PCC data
def getGenePCC(geneList):
    def getPCCfromData(pair):
        return genePCCMatrix.loc[pair[0]][pair[1]]
    pairs = [(gene1, gene2) for i, gene1 in enumerate(geneList) for j, gene2 in enumerate(geneList) if i < j]
    correlations = Parallel(n_jobs=-1)(delayed(getPCCfromData)(pair) for pair in pairs)
    return pairs, correlations

# Input: geneList is the list of geneset names
# Output: List of all total pairs, as well as mean JSD data for each pair
def getMeanJSD(geneList):
    def getMeanfromData(pair):
        return meanJSD_CM.loc[pair[0]][pair[1]]
    pairs = [(gene1, gene2) for i, gene1 in enumerate(geneList) for j, gene2 in enumerate(geneList) if i < j]
    correlations = Parallel(n_jobs=-1)(delayed(getMeanfromData)(pair) for pair in pairs)
    return pairs, correlations

# Input: geneList is the list of gene or geneset names
# Input: Correlations, the output of the functions above
def createCorrelationMatrix(geneList, correlations):
    correlation_matrix = np.zeros((len(geneList), len(geneList)))
    for i, _ in enumerate(geneList):
        for j, _ in enumerate(geneList):
            if i < j:
                correlation_matrix[i, j] = correlations.pop(0)
            else:
                correlation_matrix[i, j] = correlation_matrix[j, i]

    print("Correlation Matrix:")
    plt.figure(figsize=(10, 8))
    plt.imshow(correlation_matrix, cmap='coolwarm', interpolation='nearest')
    plt.colorbar(label='Correlation')
    plt.title('Correlation Matrix')
    plt.tight_layout()
    plt.show()

# Displays the network given by the createNetwork function, do not call alone.
# Outputs network.html, a interactive version of the network G.
def displayNetwork(geneList, G, node_sizes, colors, spread):
    # Graph layout variables
    # Generate the spring layout
    pos = nx.spring_layout(G, scale=2000, k=spread/np.sqrt(G.order()))

    graph_renderer = from_networkx(G, pos, scale=2, center=(0, 0))
    graph_renderer.node_renderer.data_source.data['sizes'] = node_sizes # Original Sizes
    graph_renderer.node_renderer.data_source.data['colors'] = colors # Original Colors
    graph_renderer.node_renderer.glyph = Circle( # Layout for each node
        radius=10,
        fill_color='colors'
    )
    # Define edge colors based on correlation
    edge_colors = ['grey' if data['weight'] > 0 else 'grey' for _, _, data in G.edges(data=True)]
    graph_renderer.edge_renderer.data_source.data['line_colors'] = edge_colors

    graph_renderer.edge_renderer.glyph = MultiLine(line_color='line_colors', line_alpha=0.1, line_width=0.3)

    plot = figure(title = 'Geneset Network', x_range=(-2000, 2000), y_range=(-2000, 2000))

    # Define labels before using it in ColumnDataSource
    labels = {gene: gene for gene in geneList}

    # Define source before using it in LabelSet
    source = ColumnDataSource({'x': np.array([pos[gene][0] for gene in geneList]),
                            'y': np.array([pos[gene][1] for gene in geneList]),
                            'name': [labels[gene] for gene in geneList]})
    
    # Labels appearing/disappearing javascript and code
    labels = LabelSet(x='x', y='y', text='name', source=source,
                background_fill_color='white', text_font_size='8pt')
    labels.visible = False # Set as false so it becomes true once its zoomed in
    zoom_cb = CustomJS(args=dict(labels = labels, plot = plot), code="""
            var xr = [plot.x_range.start, plot.x_range.end]
            var yr = [plot.y_range.start, plot.y_range.end]

            if(Math.abs(xr[1] - xr[0]) < 300 || Math.abs(yr[1] - yr[0]) < 300) { // if zoomed in, show the labels
            labels.visible = true;
            }
            else {labels.visible = false;}
        """)
    plot.js_on_event(events.MouseWheel, zoom_cb)
    plot.add_layout(labels)
    
    # Disease selector javascript and code
    s = Select(title="Select Disease", value = '', options=list(diseases.keys()))
    s.js_on_change('value', CustomJS(args=dict(s=s, plot=plot, diseases=diseases, colors=colors), code=""" 
        // Takes the disease dictionary and changes the colors based on geneset regulation for a disease
        var disease = s.value;
        var genesets = diseases[disease];
        var nodes = plot.renderers[0].node_renderer.data_source.data;
        
        for (var i = 0; i < nodes['index'].length; i++) {
            var gene = nodes['index'][i];
            if (genesets.hasOwnProperty(gene) == true) {
                var regulation = genesets[gene];
                if (regulation == -1) {
                    nodes['colors'][i] = 'lightblue';
                } else if (regulation == 0) {
                    nodes['colors'][i] = 'lime';
                } else if (regulation == 1) {
                    nodes['colors'][i] = 'orangered';
                }
            }
            else {
                nodes['colors'][i] = 'gray';
            }
        }
        plot.renderers[0].node_renderer.data_source.change.emit();
    """))
    plot.renderers.append(graph_renderer)
    # Add the HoverTool to the plot
    ht = HoverTool(tooltips=[("Gene", "@index"), ("Connectivity", "@sizes")])
    plot.add_tools(ht)
    # Output final network
    output_file('network.html')
    show(row(s, plot))

# Threshold is the percent correlation needed to create an edge between two nodes
# geneList is the list of gene or geneset names inputted to be put into the network
# Pairs and Correlations are returned by any of the correlation functions above
# Outputs a network G that is displayed using the displayNetwork() function
def createNetwork(geneList, threshold, pairs, correlations):
    #Creates a networkx graph, addes nodes and edges based of a threshold
    G = nx.Graph()
    for gene in geneList:
        G.add_node(gene)
    for pair, correlation in zip(pairs, correlations):
        if abs(correlation) > threshold:
            G.add_edge(pair[0], pair[1], weight=abs(correlation))

    # Gets the top 10 connected genes
    degrees = dict(G.degree())
    most_connected_genes = sorted(degrees, key=degrees.get, reverse=True)[:100]
    print("\nThe most connected gene sets are:")
    for gene in most_connected_genes:
        print(f"{gene}, Connectivity: {degrees[gene]}")

    # Sets colors
    node_sizes = [(1 * degrees[gene]) for gene in geneList]
    colors = ['red' if gene in most_connected_genes else 'blue' for gene in geneList]
    
    # Display the network
    displayNetwork(geneList, G, node_sizes, colors, 20)
    return G

# Threshold is the mean JSD needed to create an edge between two nodes
# geneList is the list of geneset names inputted to be put into the network
# Pairs and Correlations are returned by any of the correlation functions above
# Outputs a network G that is displayed using the displayNetwork() function
def createJSDNetwork(geneList, threshold, pairs, correlations):
    #Creates a networkx graph, addes nodes and edges based of a threshold
    G = nx.Graph()
    for gene in geneList:
        G.add_node(gene)
    for pair, correlation in zip(pairs, correlations):
        if correlation < threshold and pvalJSD_CM[pair[0]][pair[1]] > (3.127*(10**-7)):
            G.add_edge(pair[0], pair[1], weight=abs(correlation))

    # Gets the top 10 connected genes
    degrees = dict(G.degree())
    most_connected_genes = sorted(degrees, key=degrees.get, reverse=True)[:100]
    print("\nThe most connected gene sets are:")
    for gene in most_connected_genes:
        print(f"{gene}, Connectivity: {degrees[gene]}")

    # Sets colors
    node_sizes = [(1 * degrees[gene]) for gene in geneList]
    colors = ['red' if gene in most_connected_genes else 'blue' for gene in geneList]
    
    # Display the network
    displayNetwork(geneList, G, node_sizes, colors, 8)
    return G

# Input: HOA_data and a network graph
# Output: Edge sums and HOA sums for each geneset
def calculate_hoa_stats(hoa_data, G):
    hoa_gene_set_sum = {}
    for hoa in hoa_data.columns:
        hoa_gene_set_sum[hoa] = sum(G.degree(gene_set) for gene_set in G.nodes() if hoa in hoa_data.loc[gene_set])

    top_nodes = sorted(G.nodes(), key=lambda node: G.degree(node), reverse=True)[:]
    top_nodes_sum = sum(hoa_data.loc[node] for node in top_nodes)

    top_nodes_edge_hoa_sum = sum(hoa_data.loc[node] * G.degree(node) for node in top_nodes)

    return hoa_gene_set_sum, top_nodes_sum, top_nodes_edge_hoa_sum

# Input: HOA stats from calculate_hoa_stats(), text file name
# Output: hoa_component.txt, gene_normalized_hoa_component.txt, and adjusted_normalized_component.txt
def write_hoa_stats(hoa_stats, name):
    hoa_sum, top_nodes_sum, top_nodes_edge_hoa_sum = hoa_stats
    with open(os.path.join('hoa_results', f'{name}.txt'), 'w') as f:
        f.write("Sum of HOA components for each edge in the network:\n")
        for key, value in hoa_sum.items():
            f.write(f"{key}: {value}\n")
        f.write("\nSum of the row values for the top 100 nodes:\n")
        f.write(str(top_nodes_sum))
        f.write("\nEdge-HOA calculation for the top 100 nodes:\n")
        f.write(str(top_nodes_edge_hoa_sum))

# Sample Usage - Get network from Geneset JCC data
# geneset_sorted = geneset_data.sort_values(by='coefficient_estimate', ascending=False).head(1000)
# geneList = geneset_sorted.index.tolist()
# pairs, correlations = getGeneCorrelation_JCC(geneList)
# threshold = np.percentile(correlations, 90)
# createNetwork(geneList, threshold, pairs, correlations)

# Sample Usage - Finding Network from geneset_pcc_data
# geneList = geneset_pcc_data.index.tolist()
# pairs, correlations = getGenesetPCC(geneList)
# createNetwork(geneList, 0.8, pairs, correlations)

# Sample Usage - Finding Network from gene PCC data
# geneList = genePCCMatrix.index.tolist()
# pairs, correlations = getGenePCC(geneList)
# threshold = np.percentile(correlations, 90)
# createNetwork(geneList, threshold, pairs, correlations)

# Sample Usage - Finding Network from JSD data
# geneList = meanJSD_CM.index.tolist()
# pairs, correlations = getMeanJSD(geneList)
# createJSDNetwork(geneList, 0.2, pairs, correlations)

# Sample Usage - Using geneset PCC data, finds all connections for a certain geneset
# describeGeneset('YAGUE_PRETUMOR_DRUG_RESISTANCE_UP', 0.8)

# FINAL USAGE - Network from Geneset PCC data, verifies geneset is located within HOA data.
# Outputs the network to network.html
geneList = hoa_component.index.tolist()
pairs, correlations = getGenesetPCC(geneList)
G = createNetwork(geneList, .825, pairs, correlations)

# FINAL USAGE - Creates correlation matrix of final HOA genesets
createCorrelationMatrix(geneList, correlations)

# FINAL USAGE - Calculates Hallmarks of aging stats based off the network, outputs into...
# hoa_component.txt, gene_normalized_hoa_component.txt, and adjusted_normalized_component.txt
hoa_component_stats = calculate_hoa_stats(hoa_component, G)
gene_normalized_hoa_component_stats = calculate_hoa_stats(gene_normalized_hoa_component, G)
adjusted_hoa_component_stats = calculate_hoa_stats(adjusted_hoa_component, G)
write_hoa_stats(hoa_component_stats, "hoa_component")
write_hoa_stats(gene_normalized_hoa_component_stats, "gene_normalized_hoa_component")
write_hoa_stats(adjusted_hoa_component_stats, "adjusted_hoa_component")
