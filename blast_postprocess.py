import pandas as pd
import os
from Bio import SearchIO
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx

def get_blast_score_matrix(blast_dir):
    files = os.listdir(blast_dir)
    num_files =len(files)
    scores = np.zeros((num_files,num_files), dtype=int)
    for i,f in enumerate(files):
        f_path = os.path.join(blast_dir, f)
        if os.stat(f_path).st_size != 0:
            result = SearchIO.read(f_path, 'blast-tab')
            query_i = int(result.id.split("_")[0]) - 1
            for hit in result:
                seq_i = int(hit.id.split("_")[0]) - 1
                if seq_i != query_i:
                    scores[query_i][seq_i] = int(hit[0].bitscore)
    return scores

def score_matrix_to_file(scores, output_path):
    df = pd.DataFrame(data=scores,
          index=np.arange(1, len(scores) + 1),
          columns=np.arange(1, len(scores) + 1))
    df.to_csv(output_path)

def show_heatmap(scores):
    fig, ax = plt.subplots(figsize=(100,100))
    title = "Blast Bitscore Heatmap"
    plt.title(title, fontsize=12)
    ttl = ax.title
    ttl.set_position([0.5, 1.05])
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis('off')
    hm = sns.heatmap(scores, annot = True, fmt="", cmap='YlGnBu',annot_kws={"size": 5}, ax=ax, vmin = 0, vmax =1100)
    plt.show()

def show_networkx(scores, output_path = "", show_plot=True):
    G = nx.Graph()
    for i in range(0,len(scores)):
        G.add_node(i)
    for i in range(0, len(scores)-1):
        for j in range(i, len(scores[0])):
            if i != j and scores[i][j] != 0:
                G.add_edge(i,j, weight=scores[i][j])
                print("({}, {}): {}".format(i, j, scores[i,j]))

    pos = nx.spring_layout(G)
    print(pos)
    for n,p in pos.items():
        G.nodes[n]['x'] = p[0]
        G.nodes[n]['y'] = p[1]
    if output_path.endswith(".graphml"):
         nx.write_graphml(G, output_path)
    nx.draw_networkx(G, pos, with_labels=True)
    if show_plot:
        plt.show()
    return G
