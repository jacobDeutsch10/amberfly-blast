import pandas as pd
import os
from Bio import SearchIO
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from mpl_toolkits.mplot3d.axes3d import Axes3D

def get_blast_score_matrix(blast_dir):
    files = os.listdir(blast_dir)
    num_files =len(files)
    scores = np.zeros((num_files,num_files), dtype=int)
    for i,f in enumerate(files):
        f_path = os.path.join(blast_dir, f)
        if os.stat(f_path).st_size != 0:
            result = SearchIO.read(f_path, 'blast-tab')
            query_i = int(result.id.split("_")[0])
            for hit in result:
                seq_i = int(hit.id.split("_")[0])
                if seq_i != query_i and int(hit[0].bitscore) > scores[query_i][seq_i] :
                    scores[query_i][seq_i] = int(hit[0].bitscore)
    return scores

def score_matrix_to_file(scores, output_path):
    df = pd.DataFrame(data=scores,
          index=np.arange(1, len(scores) + 1),
          columns=np.arange(1, len(scores) + 1))
    df.to_csv(output_path)

def read_score_matrix(input_path):
    df = pd.read_csv(input_path, header=0, index_col=0)
    return df.to_numpy(dtype=int)

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

def show_networkx(scores, output_path = "", show_plot=True, threshold=400):
    G = nx.Graph()
    for i in range(0,len(scores)):
        G.add_node(i)
    for i in range(0, len(scores)-1):
        for j in range(i, len(scores[0])):
            if i != j and scores[i][j] != 0 and scores[i][j] > threshold:
                G.add_edge(i,j, weight=(scores[i][j]**2))
                print("({}, {}): {}".format(i, j, scores[i,j]))

    pos = nx.spring_layout(G,weight='weight', scale=50, iterations = 500)
    for n,p in pos.items():
        G.nodes[n]['x'] = p[0]
        G.nodes[n]['y'] = p[1]
    if output_path.endswith(".graphml"):
         nx.write_graphml(G, output_path)
    nx.draw_networkx(G, pos, with_labels=True)
    if show_plot:
        plt.show()
    else:
        plt.clf()
    return G
def plot_with_key(data_dir, key, ax, ax_v):
    
    nums = key.split("_")
    start = float(nums[1])
    end = float(nums[2])
    csv_name = "_".join(nums[3:])+ ".csv"
    df = pd.read_csv(os.path.join(data_dir, csv_name))
    offset = df['dt'][0]
    df['dt'] = df['dt'].apply(lambda x: x - offset)
    start_ind = 0
    i = 0
    while df['dt'][i] < start:
        start_ind +=1
        i += 1
    end_ind = start_ind

    while i < len(df['dt']) and end >= df['dt'][i] :
        end_ind += 1
        i += 1
    cols = df.columns
    
    """print(key, end=":")
    print("({}, {})".format(start_ind, end_ind))
    print(df.head(15))
    print(df.tail(15))
    print(key, end=":")
    print("({}, {})".format(start_ind, end_ind))
    print(df[cols[1]][start_ind:end_ind])
    print(key, end=":")
    print("({}, {})".format(start_ind, end_ind))
    print(df[cols[2]][start_ind:end_ind])
    print(key, end=":")
    print("({}, {})".format(start_ind, end_ind))
    print(df[cols[3]][start_ind:end_ind])"""
    ax.scatter3D(df[cols[5]][start_ind:end_ind], df[cols[6]][start_ind:end_ind], df[cols[13]][start_ind:end_ind], cmap=df[cols[23]][start_ind:end_ind])
    ax.scatter3D([df[cols[5]][start_ind]], [df[cols[6]][start_ind]], [df[cols[13]][start_ind]], c='r')
    #ax.scatter3D([df[cols[5]][end_ind-1]], [df[cols[6]][end_ind-1]], [df[cols[13]][end_ind-1]], c='r')
    ax.set_title(key)
    ax.set_xlabel(cols[5])
    ax.set_ylabel(cols[6])
    ax.set_zlabel(cols[13])
    ax_v.scatter(df[cols[1]][start_ind:end_ind], df[cols[23]][start_ind:end_ind])
    ax_v.set_xlabel(cols[1])
    ax_v.set_ylabel(cols[23])
