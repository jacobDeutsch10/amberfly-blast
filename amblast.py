from PyInquirer import style_from_dict, Token, prompt, Separator
from pprint import pprint
from cli_help import questions, DirectoryValidator
from GenomicsRAD import MultiFrame
import os
import subprocess
import pandas as pd
from Bio import SearchIO
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from blastdb import create_db_from_csv, blast_all
from blast_postprocess import *
from datetime import datetime

answers = prompt(questions)
pprint(answers)

try:
    os.mkdir(os.path.join(os.getcwd(), "RAD_sequences"))
except Exception:
    pass
RAD_path =  "./RAD_sequences/" + answers['RADOUT']
multi = MultiFrame.MultiFrame(behaviors=answers['behaviors'])
multi.filename = answers['PATMOS']
multi.get_keys()
multi.read_from_xl_FULL(answers['PATMOS'])
multi.create_behavior_bins()
multi.avg_over_time_step(0.1)
multi.print_multi()
multi.create_atom_codes(num=5)
multi.assign_atom_codes()
multi.convert_frames_to_rad(RAD_path)
del multi

print("RAD sequences generated @ " + RAD_path)
fastas_dir = os.path.join(os.getcwd(), answers["fastas"])
output_dir = os.path.join(os.getcwd(), answers["outputs"])
db_dir = os.path.join(os.getcwd(), "databases/")
global_fsa = os.path.join(db_dir, answers["db"]+".fsa")
db_path = os.path.join(db_dir, answers["db"])
try:
    os.mkdir(db_dir)
except Exception:
    pass
try:
    os.mkdir(fastas_dir)
except Exception:
    pass
try:
    os.mkdir(output_dir)
except Exception:
    pass

create_db_from_csv(RAD_path, global_fsa, answers["fastas"], answers["db"])
blast_all(os.path.join(os.getcwd(), answers["fastas"]), output_dir, db_path)
scores = get_blast_score_matrix(output_dir)

yes = ['y', 'Y', "yes", "Yes"]

heatmap = input("Do you want to see the blast score heatmap?(Y/n)")
if heatmap in yes:
    show_heatmap(scores)

graph = input("Do you want to see the fruchterman reingold graph?(Y/n)")
if graph in yes:
    try:
        os.mkdir("graphs")
    except Exception:
        pass
    graph_path = os.path.join("./graphs/", "{}-{}.graphml".format("graph",datetime.now().strftime("%d-%m-%H-%M")))
    show_networkx(scores, graph_path)
