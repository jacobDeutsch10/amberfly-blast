from PyInquirer import style_from_dict, Token, prompt, Separator
from pprint import pprint
from cli_help import *
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
import math
yes = ['y', 'Y', "yes", "Yes"]
init_answers = prompt(init_question)
pprint(init_answers)
G = None
if init_answers['mode'] == init_choices[0]:
    answers = prompt(questions)
    pprint(answers)

    base_path = os.path.join(os.getcwd(), answers["DB"])
    data_dir = os.path.join(base_path, "data")
    print(base_path)
    try:
        os.mkdir(base_path)
    except Exception:
        pass
    try:
        os.makedirs(data_dir)
    except Exception:
        pass
    atom_length = 5
    time_step = 0.1
    num_behaviors = len(answers['behaviors'])
    RAD_path = os.path.join(base_path, answers['DB']+".csv")
    multi = MultiFrame.MultiFrame(behaviors=answers['behaviors'])
    multi.filename = answers['PATMOS']
    multi.get_keys()
    multi.read_from_xl_FULL(answers['PATMOS'])
    multi.write_csvs(data_dir)
    multi.drop_columns()
    multi.create_behavior_bins()
    multi.avg_over_time_step(time_step)
    multi.print_multi()
    multi.create_atom_codes(num=atom_length)
    multi.assign_atom_codes()
    multi.convert_frames_to_rad(RAD_path)
    del multi
    
    print("RAD sequences generated @ " + RAD_path)
    fastas_dir = os.path.join(base_path, "fastas")
    output_dir = os.path.join(base_path, "blast-outputs")
    db_dir = os.path.join(base_path, "db/")
    global_fsa = os.path.join(db_dir, answers["DB"]+".fsa")
    db_path = os.path.join(db_dir, answers["DB"])
    data_dir = os.path.join(base_path, "data")

    try:
        os.makedirs(db_dir)
    except Exception:
        pass
    try:
        os.makedirs(fastas_dir)
    except Exception:
        pass
    try:
        os.makedirs(output_dir)
    except Exception:
        pass

    
    create_db_from_csv(RAD_path, global_fsa, fastas_dir, db_path, behaviors=num_behaviors)
    blast_all(fastas_dir, output_dir, db_path)
    scores = get_blast_score_matrix(output_dir)
    score_matrix_to_file(scores, os.path.join(base_path, "scores.csv"))
    
    
    
    heatmap = input("Do you want to see the blast score heatmap?(Y/n)")
    if heatmap in yes:
        show_heatmap(scores)
    
    graph = input("Do you want to see the fruchterman reingold graph?(Y/n)")
    if graph in yes:
        try:
            os.mkdir(os.path.join(base_path, "graphs"))
        except Exception:
            pass
        graph_path = os.path.join(answers["DB"],"graphs/", "{}-{}.graphml".format("graph",datetime.now().strftime("%d-%m-%H-%M")))
        thresh = get_numeric_input("Enter bitscore threshold for edges on graph")
        G = show_networkx(scores, graph_path, threshold=thresh)

data_mode = input("Do you want to enter data exploration mode?(Y/n)")
if init_answers['mode'] != init_choices[0]:
    db = prompt({
            'type': 'input',
            'name': 'db',
            'message': 'name of existing db',
            'validate': DirectoryValidator})
    base_path = os.path.join(os.getcwd(), db['db'])
    scores = read_score_matrix(os.path.join(base_path, 'scores.csv'))
    print(scores)

if data_mode in yes:
    x = -1
    y = -1
    if G is None:
        try:
            os.mkdir(os.path.join(base_path, "graphs"))
        except Exception:
            pass
        graph_path = os.path.join(base_path,"graphs/", "{}-{}.graphml".format("graph",datetime.now().strftime("%d-%m-%H-%M")))
        thresh = get_numeric_input("Enter bitscore threshold for edges on graph: ")
        G = show_networkx(scores, graph_path, show_plot=False, threshold=thresh)
    done = False
    while not done:
        msg = "Enter an integer between 0 & {}: ".format(len(scores)-1)
        num = input(msg)
        if not num.isnumeric:
            print("must enter a numerical value")
            continue
        else:
            x = int(num)
            num = input(msg)
        if not num.isnumeric:
            print("must enter a numerical value")
            continue
        else:
            y = int(num)
        if x > (len(scores)-1) or x < 0 or y > (len(scores)-1) or y < 0:
            print("indexes are our of bounds")
            continue
        else:
            score = scores[x][y]
            distance = math.sqrt((G.nodes[x]['x']-G.nodes[y]['x'])**2 + (G.nodes[x]['y']-G.nodes[y]['y'])**2 )
            print("({}, {}): score: {:d} distance: {:0.8f}".format(x, y, score, distance))
            print("node {} @ ({}, {})".format(x, G.nodes[x]['x'], G.nodes[x]['y']))
            print("node {} @ ({}, {})".format(y, G.nodes[y]['x'], G.nodes[y]['y']))
        is_done = input("say yes if you want to exit: ")
        done = is_done in yes
print("amblast exit")
