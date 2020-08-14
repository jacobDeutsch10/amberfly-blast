import pandas as pd
import os
import subprocess
import statistics
# makeblastdb -in ./amberfly.fsa -dbtype nucl -out ./amberflydb

def create_db_from_csv(csv_path, fsa_db_path, fastas_dir, db_path, atom_length=5, behaviors=4):
    """
    csv_path: path to csv file containing RAD sequences
    fsa_db_path: path to global fsa file used to generate the database
    fastas_dir: directory to place individual fasta files in
    db_path: path to db
    """
    df = pd.read_csv(csv_path, header=0)
    lengths = []
    count = 0
    nucs_per_seq = atom_length*(behaviors+2)
    dbfile = open(fsa_db_path, "w")
    stats_file = open("seq_stats.txt", "w")

    for i, row in df.iterrows():
        start = 0.0
        key = row["keys"]
        seq = row["Strings"]
        fasta_filepath = os.path.join(fastas_dir, str(count) + "_" + key + ".fsa")
        seq_len = len(seq)
        if seq_len < 9*nucs_per_seq:
            continue
        elif seq_len < 20*nucs_per_seq:
            end = len(seq)/nucs_per_seq * 0.1
            times = "_{:0.1f}_{:0.1f}".format(start, end)
            ofile = open(fasta_filepath, "w" )
            ofile.write(">" + str(count) + times + "_" + key + "\n" + seq + "\n")
            dbfile.write(">" + str(count) + times + "_" + key + "\n" + seq + "\n")
            stats_file.write(str(count) + times + "_" + key + ": "+ str(len(seq))+"\n")
            ofile.close()
            count +=1
        elif seq_len < 30*nucs_per_seq:
            mid = seq_len//2
            mod = mid % nucs_per_seq
            mid = mid + mod
            end = (mid / nucs_per_seq) * 0.1
            times = "_{:0.1f}_{:0.1f}".format(start, end)
            ofile = open(fasta_filepath, "w" )
            ofile.write(">" + str(count) + times + "_" + key + "\n" + seq[0:mid] + "\n")
            dbfile.write(">" + str(count) + times + "_" + key + "\n" + seq[0:mid] + "\n")
            stats_file.write(str(count) + times + "_" + key + ": "+ str(len(seq[0:mid]))+"\n")
            ofile.close()
            count += 1
            start = end
            end = ((seq_len - mid)/ nucs_per_seq) * 0.1 + start
            times = "_{:0.1f}_{:0.1f}".format(start, end)
            fasta_filepath = os.path.join(fastas_dir, str(count) + "_" + key + ".fsa")
            ofile = open(fasta_filepath, "w" )
            ofile.write(">" + str(count) + times + "_" + key + "\n" + seq[mid:] + "\n")
            dbfile.write(">" + str(count) + times + "_" + key + "\n" + seq[mid:] + "\n")
            stats_file.write(str(count) + times + "_" + key + ":" + str(len(seq[mid:]))+"\n")
            ofile.close()
            count+=1
        else:
            while seq_len > 20*nucs_per_seq:
                temp = seq[0:15*nucs_per_seq]
                seq = seq[15*nucs_per_seq:]
                seq_len = len(seq)
                end = 15*0.1 + start
                times = "_{:0.1f}_{:0.1f}".format(start, end)
                fasta_filepath = os.path.join(fastas_dir, str(count) + "_" + key + ".fsa")
                ofile = open(fasta_filepath, "w" )
                ofile.write(">" + str(count) + times + "_" + key + "\n" + temp + "\n")
                dbfile.write(">" + str(count) + times +"_" + key + "\n" + temp + "\n")
                stats_file.write(str(count) + times +"_" + key + ": " + str(len(temp))+"\n")
                ofile.close()
                count +=1
                start = end
            if seq_len > 8*nucs_per_seq:
                end = (seq_len/nucs_per_seq)*0.1 + start
                fasta_filepath = os.path.join(fastas_dir, str(count) + "_" + key + ".fsa")
                times = "_{:0.1f}_{:0.1f}".format(start, end)
                ofile = open(fasta_filepath, "w" )
                ofile.write(">" + str(count) +times + "_" + key + "\n" + seq + "\n")
                dbfile.write(">" + str(count) + times +"_" + key + "\n" + seq + "\n")
                stats_file.write(str(count) + times +"_" + key + ": " + str(len(seq))+"\n")
                ofile.close()
                count +=1 
        lengths.append(len(row["Strings"]))
        print(key + ": " + str(len(row["Strings"]))) 
    dbfile.close()
    stats_file.close()
    subprocess.call(['makeblastdb', '-in', fsa_db_path, '-dbtype', 'nucl', '-out', db_path])
    print("average: " + str(statistics.mean(lengths)))
    print("meadian: " + str(statistics.median(lengths)))
    print("min: " + str(min(lengths)))
    print("max: " + str(max(lengths)))

def blast_all(fastas_dir, output_dir, db_path):
    for f in os.listdir(fastas_dir):
        key_name = f.split(".")[0]
        fasta_filepath = os.path.join(fastas_dir, f)
        output_filepath = os.path.join(output_dir, "out_"+ key_name + ".txt")
        subprocess.call(['blastn', '-query', fasta_filepath, '-db', db_path, '-outfmt', '6','-out', output_filepath])
