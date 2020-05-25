### amberfly-blast
### Step 1: clone this repo
```
git clone https://github.com/jacobDeutsch10/amberfly-blast.git
```
### Step 2: update the submodule for GenomicsRAD. Type the following commands in the top directory for amberfly-blast
```
git submodule init
git submodule update
```
### Step 3:Install blast
```
sudo apt-get install ncbi-blast+
````
### Step 4: install dependencies with pip:
- numpy
- scipy
- matplotlib
- pandas
- networkx
- biopython
- seaborn
- PyInquirer

### Step 5: Run the program:
```
python amblast.py
```
### Step 6: Follow prompts. Make sure to type full path of PATMOS output file.
