import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import math
import sys
import yaml
from pathlib import Path
import pandas as pd

if (len(sys.argv) != 2):
    print("Invalid arguments. Usage: python3 ../analysis/van_krevelen.py <path to yaml config>")
    sys.exit(1)

with open(sys.argv[1], 'r') as f:
    config = yaml.safe_load(f)


chems = { 'H3O': 0, 'NH4': 1, 'NH3': 2, 'H2O': 3, 'C': 4, 'H': 5, 'O': 6, 'N': 7, 'S': 8, 'F': 9, 'Si': 10, 'Cl': 11, 'K': 12, 'Ar': 13, 'Ca': 14, 'Na': 15, 'Al': 16, 'Br': 17, 'P': 18, '+': 19  }
chem_masses = { 'H3O': 19.01838971, 'NH4': 18.03437412, 'NH3': 17.02654909, 'H2O': 18.010564679999998, 'C': 12, 'H': 1.00782503, 'O': 15.99491462, 'N': 14.00307400, 'S': 31.97207117, 'F': 18.998403, 'Si': 28.0855, 'Cl': 35.45, 'K': 38.96370649, 'Ar': 39.96238312, 'Ca': 39.96259086, 'Na': 22.98976928, 'Al': 26.98153853, 'Br': 78.9183376, 'P': 30.973761998, '+': -0.000548 }
n_elements = len(chems)

def encodeCompounds(chemicals):
    '''
        Function:
            Encode a list chemical compounds
        
        Parameters:
            chemicals (list(str)): Unencoded compound strings

        Returns:
            df (pd.DataFrame): Encoded compounds in a dataframe
            encoded_compounds (list(int)): Encoded compounds in a list
    '''
    encoded_compounds = []
    zero = [0 for _ in range(len(chems))]

    ctr = 1
    for chemical in chemicals:
        ctr += 1
        encoding = list(zero)
        while '(' in chemical:
            end = chemical.index(")")
            a = chemical[chemical.index("(") + 1 : end]
            encoding[chems[a]] = .01
            chemical = chemical[end + 1:]
        
        for i in range(len(chemical)):
            c = chemical[i]
            if i < len(chemical) - 1 and chemical[i + 1].isalpha() and chemical[i + 1].islower():
                c += chemical[i + 1]

            if (not (c[0].isalpha() and c[0].isupper())):
                continue

            num = ""
            for j in range(i + len(c), len(chemical)):
                if (not chemical[j].isdigit()):
                    break

                num += chemical[j]

            if num == '':
                num = 1
            
            encoding[chems[c]] = float(num) / float(100)
                
        if encoding != list(zero):
            encoded_compounds.append(encoding)
        else:
            print(ctr)
            print(chemical)
            
    df = pd.DataFrame()
    ctr = 0
    for key in chems:
        df.insert(loc=ctr, column=key, value=[i[ctr] for i in encoded_compounds])
        ctr += 1
        
    return df, encoded_compounds

    
oc = []
hc = []
colors = []

area1 = ""
area1_line = 0
area2 = ""
area2_line = 0

run_dir = config['paths']['core']['run_root'] + config['run']['run_id'] + '/'
plots_dir = run_dir + config['paths']['analysis']['plot_dir']
Path(plots_dir).mkdir(parents=True, exist_ok=True)
peak_list_path = config['paths']['core']['in_peak_list_dir'] + config['paths']['core']['in_peak_list_filename']

with open(peak_list_path, 'r') as f:
    lines = f.readlines()[1:]

    for (idx, i) in enumerate(lines):
        if i.split('\t')[5] == '':
            continue
            
        _, ec = encodeCompounds([i.split('\t')[5]])

        if (ec[0][4] == 0):
            continue

        sample_oc = ec[0][6] / ec[0][4]
        sample_hc = ec[0][5] / ec[0][4]
    
        oc.append(ec[0][6] / ec[0][4])
        hc.append(ec[0][5] / ec[0][4])
        
        if area1 == "" and sample_oc <= .3 and sample_hc >= 1.5 and sample_hc <= 2.5:
            area1 = i.split('\t')[5]
            area1_line = idx

        if area2 == "" and sample_oc > .3 and sample_oc <= .75 and sample_hc >= 1 and sample_hc <= 2:
            area2 = i.split('\t')[5]
            area2_line = idx
    
        colors.append(math.floor(float(i.split('\t')[3])))
        
    f.close()

fig = plt.figure(figsize=(8, 5))
gs = fig.add_gridspec(
    1, 2,
    width_ratios=[1, 0.035],
    wspace=0.08
)

ax = fig.add_subplot(gs[0])
cax = fig.add_subplot(gs[1])

a = ax.scatter(oc, hc, cmap='viridis', alpha=0.7, c=colors, s=[i / 100 for i in colors])

ax.set_xlim(-.1, 3)   # O:C typical range
ax.set_ylim(-.1, 4)   # H:C typical range

#plt.gca().add_patch(Rectangle((0.2, 0.7), 0.3, 0.5,
#                              edgecolor='red', facecolor='none', lw=1.5, linestyle='--'))
#plt.text(0.22, 1.15, 'Lignin-like', fontsize=10, color='red')
ax.add_patch(plt.Rectangle(
    (0, 1.5),  # (O:C start, H:C start)
    0.4,         # O:C width
    1.0,         # H:C height
    edgecolor='orange',
    facecolor='none',
    lw=4,
    label='Primary VOC emissions'
))


ax.add_patch(plt.Rectangle(
    (.4, 1.0),  # (O:C start, H:C start)
    0.4,         # O:C width
    1.0,         # H:C height
    edgecolor='yellow',
    facecolor='none',
    lw=4,
    label='Oxidized gas-phase products'
))

ax.legend(loc='upper right', fontsize=17, frameon=True)

ax.set_title("Van Krevelen plot", fontsize=27)
ax.set_xlabel("O:C", fontsize=25)
ax.set_ylabel("H:C", fontsize=25)
ax.tick_params(axis='both', labelsize=20)

cb = fig.colorbar(a, cax=cax)
cb.set_label('m/z', fontsize=20)
cb.ax.tick_params(labelsize=20)

fig.subplots_adjust(
    left=0.18,
    right=0.88,
    bottom=0.17,
    top=0.88
)

fig.savefig(f'{plots_dir}van_krevelen.pdf', format='pdf', dpi=1000)
plt.show()
