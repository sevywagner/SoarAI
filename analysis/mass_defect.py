import matplotlib.pyplot as plt
import numpy as np
import math
import sys
import yaml
from pathlib import Path

x = []
y = []
colors = []

if (len(sys.argv) != 2):
    print("Invalid arguments. Usage: python3 ../analysis/mass_defect.py <path to yaml config>")
    sys.exit(1)

with open(sys.argv[1], 'r') as f:
    config = yaml.safe_load(f)
    
run_dir = config['paths']['core']['run_root'] + config['run']['run_id'] + '/'
plots_dir = run_dir + config['paths']['analysis']['plot_dir']
Path(plots_dir).mkdir(parents=True, exist_ok=True)
peak_list_path = config['paths']['core']['in_peak_list_dir'] + config['paths']['core']['in_peak_list_filename']

with open(peak_list_path, 'r') as f:
    lines = f.readlines()[1:]
    for i in range(len(lines)):
        mz = float(lines[i].split('\t')[3])
        x.append(mz)
        y.append(mz - round(mz))
        colors.append(mz)
    f.close()

fig = plt.figure(figsize=(8, 5))
gs = fig.add_gridspec(
    1, 2,
    width_ratios=[1, 0.035],
    wspace=0.08
)

ax = fig.add_subplot(gs[0])
#cax = fig.add_subplot(gs[1])

sc = ax.scatter(x, y, s=2, alpha=0.7)

ax.annotate(
    "",
    xy=(435, 0.43),
    xytext=(270, 0.235),
    arrowprops=dict(
        arrowstyle="->",
        linestyle="-",
        linewidth=3,
        edgecolor="black",
        facecolor="white",
        mutation_scale=18,
    ),
    zorder=5,
)

ax.text(
    340,
    0.285,
    r"$+\mathrm{CH}_2$",
    fontsize=22,
    rotation=32,
    ha="left",
    va="center",
    color="black",
    zorder=6,
)

ax.set_title("Mass defect plot", fontsize=27)
ax.set_xlabel("m/z", fontsize=25)
ax.set_ylabel("Defect", fontsize=25)
ax.tick_params(axis='both', labelsize=20)


#plt.rcParams.update({'font.size': 14, 'axes.labelweight': 'bold'})

fig.subplots_adjust(
    left=0.18,
    right=0.88,
    bottom=0.17,
    top=0.88
)

fig.savefig(f'{plots_dir}mass_defect.pdf', format='pdf', dpi=1000)
plt.show()
