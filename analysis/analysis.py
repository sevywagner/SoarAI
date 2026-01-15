from sklearn.metrics import roc_auc_score, f1_score, confusion_matrix, roc_curve
import matplotlib.pyplot as plt
import yaml
import numpy as np
import math
import sys
from pathlib import Path

if (len(sys.argv) != 3):
    print("Invalid arguments. Usage: python3 ../analysis/analysis.py <path to yaml config> <reagent ion (NH4|NO)>")
    sys.exit(1)

with open(sys.argv[1], 'r') as f:
    config = yaml.safe_load(f)

reagent_ion = sys.argv[2]
    
run_dir = config['paths']['core']['run_root'] + config['run']['run_id'] + '/'
preds_dir = run_dir + config['paths']['core']['pred_output_dir']
    
files = [
    preds_dir + reagent_ion + '_xgboost_preds.txt',
    preds_dir + reagent_ion + '_lgbm_preds.txt',
    preds_dir + reagent_ion + '_gbdt_preds.txt',
    preds_dir + reagent_ion + '_nn_preds.txt',
]

ys = []
y_preds = []

for path in files:
    y = []
    y_pred = []
    
    with open(path, 'r') as fr:
        lines = fr.readlines()
    
        for line in lines:
            s = line.split(' ')
            
            y.append(int(s[0]))
            y_pred.append(int(s[1]))

        fr.close()

    ys.append(y)
    y_preds.append(y_pred)
    test = np.array(y) == 1

colors = ['red', 'blue', 'green', 'yellow']
labels = ['XGBoost', 'LightGBM', 'GBDT', 'Neural Network']
labels = [sys.argv[2] + ' ' + x for x in labels]

fig, ax = plt.subplots(1, 1)
cms = []

for i in range(len(ys)):
    print('F1: ', f1_score(ys[i], y_preds[i]))
    print('ROC AUC: ', roc_auc_score(ys[i], y_preds[i]))

    accuracy = 0
    for y, y_pred in zip(ys[i], y_preds[i]):
        if y == y_pred:
            accuracy += 1
    
    print('Accuracy: ', accuracy / len(ys[i]))
    
    cms.append(confusion_matrix(ys[i], y_preds[i]))

    fpr, tpr, _ = roc_curve(ys[i], y_preds[i])
    ax.plot(fpr, tpr, color=colors[i], linestyle='--', label=labels[i])

ax.set_xlabel('False Positive Rate')
ax.set_ylabel('True Positive Rate')

ax.set_title('ROC Curve')
ax.legend()

plots_dir = run_dir + config['paths']['analysis']['plot_dir']
Path(plots_dir).mkdir(parents=True, exist_ok=True)
fig.savefig(f'{plots_dir}{sys.argv[2]}_ROC_AUC.png', format='png', dpi=600)

plt.show()

for x, cm in enumerate(cms):
    fig, ax = plt.subplots(1, 1)
    im = ax.imshow(cm, cmap="Blues")
    plt.colorbar(im, ax=ax)

    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            ax.text(j, i, cm[i, j], ha='center', va='center', color=('white' if i == 0 and j == 0 else 'black'))

    ax.set_xlabel('Predicted Label')
    ax.set_ylabel('True Label')
    ax.set_title(f'{labels[x]} Confusion Matrix')
    ax.set_xticks(np.arange(cm.shape[1]))
    ax.set_yticks(np.arange(cm.shape[0]))
    ax.set_xticklabels(['0', '1'])
    ax.set_yticklabels(['0', '1'])

    fig.savefig(f'{plots_dir}{labels[x]}.png', format='png', dpi=600)
    
    plt.show()
