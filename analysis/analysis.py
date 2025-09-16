from sklearn.metrics import roc_auc_score, f1_score, confusion_matrix, roc_curve
import matplotlib.pyplot as plt
import numpy as np
import sys

if (len(sys.argv) < 3):
    print("Must pass in path(s) preds and output plots directory ie python3 this_program.py /path/to/plots/dir/ /path/to/preds1.txt /path/to/preds2.txt ...")
    sys.exit(1)

if sys.argv[1][-1] != '/':
    sys.argv[1] += '/'

files = sys.argv[2:]

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

colors = ['red', 'blue', 'green', 'yellow']
labels = ['XGBoost', 'LightGBM', 'GBDT', 'Neural Network']
    
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

fig.savefig(f'{sys.argv[1]}ROC_AUC.png', format='png', dpi=600)

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

    fig.savefig(f'{sys.argv[1]}{labels[x]}.png', format='png', dpi=600)
    
    plt.show()
