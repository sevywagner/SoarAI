from sklearn.metrics import roc_auc_score, precision_score, recall_score, confusion_matrix, roc_curve
import matplotlib.pyplot as plt
import yaml
import numpy as np
import math
import sys
from pathlib import Path

if (len(sys.argv) != 3):
    print("Invalid arguments. Usage: python3 ../analysis/analysis.py <path to yaml config> <logit dump id>")
    sys.exit(1)

with open(sys.argv[1], 'r') as f:
    config = yaml.safe_load(f)

dump_logits = sys.argv[2]

reagent_ions = [
    'NH4',
    'NO',
    'H3O'
]

run_dir = config['paths']['core']['run_root'] + config['run']['run_id'] + '/'
preds_dir = run_dir + config['paths']['core']['pred_output_dir']
logits_dir = preds_dir + config['paths']['core']['pred_logit_output_dir']

logit_files = [[logits_dir + 'xgboost_logits_' + dump_logits + '.txt', logits_dir + 'nn_logits_' + dump_logits + '.txt']]
files = [[preds_dir + 'xgboost_preds.txt', preds_dir + 'nn_preds.txt']]

for ion in reagent_ions:
    logit_files.append([
        logits_dir + ion + '_xgboost_logits_' + dump_logits + '.txt',
        logits_dir + ion + '_nn_logits_' + dump_logits + '.txt',
    ])

    files.append([
        preds_dir + ion + '_xgboost_preds.txt',
        preds_dir + ion + '_nn_preds.txt',
    ])


for idx, (reagent_files, reagent_logit_files) in enumerate(zip(files, logit_files)):
    ys = []
    y_preds = []
    y_preds_raw = []

    print('========================================')
    print('All reagents' if idx == 0 else reagent_ions[idx - 1])
    print('========================================')

    for path in reagent_files:
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

    for path in reagent_logit_files:
        y_pred_raw = []

        with open(path, 'r') as fr:
            lines = fr.readlines()

            for line in lines:
                y_pred_raw.append(float(line))

            fr.close()

            y_preds_raw.append(y_pred_raw)

        colors = ['red', 'blue']
        labels = []
        if idx == 0:
            labels = ['XGBoost', 'NN']
        else:
            labels = ['XGBoost_' + reagent_ions[idx - 1], 'NN_' + reagent_ions[idx - 1]]

        fig, ax = plt.subplots(1, 1)
        cms = []

    for i in range(len(ys)):
        print(labels[i])
        print('----------------------------------------')
        print('Precision: ', precision_score(ys[i], y_preds[i]))
        print('Recall: ', recall_score(ys[i], y_preds[i]))
        print('ROC AUC: ', roc_auc_score(ys[i], y_preds_raw[i]))

        accuracy = 0
        for y, y_pred in zip(ys[i], y_preds[i]):
            if y == y_pred:
                accuracy += 1

        print('Accuracy: ', accuracy / len(ys[i]))
        print('----------------------------------------')

        cms.append(confusion_matrix(ys[i], y_preds[i]))

        fpr, tpr, _ = roc_curve(ys[i], y_preds_raw[i])
        ax.plot(fpr, tpr, color=colors[i], linestyle='--', label=labels[i])

    ax.set_xlabel('False positive rate')
    ax.set_ylabel('True positive rate')

    ax.set_title('ROC curve')
    ax.legend()

    plots_dir = run_dir + config['paths']['analysis']['plot_dir']
    Path(plots_dir).mkdir(parents=True, exist_ok=True)
    fig.savefig(f"{plots_dir}{reagent_ions[idx - 1] + '_' if idx > 0 else ''}ROC_AUC.png", format='png', dpi=600)

    plt.show()

    for x, cm in enumerate(cms):
        fig, ax = plt.subplots(1, 1, figsize=(7, 6))
        im = ax.imshow(cm, cmap="Blues")
        cb = plt.colorbar(im, ax=ax, pad=0.08)
        cb.ax.tick_params(labelsize=17)

        cb.ax.yaxis.get_offset_text().set_visible(False)

        fig.canvas.draw()

        if cb.formatter.orderOfMagnitude == 6:
            cb.ax.yaxis.get_offset_text().set_visible(False)
            cb.ax.text(
                0.5, 1.02, r'$10^6$',
                transform=cb.ax.transAxes,
                ha='center', va='bottom', fontsize=20
            )
        else:
            cb.ax.yaxis.get_offset_text().set_visible(True)

        for i in range(cm.shape[0]):
            for j in range(cm.shape[1]):
                ax.text(j, i, cm[i, j], ha='center', va='center', fontsize=25, color=('white' if i == 0 and j == 0 else 'black'))

        ax.set_xlabel('Predicted label', fontsize=25)
        ax.set_ylabel('True label', fontsize=25)
        ax.set_title(f'{labels[x]} confusion matrix', fontsize=27)
        ax.set_xticks(np.arange(cm.shape[1]))
        ax.set_yticks(np.arange(cm.shape[0]))
        ax.set_xticklabels(['0', '1'])
        ax.set_yticklabels(['0', '1'])
        ax.tick_params(axis='both', labelsize=17)

        fig.savefig(f'{plots_dir}{labels[x]}.pdf', pad_inches=.2, format='pdf', dpi=600)

        plt.show()
