import numpy as np
from sklearn.metrics import recall_score, precision_score, f1_score, confusion_matrix, roc_auc_score
from sklearn.model_selection import train_test_split
import pickle
import argparse
import yaml

parser = argparse.ArgumentParser(
    prog='CreateModels',
    description='Create TF, GB, and/or LGBM models'
)

parser.add_argument('-n', '--neural-net', action='store_true')
parser.add_argument('-g', '--gradient-booster', action='store_true')
parser.add_argument('-l', '--lgbm', action='store_true')
parser.add_argument('-d', '--dump-logits')
parser.add_argument('-r', '--reagent-ion', required=True)
parser.add_argument('-c', '--config', required=True)

args = parser.parse_args()

with open(args.config, 'r') as f:
    config = yaml.safe_load(f)

RANDOM_SEED = config['run']['py_model_train_seed']

def getData(comboFile, split=False):
    '''
        Function:
            Get data from labeled "combo" file and seperate train and test sets
        
        Parameters:
            comboFile (str): name of the combo file
            split (boolean): whether or not to split the data into train or test sets

        Returns:
            if split
                x_train, x_test, y_train, y_test (list(list(float)), list(list(float)), list(int), list(int)):
                    training data samples,
                    test data samples,
                    train labels,
                    test labels
            else
                x, y (list(list(float))), list(int)):
                    data samples,
                    labels
    '''

    X = []
    y = []
    with open(comboFile) as f:
        lines = f.readlines()
        for line in lines:
            segments = line.split(',')

            sample = [float(i) for i in segments[:-1]]
            label = float(segments[-1])

            X.append(sample)
            y.append(label)

        f.close()

    if split:
        return train_test_split(X, y,
                                test_size=.3,
                                random_state=RANDOM_SEED,
                                shuffle=True)
    else:
        return np.array(X, dtype=np.float32), np.array(y, dtype=np.int8)

def output_preds(output_path, preds, labels):
    '''
        Function:
            Output model predictions and true labels side by side to a file

        Parameters:
            output_path (str): path to the file to output the preds to
            preds (list(int) || list(float)): model predictions
            labels (list(int)): true data labels
    '''
    with open(output_path, 'w+') as fw:
        for i in range(len(preds)):
            fw.write(str(labels[i]) + " " + str(int(preds[i] >= 0.5)) + '\n')
        fw.close()

def output_logits(output_path, logits):
    '''
    Function:
            Output raw prediction values to a file
    
    Parameters:
            output_path (str): path to the file to output the preds to
    '''
    with open(output_path, 'w+') as fw:
        for value in logits:
            fw.write(str(value) + '\n')
        fw.close()

def examine_preds(y_pred, y_test):
    '''
        Function:
            Print analysis of model predictions with binary classification metrics

        Parameters:
            y_pred (list(int)): model predictions
            y_test (list(int)): true data labels
    '''
    print(f'Precision: {precision_score(y_test, y_pred)}')
    print(f'Recall: {recall_score(y_test, y_pred)}')
    print(f'F1: {f1_score(y_test, y_pred)}')
    print(f'ROC AUC: {roc_auc_score(y_test, y_pred)}')
    print(f'Confusion Matrix: \n{confusion_matrix(y_test, y_pred)}')


def create_lgbm_model(model_output_path,
                      preds_output_path,
                      logit_output_path,
                      x_train, y_train,
                      x_test, y_test,
                      dump_logits):
    '''
        Function:
            Create a lightGBM classifier

        Parameters:
            model_output_path (string): path to model output file
            preds_output_path (string): path to prediction output file
            logit_output_path (string): path to logit output file
            x_train (list(list(float))): train data
            y_train (list(int)): train data labels
            x_test (list(list(float))): test data
            y_test (list(int)): test data labels
            dump_logits (bool): whether or not to output logits (used in determinism test)
    '''
    from flaml.default import LGBMClassifier

    model = LGBMClassifier(random_state=RANDOM_SEED)
    model.fit(x_train, y_train)
    
    preds = model.predict(x_test)
    examine_preds(preds, y_test);
    
    output_preds(preds_output_path, preds, y_test)
    if (dump_logits):
        output_logits(logit_output_path, preds)
        
    with open(model_output_path, 'wb') as f:
        pickle.dump(model.booster_, f)
        f.close()

def create_neural_network(model_output_path,
                          preds_output_path,
                          logit_output_path,
                          X_train, y_train,
                          X_test, y_test,
                          dump_logits,
                          epochs):
    '''
        Function: 
            Create a neural network

        Parameters:
            model_output_path (string): path to model output file
            preds_output_path (string): path to prediction output file
            logit_output_path (string): path to logit output file
            x_train (list(list(float))): train data
            y_train (list(int)): train data labels
            x_test (list(list(float))): test data
            y_test (list(int)): test data labels
            dump_logits (bool): whether or not to output logits (used in determinism test)
            epochs (int): number of epochs to train for
    '''
    import numpy as np
    import tensorflow as tf
    from tensorflow import keras
    from tensorflow.keras import layers

    tf.random.set_seed(RANDOM_SEED);
    tf.keras.utils.set_random_seed(RANDOM_SEED)
    
    X_train = tf.constant(X_train)
    X_test = tf.constant(X_test)
    
    normalizer = layers.Normalization(axis=-1)
    normalizer.adapt(X_train)

    def binary_focal_loss(gamma=2.0, alpha=0.25):
        def loss(y_true, y_pred):
            eps = keras.backend.epsilon()
            y_pred = keras.backend.clip(y_pred, eps, 1 - eps)
            pt = tf.where(tf.equal(y_true, 1), y_pred, 1 - y_pred)
            w = tf.where(tf.equal(y_true, 1), alpha, 1 - alpha)
            return -tf.reduce_mean(w * tf.pow(1 - pt, gamma) * tf.math.log(pt))
        return loss

    USE_FOCAL = False
    def make_model(input_dim=18, l2=1e-4, dropout=0.25):
        inp = keras.Input(shape=(input_dim,), name="features")
        x = normalizer(inp)
        x = layers.Dense(64, activation=tf.nn.leaky_relu, kernel_regularizer=keras.regularizers.l2(l2))(x)
        x = layers.Dropout(dropout)(x)
        x = layers.Dense(32, activation=tf.nn.leaky_relu, kernel_regularizer=keras.regularizers.l2(l2))(x)
        x = layers.Dropout(dropout)(x)
        out = layers.Dense(1, activation="sigmoid", name="out")(x)
        model = keras.Model(inp, out)
        return model

    model = make_model(18)

    opt = keras.optimizers.Adam(learning_rate=3e-4)
    loss_fn = binary_focal_loss(gamma=2.0, alpha=0.25) if USE_FOCAL else keras.losses.BinaryCrossentropy()

    model.compile(
        optimizer=opt,
        loss=loss_fn,
        metrics=[
            keras.metrics.AUC(curve="ROC", name="auc"),
            keras.metrics.AUC(curve="PR", name="pr_auc"),
            keras.metrics.Precision(name="precision"),
            keras.metrics.Recall(name="recall"),
        ],
    )

    neg = (y_train == 0).sum()
    pos = (y_train == 1).sum()
    pos_weight = neg / max(pos, 1)
    class_weight = {0: 1.0, 1: pos_weight}

    history = model.fit(
        X_train, y_train,
        epochs=epochs,
        class_weight=class_weight
    )

    preds = np.array(model.predict(X_test, verbose=False).ravel())
            
    if (dump_logits):
        output_logits(logit_output_path, preds)

    preds = np.round(preds)
    output_preds(preds_output_path, preds, y_test)
    examine_preds(preds, y_test);

    model.save(model_output_path)

def tune_gbdt(model_output_path,
              preds_output_path,
              logit_output_path,
              x_train, y_train,
              x_test, y_test,
              dump_logits):
    '''
        Function:
            Hyperparameter tuning job for a GDBT model

        Parameters:
            model_output_path (string): path to model output file
            preds_output_path (string): path to prediction output file
            logit_output_path (string): path to logit output file
            x_train (list(list(float))): train data
            y_train (list(int)): train data labels
            x_test (list(list(float))): test data
            y_test (list(int)): test data labels
            dump_logits (bool): whether or not to output logits (used in determinism test)
    '''

    from sklearn.model_selection import GridSearchCV
    from sklearn.ensemble import GradientBoostingClassifier
    
    param_grid = {
        'n_estimators': [200, 400, 600, 800, 1000],
        'learning_rate': [0.1, 0.2, 0.3],
        'subsample': [.1, .25, 1.0],
        'max_depth': [5, 6]
    }

    grid_search = GridSearchCV(estimator=GradientBoostingClassifier(random_state=RANDOM_SEED),
                               param_grid=param_grid,
                               cv=5,
                               scoring='recall',
                               n_jobs=-1)

    grid_search.fit(x_train, y_train)
    best_model = grid_search.best_estimator_
    print(best_model.get_params())
    
    preds = best_model.predict(x_test)
    examine_preds(preds, y_test);
    
    output_preds(preds_output_path, preds, y_test)
    if (dump_logits):
        output_logits(logit_output_path, preds)
    
    with open(model_output_path, 'wb') as f:
        pickle.dump(best_model, f)
        f.close()

def create_gbdt_model(model_output_path,
                      preds_output_path,
                      logit_output_path,
                      x_train, y_train,
                      x_test, y_test,
                      dump_logits):
    '''
        Function:
            Create a GDBT model

        Parameters:
            model_output_path (string): path to model output file
            preds_output_path (string): path to prediction output file
            logit_output_path (string): path to logit output file
            x_train (list(list(float))): train data
            y_train (list(int)): train data labels
            x_test (list(list(float))): test data
            y_test (list(int)): test data labels
            dump_logits (bool): whether or not to output logits (used in determinism test)
    '''

    from sklearn.ensemble import GradientBoostingClassifier

    model = GradientBoostingClassifier(n_estimators=600,
                                       subsample=1.0,
                                       learning_rate=0.3,
                                       random_state=RANDOM_SEED)
    model.fit(x_train, y_train)
    
    preds = model.predict(x_test)
    examine_preds(preds, y_test)

    output_preds(preds_output_path, preds, y_test)
    if (dump_logits):
        output_logits(logit_output_path, preds)
    
    with open(model_output_path, 'wb') as f:
        pickle.dump(model, f)
        f.close()

run_dir = config['paths']['core']['run_root'] + config['run']['run_id'] + '/'
combo_dir = run_dir + config['paths']['core']['data_dir'] + config['paths']['core']['data_combo_files_dir']
train_combos_path = combo_dir + config['paths']['core']['data_combo_' + args.reagent_ion + '_train']
test_combos_path = combo_dir + config['paths']['core']['data_combo_' + args.reagent_ion + '_test']

ion = args.reagent_ion
x_train, y_train = getData(train_combos_path)
x_test, y_test = getData(test_combos_path)

model_output_dir = run_dir + config['paths']['core']['model_output_dir']
pred_output_dir = run_dir + config['paths']['core']['pred_output_dir']
logit_output_dir = pred_output_dir + config['paths']['core']['pred_logit_output_dir']

dump_logits = args.dump_logits != None
logit_id = args.dump_logits if dump_logits else ''

if (args.lgbm):
    create_lgbm_model(f"{model_output_dir}{args.reagent_ion}_lgbm.pkl",
                      f"{pred_output_dir}{args.reagent_ion}_lgbm_preds.txt",
                      f"{logit_output_dir}{args.reagent_ion}_lgbm_logits_{logit_id}.txt",
                      x_train, y_train,
                      x_test, y_test,
                      dump_logits)
 
if (args.neural_net):
    create_neural_network(f"{model_output_dir}{args.reagent_ion}_nn.keras",
                          f"{pred_output_dir}{args.reagent_ion}_nn_preds.txt",
                          f"{logit_output_dir}{args.reagent_ion}_nn_logits_{logit_id}.txt",
                          x_train, y_train,
                          x_test, y_test,
                          dump_logits,
                          int(config['core']['model_hyperparams']['nn'][args.reagent_ion + '_reagent']['epochs']))

if (args.gradient_booster):
    create_gbdt_model(f"{model_output_dir}{args.reagent_ion}_gbdt.pkl",
                      f"{pred_output_dir}{args.reagent_ion}_gbdt_preds.txt",
                      f"{logit_output_dir}{args.reagent_ion}_gbdt_logits_{logit_id}.txt",
                      x_train, y_train,
                      x_test, y_test,
                      dump_logits)
