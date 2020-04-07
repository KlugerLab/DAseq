import os, sys
import numpy as np
import pandas as pd

from keras import backend as K
from keras.models import Model
from keras.layers import Input, Dense
from keras.activations import relu
from keras.callbacks import EarlyStopping

def k_fold_split(x, p, k):
    return [x[p_] for p_ in np.split(p, k)]

def rev_split(xs, p):
    n = len(p)
    inv_p = np.empty(n)
    inv_p[p] = np.arange(n)
    inv_p = inv_p.astype(int)
    return np.concatenate(xs)[inv_p]

def make_splits(n, k):
    sizes = int(n / k)
    splits = [sizes] * (k - 1) + [sizes + n % k]
    for i in range(k - 1):
        splits[i+1] = splits[i] + splits[i+1]
    return splits[:-1]

def k_fold_predict(data, labels, k_folds, architecture=[8]*8, activations='relu', end_activation='sigmoid'):
    # os.environ['CUDA_VISIBLE_DEVICES'] = '4'
    # build layers

    layers = []

    for width in architecture:

        layers.append(Dense(width, activation=activations))

    layers.append(Dense(1, activation=end_activation))



    # build neural network

    input_shape = data.shape[1:]

    x = x0 = Input(shape=input_shape)

    for layer in layers:

        x = layer(x)



    model = Model(inputs=x0, outputs=x)



    y_tests = []



    p = np.random.permutation(len(data))



    for i in range(k_folds):

        val_idx = (i - 1) % k_folds

        test_idx = i

        k_folds_ = make_splits(len(data), k_folds)

        x_full, y_full = k_fold_split(data, p, k_folds_), k_fold_split(labels, p, k_folds_)

        x_val_, y_val_ = x_full[val_idx], y_full[val_idx]

        x_test_, y_test_ = x_full[test_idx], y_full[test_idx]
	# remove VALIDATION AND TEST SETS from the training set
        if k_folds > 2:
            del x_full[val_idx]
            del y_full[val_idx]
        # if val_idx came before test_idx, we have to remove the (test_idx - 1)th element (as we have already deleted val_idx so the index corresponding to the test set has changed)
        if k_folds > 1 and val_idx < test_idx:
            del x_full[(test_idx-1)]
            del y_full[(test_idx-1)]
        elif k_folds > 1:
        # otherwise, simply remove the (test_idx)th element
            del x_full[test_idx]
            del y_full[test_idx]

        x_train_, y_train_ = np.concatenate(x_full), np.concatenate(y_full)

        model.compile('adam', loss='binary_crossentropy', metrics=['acc'])



        epochs = 1000

        batch_size = len(x_train_)



        model.fit(

            x=x_train_,

            y=y_train_,

            epochs=epochs,

            batch_size=batch_size,

            validation_data=[x_val_, y_val_],

            callbacks=[EarlyStopping(patience=10)],

            verbose=0)



#        print("Finished {} / {} folds.".format(i + 1, k_folds))



        y_tests.append(model.predict(x_test_).reshape((-1,)))



    y_full = rev_split(y_tests, p)



    return y_full



def k_fold_predict_linear(data, labels, k_folds):
    return k_fold_predict(data, labels, k_folds, architecture=[], activations=None, end_activation='sigmoid')


