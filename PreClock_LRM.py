#!/usr/bin/env python

import subprocess
import os
import sys
from ete3 import Tree
import dendropy
from sklearn.linear_model import LinearRegression
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse


def read_args(args):

    parser = argparse.ArgumentParser()
    parser.add_argument('raxml_bi_tree',
                        nargs = "?",
                        metavar = "raxml_bipartition_tree",
                        help = 'Input raxml bipartition tree file.',
                        type = str)
    parser.add_argument('mp_file',
                        nargs = "?",
                        help = "Input the mapping file delimited by tab. (1st col is tip label and 2nd col is time stamps. [ya]))",
                        metavar = 'mapping_file',
                        type = str)
    parser.add_argument('opt_png',
                        nargs = '?',
                        help = 'Specify the output name with suffix [.png].',
                        metavar = 'output',
                        type = str)

    return vars(parser.parse_args())



def tip_to_root(tre, tips_lst, mp):
    x_axis = []
    y_axis = []
    root = tre.get_tree_root()
    for tip in tips_lst:
        x_axis.append(float(mp[tip]) * -1)
        y_axis.append(float(tre.get_distance(tip, root)))
    return x_axis, y_axis

def lm_model(x, y):
    x = np.asarray(x)
    X = np.asarray(x)
    y = np.asarray(y)
    x = np.array(x.astype(np.float64)).reshape((-1, 1))
    y = np.array(y.astype(np.float64))
    model = LinearRegression().fit(x, y)
    r_sq = model.score(x, y)
    slope = model.coef_
    root_age = model.intercept_/slope
    residuals = np.square(y - model.predict(x))
    y_predict = model.predict(x)
    return r_sq, root_age[0], slope[0], residuals.sum(), y_predict, X, y


def rooting_trees(tre, tips_lst, mp):
    lm_dict = {}
    rooted_trees = {}
    tree = 0
    for n in tre.traverse():
        if n != tre:
            tree += 1
            tre.set_outgroup(n)
            tre.name = "tree_" + str(tree)
            lm_dict[tre.name] = tip_to_root(tre, tips_lst, mp)
            rooted_trees[tre.name] = tre.write()
    return lm_dict, rooted_trees

def lm_trees(lm_dict):
    trees, r_sq, root_age, slope, residuals_sum, y_predict, x_axis, y_axis =[], [], [], [], [], [], [], []
    for i in lm_dict:
        x = lm_dict[i][0]
        y = lm_dict[i][1]
        LM_model = lm_model(x, y)
        trees.append(i)
        r_sq.append(LM_model[0])
        root_age.append(LM_model[1])
        slope.append(LM_model[2])
        residuals_sum.append(LM_model[3])
        y_predict.append(LM_model[4])
        x_axis.append(LM_model[5])
        y_axis.append(LM_model[6])

    result_dict = {'trees': trees, 'r_sq': r_sq, 'root_age': root_age, 'slope': slope,
     'sq_residual_sum': residuals_sum, 'predicts': y_predict, 'x_axis': x_axis, 'y_axis': y_axis}

    return pd.DataFrame.from_dict(result_dict)  

def temp_est(ml_tree, mp_file, opt_png):

    tre = Tree(open(ml_tree).read().rstrip(), format = 2)
    mp = {i.strip().split('\t')[0]: i.strip().split('\t')[1] for i in open(mp_file).readlines()}
    tips_lst = [i.strip().split('\t')[0] for i in open(mp_file).readlines()]
    tree_lm = rooting_trees(tre, tips_lst, mp) #1. lm model for all trees, 2. tree in format2
    trees_docs = lm_trees(tree_lm[0])
    best_fit_tree = trees_docs[trees_docs.sq_residual_sum == trees_docs.sq_residual_sum.min()]
    x_axis = best_fit_tree.iloc[-1]['x_axis']
    y_predict = best_fit_tree.iloc[-1]['predicts']
    y_axis = best_fit_tree.iloc[-1]['y_axis']
    root_age = best_fit_tree.iloc[-1]['root_age']
    slope = best_fit_tree.iloc[-1]['slope']
    r_sq = best_fit_tree.iloc[-1]['r_sq']
    
    fig, ax = plt.subplots()
    ax.plot(x_axis, y_axis, 'ro')
    ax.plot(x_axis, y_predict, 'b')
    ax.set_xlabel('Years ago')
    ax.set_ylabel('Tip-to-Root distance')
    ax.text(0.05, 0.9, r'$r^2$: {}'.format(r_sq),
     horizontalalignment='left',
     verticalalignment='center',
     transform = ax.transAxes)
    ax.text(0.05, 0.85, r'Root age (x_axis intercept): {}'.format(root_age),
     horizontalalignment='left',
     verticalalignment='center',
     transform = ax.transAxes)
    ax.text(0.05, 0.80, r'slope: {}'.format(slope),
     horizontalalignment='left',
     verticalalignment='center',
     transform = ax.transAxes)  
    plt.savefig(opt_png)

def main():
    pars = read_args(sys.argv)
    temp_est(pars['raxml_bi_tree'], pars['mp_file'], pars['opt_png'])
    
 

if __name__ == '__main__':
    main()








