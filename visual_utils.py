#!/usr/bin/env python


import matplotlib.pyplot as plt
import numpy as numpyp

def hist_plot(ax, data, param_dict):

	out = ax.hist(data, **param_dict)

	return out
	