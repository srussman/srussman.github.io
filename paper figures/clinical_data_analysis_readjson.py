#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  7 12:40:59 2021

@author: srussman@ucsd.edu

This code reads .json data from the clinical IONM equipment
"""

import json
import os

os.chdir('/Users/samantha/Downloads')
data = open('yourfile.json')

jsondata = data.read()

obj = json.loads(jsondata)



cases=obj['Cases']

# Run through all modes and all names
m = 0
n=2

trials = cases[0]['Modes'][m]['Trials']
modes = cases[0]['Modes']
mode_number = len(cases[0]['Modes']) #choice of L UE/LE EMG etc in modes list
mode = cases[0]['Modes'][m]['Name']
name_number = len(cases[0]['Modes'][m]['Trials'][0]['Traces'])
name = cases[0]['Modes'][m]['Trials'][0]['Traces'][n]['Channel']['Name']

trace_data = []
timestamp_list = []
events = []
trace_scalar = [];
for trial in trials:
    trace = trial['Traces'][n]
    trace_scalar.append(trace['TraceDataScalar'])
    timestamp_list.append(trial['Timestamp'])
    trace_data.append(trace['TraceData'])
    # Find if labels are on a given trial

# dict_events = cases[0]['Events']
# messages = []
# timestamp_labels = []
# for event in dict_events:
#     message = event['Message']
#     timestamp = event["Timestamp"]
#     messages.append(message)
#     timestamp_labels.append(timestamp)

#
import numpy as np
import matplotlib.pyplot as plt

trace_data_np = np.array(trace_data)
#fig1, ax1 = plt.subplots(1,1)
#ax1.plot(timestamp,trace_data_np[:,0])

fig2, ax2 = plt.subplots(1,1)
for i in range(0,20): #len(trace_data)):
 #mean = np.mean(trace_data_np,0)

    plt.plot(-trace_data_np[i,:])

#plt.plot(-mean)
#%%

# Save array as csv file
filename = "0319_%s.csv" % (mode + name + "_trace_data")
os.chdir('/Users/samantha/Desktop/Clinical_data')
np.savetxt(filename, trace_data_np, delimiter=",")

# filename = "%s.csv" % (mode + name + "_messages")
# os.chdir('/Users/samantha/Desktop/Clinical_data')
# np.savetxt(filename, messages, delimiter=",", fmt='%s')

# filename = "0319_%s.csv" % (mode + name + "_timestamp_labels")
# os.chdir('/Users/samantha/Desktop/Clinical_data')
# np.savetxt(filename, timestamp_labels, delimiter=",", fmt='%s')

filename = "0319_%s.csv" % (mode + name + "_trace_data_scalar")
os.chdir('/Users/samantha/Desktop/Clinical_data')
np.savetxt(filename, trace_scalar, delimiter=",")

filename = "0319_%s.csv" % (mode + name + "_timestamp_list")
os.chdir('/Users/samantha/Desktop/Clinical_data')
np.savetxt(filename, timestamp_list, delimiter=",")