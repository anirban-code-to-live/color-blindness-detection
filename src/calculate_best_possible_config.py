
# coding: utf-8

# In[39]:

import numpy as np
import pickle
import operator as op
import math


# In[40]:

def ncr(n, r):
    r = min(r, n-r)
    numer = reduce(op.mul, xrange(n, n-r, -1), 1)
    denom = reduce(op.mul, xrange(1, r+1), 1)
    return numer//denom


# In[41]:

red_exon_count = [97.0, 237.0, 107.5, 167.5, 303.5, 235.0]
green_exon_count = [97.0, 156.5, 99.5, 87.5, 217.5, 235.0]
red_exon_count = [int(count) for count in red_exon_count]
green_exon_count = [int(count) for count in green_exon_count]
total_count = []
assert len(red_exon_count) == len(green_exon_count)
for i in range(len(red_exon_count)):
    count = red_exon_count[i] + green_exon_count[i]
    total_count.append(count)
# print(total_count)


# In[42]:

config_1_red = [0.33, 0.33, 0.33, 0.33]
config_1_green = [0.66, 0.66, 0.66, 0.66]

config_2_red = [0.5, 0.5, 0, 0]
config_2_green = [0.5, 0.5, 1.0, 1.0]

config_3_red = [0.25, 0.25, 0.5, 0.5]
config_3_green = [0.75, 0.75, 0.5, 0.5]

config_4_red = [0.25, 0.25, 0.25, 0.5]
config_4_green = [0.75, 0.75, 0.75, 0.5]


# In[51]:

probabilities = []
prob_1 = 0
for i in range(len(config_1_red)):
    prob_1 += red_exon_count[i+1]*math.log(config_1_red[i]) + green_exon_count[i+1]*math.log(config_1_green[i]) + math.log(ncr(total_count[i+1], red_exon_count[i+1]))

prob_3 = 0
for i in range(len(config_3_red)):
    prob_3 += red_exon_count[i+1]*math.log(config_3_red[i]) + green_exon_count[i+1]*math.log(config_3_green[i]) + math.log(ncr(total_count[i+1], red_exon_count[i+1]))

prob_4 = 0
for i in range(len(config_4_red)):
    prob_4 += red_exon_count[i+1]*math.log(config_4_red[i]) + green_exon_count[i+1]*math.log(config_4_green[i]) + math.log(ncr(total_count[i+1], red_exon_count[i+1]))

print('Log probability for config-1 :: ' + str(prob_1))
print('Log probability for config-3 :: ' + str(prob_3))
print('Log probability for config-4 :: ' + str(prob_4))

print('')
print('Probability for config-1 :: ' + str(math.exp(prob_1)))
print('Probability for config-3 :: ' + str(math.exp(prob_3)))
print('Probability for config-4 :: ' + str(math.exp(prob_4)))

