# -*- coding: utf-8 -*-
"""
Created on Sat Mar 25 22:53:57 2023

@author: Anderson Almeida
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

dir = r'S:\Área de Trabalho\EmiliyHunt'
dir_our = r'S:\Área de Trabalho\Artigo Massa DR3'
version = '2023' 


cluster_name = 'Melotte_22'

###############################################################################
#load clusters
clusters = pd.read_parquet(dir + r'\data\parquet\clusters.parquet')

#select only open clusters
mask_oc = clusters['kind'] == 'o'
clusters = clusters[mask_oc]

#aply filter
clusters = clusters[clusters['name'] == cluster_name]
###############################################################################
#load members
members = pd.read_parquet(dir + r'\data\parquet\members.parquet')

#aply filter name
members = members[members['name'] == cluster_name]

#aply fiter probability and mag
members = members[(members['probability'] > 0.5)&(members['phot_g_mean_mag'] <19.)]

###############################################################################
# Emily fundamental parameters
age = clusters['log_age_84']
dist = clusters['distance_84']/1000 #kpc
Av = clusters['a_v_84']





###############################################################################
plt.figure()
plt.scatter(members['ra'], members['dec'])
plt.xlabel('RA')
plt.ylabel('DEC')

plt.figure()
plt.scatter(members['phot_bp_mean_mag'] - members['phot_rp_mean_mag'], members['phot_g_mean_mag'])
plt.gca().invert_yaxis()
plt.ylabel(r'$G_{mag}$')
plt.xlabel(r'$G_{BP}-G_{RP}$')


###############################################################################
#lendo membros
data_dir = dir_our + r'\results_eDR3_likelihood_{}\membership_data_edr3\\'.format(version)
file = '{}_data_stars.dat'.format(cluster_name)
data_obs = np.genfromtxt(data_dir + file, delimiter=';', names=True, dtype=None)

#corte de magnitude
data_obs = data_obs[(data_obs['Pmemb'] > 0.5)&(data_obs['Gmag'] <19.)]


###############################################################################
# lendo do log results os parametros fundamentais
cluster = np.genfromtxt(dir_our + r'\results_eDR3_likelihood_{}\results\log-results-eDR3.txt'.format(version), 
                                        delimiter=';', names = True, dtype=None, encoding=None, autostrip=True)

#aply filter name
cluster = cluster[cluster['name'] == cluster_name]

###############################################################################
#our fundamental parameters

age_our = cluster['age'][0]
dist_our = cluster['dist'][0]/1000 # em kpc
# FeH_our = cluster['FeH'][0]
Av_our = cluster['Av'][0]



x = ['log(age)', 'Dist. (kpc)', 'Av.(mag)']
x1 = [1,2,3]
y1 = [age_our, dist_our, Av_our]
y2 = [age, dist, Av]

# criar subplots individuais
fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(10, 4))

# definir a largura das barras
bar_width = 0.35

for i in range(3):

    pos = np.arange(len(x1)) + bar_width
    
    axs[i].bar(x1[i], y1[i], bar_width, alpha=0.5)
    axs[i].bar(x1[i]+ bar_width, y2[i], bar_width, alpha=0.5)
    axs[i].set_xlabel(x[i])
    axs[i].set_xticks([])
    
fig.legend(['Our', 'Emily'])





































