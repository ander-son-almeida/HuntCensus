# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 13:24:58 2023

@author: Anderson Almeida
"""

import numpy as np
import pandas as pd 
from astropy import units as u
from astropy.coordinates import SkyCoord
from scipy.optimize import curve_fit 
import streamlit as st
import plotly.express as px
import plotly.graph_objects as go
from oc_tools_padova_edr3 import *
import matplotlib.pyplot as plt

###############################################################################	
# clusters in common
cluster01 = np.genfromtxt('data/log-results-eDR3.txt', delimiter=';', names = True, 
                        dtype=None, encoding=None, autostrip=True)

cluster02 = pd.read_parquet('data/parquet/clusters.parquet')
cluster02 = cluster02.to_records()

ab, a_ind, b_ind = np.intersect1d(cluster01['name'],cluster02['name'], 
                                  return_indices=True)
cluster01 = cluster01[a_ind]
cluster02 = cluster02[b_ind]

#dist
dist_dt = pd.DataFrame({'dist': cluster01['dist']/1000, 'distance_84': cluster02['distance_84']/1000})
scatter_dist = px.scatter(dist_dt, x='dist', y='distance_84', opacity=0.3)

diag_line_dist = px.line(x = [(cluster01['dist']/1000).min(),(cluster02['distance_84']/1000).max()], 
                         y = [(cluster01['dist']/1000).min(),(cluster02['distance_84']/1000).max()], 
                         color_discrete_sequence=['red'])

fig_dist = go.Figure(data=scatter_dist.data + diag_line_dist.data)
fig_dist.update_layout(xaxis_title= 'Our',
                  yaxis_title="Hunt")

#age
age_dt = pd.DataFrame({'age': cluster01['age'], 'log_age_84': cluster02['log_age_84']})
scatter_age = px.scatter(age_dt, x='age', y='log_age_84', opacity=0.3)
diag_line_age = px.line(x = [cluster01['age'].min(),cluster02['log_age_84'].max()], 
                         y = [cluster01['age'].min(),cluster02['log_age_84'].max()], 
                         color_discrete_sequence=['red'])
fig_age = go.Figure(data=scatter_age.data + diag_line_age.data)
fig_age.update_layout(xaxis_title= 'Our',
                  yaxis_title="Hunt")


#av
age_dt = pd.DataFrame({'Av': cluster01['Av'], 'a_v_84': cluster02['a_v_84']})
scatter_av = px.scatter(age_dt, x='Av', y='a_v_84', opacity=0.3)
fig_av = go.Figure(data=scatter_av)
diag_line_av = px.line(x = [cluster01['Av'].min(),cluster02['a_v_84'].max()], 
                         y = [cluster01['Av'].min(),cluster02['a_v_84'].max()], 
                         color_discrete_sequence=['red'])
fig_av = go.Figure(data=fig_av.data + diag_line_av.data)
fig_av.update_layout(xaxis_title= 'Our',
                  yaxis_title="Hunt")

container1 = st.container()
col4, col5 = st.columns(2)

with container1:
    with col4:
        st.subheader('''
        The plots below compare the fundamental parameters of all clusters:
        ''')
        
    with col5:
        st.write('''
        
        ''')

container3 = st.container()
col6, col7, col8  = st.columns(3)

with container3:
    with col6:
        st.subheader("Distance")
        st.plotly_chart(fig_dist, use_container_width=True)

    with col7:
        st.subheader("age")
        st.plotly_chart(fig_age, use_container_width=True)
        
    with col8:
        st.subheader("Av")
        st.plotly_chart(fig_av, use_container_width=True)



