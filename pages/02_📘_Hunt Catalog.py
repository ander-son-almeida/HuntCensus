# -*- coding: utf-8 -*-
"""
Created on Fri Aug 12 23:15:51 2022

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


st.set_page_config(page_title="Emily Hunt Catalog",layout='wide', page_icon='📘')

st.markdown(
    """
    <style>
    body {
        background-color: #fff;
    }
    </style>
    """,
    unsafe_allow_html=True
)

###############################################################################
#load clusters
clusters = pd.read_parquet('data/parquet/clusters.parquet')

#select only open clusters
mask_oc = clusters['kind'] == 'o'
clusters = clusters[mask_oc]

# Interface: Select clusters name
list_clusters = clusters['name']
cluster_name = st.sidebar.selectbox(
    "Select open cluster:",
    (list(list_clusters)))

#aply filter
clusters = clusters[clusters['name'] == cluster_name]

###############################################################################
#load members
members_ship = pd.read_parquet('data/parquet/members.parquet')

#aply filter name
members_ship = members_ship[members_ship['name'] == cluster_name]

# #aply fiter probability and mag
# members_ship = members_ship[(members_ship['probability'] > 0.5)&(members_ship['phot_g_mean_mag'] <19.)]


RA = clusters['ra']
DEC = clusters['dec']
age = clusters['log_age_84'].iloc[0]
dist = (clusters['distance_84']/1000).iloc[0]
Av = clusters['a_v_84'].iloc[0]


# bar with fundamental parameters
st.sidebar.subheader("Fundamental parameters:")
st.sidebar.subheader("$log(age) = {}$".format(np.around(age,decimals=3)))
st.sidebar.subheader("$Dist. = {}~(kpc)$".format(np.around(dist,decimals=3)))
st.sidebar.subheader("$Av. = {}~(mag)$".format(np.around(Av,decimals=3)))


#Graphics
###############################################################################
# CMD 
# read isochrones
mod_grid, age_grid, z_grid = load_mod_grid()
filters = ['Gmag','G_BPmag','G_RPmag']
refMag = 'Gmag' 

grid_iso = get_iso_from_grid(age,(10.**0)*0.0152,filters,refMag, nointerp=False)
fit_iso = make_obs_iso(filters, grid_iso, dist, Av, gaia_ext = True) 

cor_obs = members_ship['phot_bp_mean_mag']-members_ship['phot_rp_mean_mag']
absMag_obs = members_ship['phot_g_mean_mag']

cmd_scatter = pd.DataFrame({'G_BPmag - G_RPmag': cor_obs, 'Gmag': absMag_obs})

cmd_iso = pd.DataFrame({'G_BPmag - G_RPmag': fit_iso['G_BPmag']-fit_iso['G_RPmag'], 
                        'Gmag': fit_iso['Gmag']})

fig1 = px.scatter(cmd_scatter, x = 'G_BPmag - G_RPmag', y = 'Gmag',
                  opacity=0.3)

fig2 = px.line(cmd_iso, x = 'G_BPmag - G_RPmag', y = 'Gmag', color_discrete_sequence=['red'])

fig = go.Figure(data = fig1.data + fig2.data).update_layout(coloraxis=fig1.layout.coloraxis)
fig.update_layout(xaxis_title= 'G_BP - G_RP (mag)',
                  yaxis_title="G (mag)",
                  coloraxis_colorbar=dict(title="M☉"),
                  yaxis_range=[20,5])

###############################################################################	   
# RA x DEC 
ra_dec = pd.DataFrame({'RA': members_ship['ra'], 
                       'DEC': members_ship['dec']})

fig_ra_dec = px.scatter(ra_dec, x = 'RA', y = 'DEC')

###############################################################################	

container1 = st.container()
col1, col2 = st.columns(2)


with container1:
    
    
    with col1:
        st.subheader("CMD")
        st.plotly_chart(fig, use_container_width=True)

    with col2:
        st.subheader("Distribution RA and DEC")
        st.plotly_chart(fig_ra_dec, use_container_width=True)
        

    