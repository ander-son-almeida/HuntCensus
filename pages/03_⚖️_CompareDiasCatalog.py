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
import requests
from io import BytesIO

st.set_page_config(page_title="Compare Dias Catalog",layout='wide', page_icon='⚖️')

###############################################################################
# CATALOG DIAS
###############################################################################
# read isochrones
mod_grid, age_grid, z_grid = load_mod_grid()
filters = ['Gmag','G_BPmag','G_RPmag']
refMag = 'Gmag' 

# fundamental parameters
cluster_our = np.genfromtxt('data/log-results-eDR3.txt', delimiter=';', names = True, 
                        dtype=None, encoding=None, autostrip=True)

# Interface: Select cluster_ours name
list_cluster_ours = cluster_our['name']
cluster_our_name = st.sidebar.selectbox(
    "Select open cluster_our:",
    (list(list_cluster_ours)))
  
# read memberships
url = 'https://github.com/ander-son-almeida/EmiliyHuntCensus/blob/main/data/membership_data_edr3/{}_data_stars.dat'.format(cluster_our_name)
response = requests.get(url)

# file = 'data/membership_data_edr3/{}_data_stars.npy'.format(cluster_our_name)
members_ship = np.load(BytesIO(response.content), allow_pickle=False)

# select fundamental parameters cluster_our	
ind = np.where(cluster_our['cluster_our'] == cluster_our_name)

RA_our = cluster_our['RA_ICRS'][ind]
DEC_our = cluster_our['DE_ICRS'][ind]
age_our = cluster_our['age'][ind]
e_age_our = cluster_our['e_age'][ind]
dist_our = cluster_our['dist'][ind]
e_dist_our = cluster_our['e_dist'][ind]
FeH_our = cluster_our['FeH'][ind]
e_FeH_our = cluster_our['e_FeH'][ind]
Av_our = cluster_our['Av'][ind]
e_Av_our = cluster_our['e_Av'][ind]


# bar with fundamental parameters
# st.sidebar.subheader("Fundamental parameters:")
# st.sidebar.subheader("$log(age) = {} \pm {}$".format(age[0], e_age[0]))
# st.sidebar.subheader("$Dist. = {} \pm {}~(kpc)$".format(dist[0],e_dist[0]))
# st.sidebar.subheader("$Av. = {} \pm {}~(mag)$".format(Av[0],e_Av[0]))
# st.sidebar.subheader("$FeH = {} \pm {}$".format(FeH[0],e_FeH[0]))


###############################################################################
# CATALOG EMILY HUNT
###############################################################################

#load clusters
clustersEmily = pd.read_parquet('data/parquet/clusters.parquet')

#select only open clusters
mask_oc = clustersEmily['kind'] == 'o'
clustersEmily = clustersEmily[mask_oc]

#aply filter
clustersEmily = clustersEmily[clustersEmily['name'] == cluster_our_name]


###############################################################################
#load members
members_ship_Emily = pd.read_parquet('data/parquet/members.parquet')

#aply filter name
members_ship_Emily = members_ship_Emily[members_ship_Emily['name'] == cluster_our_name]


RA = clustersEmily['ra']
DEC = clustersEmily['dec']
age = clustersEmily['log_age_84'].iloc[0]
dist = (clustersEmily['distance_84']/1000).iloc[0]
Av = clustersEmily['a_v_84'].iloc[0]


# bar with fundamental parameters
# st.sidebar.subheader("Fundamental parameters:")
# st.sidebar.subheader("$log(age) = {}$".format(np.around(age,decimals=3)))
# st.sidebar.subheader("$Dist. = {}~(kpc)$".format(np.around(dist,decimals=3)))
# st.sidebar.subheader("$Av. = {}~(mag)$".format(np.around(Av,decimals=3)))




#Graphics
###############################################################################
# CMD Emily

cor_obs_emily = members_ship_Emily['phot_bp_mean_mag']-members_ship_Emily['phot_rp_mean_mag']
absMag_obs_emily = members_ship_Emily['phot_g_mean_mag']

cmd_scatter_emily = pd.DataFrame({'G_BPmag - G_RPmag': cor_obs_emily, 'Gmag': absMag_obs_emily})

cmd_emily = px.scatter(cmd_scatter_emily, x = 'G_BPmag - G_RPmag', y = 'Gmag',
                  opacity=0.6)

fig_CMD_emily = go.Figure(data = cmd_emily.data).update_layout(coloraxis=cmd_emily.layout.coloraxis)
fig_CMD_emily.update_layout(xaxis_title= 'G_BP - G_RP (mag)',
                  yaxis_title="G (mag)",
                  coloraxis_colorbar=dict(title="M☉"),
                  yaxis_range=[20,5])



#CMD Dias

grid_iso = get_iso_from_grid(age,(10.**FeH)*0.0152,filters,refMag, nointerp=False)
fit_iso = make_obs_iso(filters, grid_iso, dist, Av, gaia_ext = True) 
cor_obs = members_ship['BPmag']-members_ship['RPmag']
absMag_obs = members_ship['Gmag']

cmd_scatter_dias = pd.DataFrame({'G_BPmag - G_RPmag': cor_obs, 'Gmag': absMag_obs})

cmd_iso = pd.DataFrame({'G_BPmag - G_RPmag': fit_iso['G_BPmag']-fit_iso['G_RPmag'], 
                        'Gmag': fit_iso['Gmag']})

fig_CMD_dias1 = px.scatter(cmd_scatter_dias, x = 'G_BPmag - G_RPmag', y = 'Gmag',
                  opacity=0.6)

fig_CMD_dias2 = px.line(cmd_iso, x = 'G_BPmag - G_RPmag', y = 'Gmag')

fig_CMD_dias = go.Figure(data = fig_CMD_dias1.data + fig_CMD_dias2.data).update_layout(coloraxis=fig_CMD_dias1.layout.coloraxis)
fig_CMD_dias.update_layout(xaxis_title= 'G_BP - G_RP (mag)',
                  yaxis_title="G (mag)",
                  coloraxis_colorbar=dict(title="M☉"),
                  yaxis_range=[20,5])

# ###############################################################################	   
# # RA x DEC 
# # the mass is arranged according to the mass of the primary
# ind = np.argsort(members_ship['mass'])

# ra_dec = pd.DataFrame({'RA': members_ship['RA_ICRS'][ind], 
#                         'DEC': members_ship['DE_ICRS'][ind], 'Mass': members_ship['mass'][ind]})

# fig_ra_dec = px.scatter(ra_dec, x = 'RA', y = 'DEC', color= 'Mass', 
#                         color_continuous_scale = 'jet_r')

# ###############################################################################	
# # Segregation Mass
# Mc = 1.0
# c1 = SkyCoord(ra=RA*u.degree, dec=DEC*u.degree, distance=dist*u.kpc) # open cluster_our center
# c2 = SkyCoord(ra=members_ship['RA_ICRS']*u.degree, dec=members_ship['DE_ICRS']*u.degree, distance=dist*u.kpc) 
# mass_members_ship = members_ship['mass'] + members_ship['comp_mass']
# star_dist = np.array(c1.separation_3d(c2)*1000)


# seg1 = pd.DataFrame({'Mc < 1M☉': star_dist[mass_members_ship < Mc]})
# seg2 = pd.DataFrame({'Mc > 1M☉': star_dist[mass_members_ship > Mc]})

# seg = pd.concat([seg1,seg2], axis=1)
# seg = seg.fillna(0)

# hist, bin_edges = np.histogram(star_dist[mass_members_ship < Mc], density=True)
# hist2, bin_edges2 = np.histogram(star_dist[mass_members_ship > Mc], density=True)

# xaxis_max = np.concatenate((bin_edges, bin_edges2), axis=0)
# yaxis_max = np.concatenate((hist, hist2), axis=0)

# seg = px.histogram(seg, histnorm='probability density', opacity=0.7)
# seg.add_vline(x=np.average(star_dist[mass_members_ship < Mc]), line_dash = 'dash', line_color = 'blue')
# seg.add_vline(x=np.average(star_dist[mass_members_ship > Mc]), line_dash = 'dash', line_color = 'red')

# seg.update_layout(xaxis_title= 'Distance (pc)',
#                   legend={'title_text':''},
#                   yaxis_title='Count',
#                   xaxis_range=[1,xaxis_max.max()],
#                   yaxis_range=[0,yaxis_max.max()+0.02])


# ###############################################################################	
# # FM Single
# ind_indv = members_ship['comp_mass'] == 0
# mass =  members_ship['mass'][ind_indv]

# (alpha_high_mass_sing, alpha_low_mass_sing, Mc_sing, offset_sing, alpha_high_mass_error_sing, \
#     alpha_low_mass_error_sing, Mc_error_sing, offset_error_sing, mass_cnt, mass_cnt_er, 
#     mass_bin_ctr, inv_mass, inv_mass_wd, popt) = fit_MF(mass, 'Sigle')

# xplot = np.linspace(mass_bin_ctr.min(),mass_bin_ctr.max(),1000)

# fm_ind = pd.DataFrame({'mass_bin_ctr': mass_bin_ctr, 'mass_cnt': mass_cnt, 'mass_cnt_er': mass_cnt_er})
# fm_ind_adj = pd.DataFrame({'xplot': xplot, 'ajuste': twosided_IMF(xplot, *popt)})

# plot_ind1 = px.scatter(fm_ind, x="mass_bin_ctr", y="mass_cnt", error_y="mass_cnt_er")
# plot_ind2 = px.line(fm_ind_adj, x = 'xplot', y = 'ajuste', color_discrete_sequence = ['orange'])
# plot_ind = go.Figure(data = plot_ind1.data + plot_ind2.data)

# plot_ind.update_layout(xaxis_title = 'log(M☉)', yaxis_title='ξ(log(M☉)')

# ###############################################################################	
# # FM Primary

# ind_bin = members_ship['comp_mass'] > 0
# mass =  members_ship['mass'][ind_bin]

# (alpha_high_mass_prim, alpha_low_mass_prim, Mc_prim, offset_prim, alpha_high_mass_error_prim, \
#     alpha_low_mass_error_prim, Mc_error_prim, offset_error_prim, mass_cnt, mass_cnt_er, 
#     mass_bin_ctr, inv_mass, inv_mass_wd, popt) = fit_MF(mass, 'Primary')

# xplot = np.linspace(mass_bin_ctr.min(),mass_bin_ctr.max(),1000)

# fm_prim = pd.DataFrame({'mass_bin_ctr': mass_bin_ctr, 'mass_cnt': mass_cnt, 'mass_cnt_er': mass_cnt_er})
# fm_prim_adj = pd.DataFrame({'xplot': xplot, 'ajuste': twosided_IMF(xplot, *popt)})

# plot_prim1 = px.scatter(fm_prim, x="mass_bin_ctr", y="mass_cnt", error_y="mass_cnt_er")
# plot_prim2 = px.line(fm_prim_adj, x = 'xplot', y = 'ajuste', color_discrete_sequence = ['orange'])
# plot_prim = go.Figure(data = plot_prim1.data + plot_prim2.data)

# plot_prim.update_layout(xaxis_title = 'log(M☉)', yaxis_title='ξ(log(M☉)')

# ###############################################################################	
# # FM Secundary

# ind_bin = members_ship['comp_mass'] > 0
# mass =  members_ship['comp_mass'][ind_bin]

# title = 'Secundárias \n 'r'$\alpha_A = {} \pm {}$; $\alpha_B = {} \pm {}$; $M_c = {} \pm {}$'

# (alpha_high_mass_sec, alpha_low_mass_sec, Mc_sec, offset_sec, alpha_high_mass_error_sec, \
#     alpha_low_mass_error_sec, Mc_error_sec, offset_error_sec, mass_cnt, mass_cnt_er, 
#     mass_bin_ctr, inv_mass, inv_mass_wd, popt) = fit_MF(mass, 'Secundary')

# xplot = np.linspace(mass_bin_ctr.min(),mass_bin_ctr.max(),1000)

# fm_sec = pd.DataFrame({'mass_bin_ctr': mass_bin_ctr, 'mass_cnt': mass_cnt, 'mass_cnt_er': mass_cnt_er})
# fm_sec_adj = pd.DataFrame({'xplot': xplot, 'ajuste': twosided_IMF(xplot, *popt)})

# plot_sec1 = px.scatter(fm_sec, x="mass_bin_ctr", y="mass_cnt", error_y="mass_cnt_er")
# plot_sec2 = px.line(fm_sec_adj, x = 'xplot', y = 'ajuste', color_discrete_sequence = ['orange'])
# plot_sec = go.Figure(data = plot_sec1.data + plot_sec2.data)

# plot_sec.update_layout(xaxis_title = 'log(M☉)', yaxis_title='ξ(log(M☉)')

# ###############################################################################	
# # FM Binary

# ind_bin = members_ship['comp_mass'] > 0
# mass =  np.concatenate((members_ship['mass'][ind_bin], members_ship['comp_mass'][ind_bin]), axis = 0)

# (alpha_high_mass_bin, alpha_low_mass_bin, Mc_bin, offset_bin, alpha_high_mass_error_bin, \
#     alpha_low_mass_error_bin, Mc_error_bin, offset_error_bin, mass_cnt, mass_cnt_er, 
#     mass_bin_ctr, inv_mass, inv_mass_wd, popt) = fit_MF(mass, 'Binary')

# xplot = np.linspace(mass_bin_ctr.min(),mass_bin_ctr.max(),1000)

# fm_bin = pd.DataFrame({'mass_bin_ctr': mass_bin_ctr, 'mass_cnt': mass_cnt, 'mass_cnt_er': mass_cnt_er})
# fm_bin_adj = pd.DataFrame({'xplot': xplot, 'ajuste': twosided_IMF(xplot, *popt)})

# plot_bin1 = px.scatter(fm_bin, x="mass_bin_ctr", y="mass_cnt", error_y="mass_cnt_er")
# plot_bin2 = px.line(fm_bin_adj, x = 'xplot', y = 'ajuste', color_discrete_sequence = ['orange'])
# plot_bin = go.Figure(data = plot_bin1.data + plot_bin2.data)

# plot_bin.update_layout(xaxis_title = 'log(M☉)', yaxis_title='ξ(log(M☉)')

container1 = st.container()
col1, col2  = st.columns(2)


with container1:
    with col1:
        st.subheader("CMD Emily")
        st.plotly_chart(fig_CMD_emily, use_container_width=True)

    with col2:
        st.subheader("CMD Dias")
        st.plotly_chart(fig_CMD_dias, use_container_width=True)
        
    # with col3:
    #     st.subheader("Segregation ratio")
    #     st.plotly_chart(seg, use_container_width=True)


# container2 = st.container()
# col4, col5 = st.columns(2)

# with container2:
    
#     st.header("Mass functions")
#     with col4:
#         st.subheader("Single")
#         st.info('$\\alpha_{{A}}={}~\pm~{};~'
#                 '\\alpha_{{B}}={}~\pm~{};~'
#                 'M_{{C}}={}~\pm~{}$'.format(np.around(alpha_high_mass_sing,decimals=2), 
#                                                 np.around(alpha_high_mass_error_sing,decimals=2),
#                                                 np.around(alpha_low_mass_sing,decimals=2),
#                                                 np.around(alpha_low_mass_error_sing,decimals=2),
#                                                 np.around(Mc_sing,decimals=2),
#                                                 np.around(Mc_error_sing,decimals=2)
#                                                 ))
#         st.plotly_chart(plot_ind, use_container_width=True)
    
#     with col5:
#         st.subheader("Primary")
#         st.info('$\\alpha_{{A}}={}~\pm~{};~'
#                 '\\alpha_{{B}}={}~\pm~{};~'
#                 'M_{{C}}={}~\pm~{}$'.format(np.around(alpha_high_mass_prim,decimals=2), 
#                                                 np.around(alpha_high_mass_error_prim,decimals=2),
#                                                 np.around(alpha_low_mass_prim,decimals=2),
#                                                 np.around(alpha_low_mass_error_prim,decimals=2),
#                                                 np.around(Mc_prim,decimals=2),
#                                                 np.around(Mc_error_prim,decimals=2)
#                                                 ))
#         st.plotly_chart(plot_prim, use_container_width=True)

# container3 = st.container()
# col6, col7 = st.columns(2)

# with container3:
#     with col6:
#         st.subheader("Secundary")
#         st.info('$\\alpha_{{A}}={}~\pm~{};~'
#                 '\\alpha_{{B}}={}~\pm~{};~'
#                 'M_{{C}}={}~\pm~{}$'.format(np.around(alpha_high_mass_sec,decimals=2), 
#                                                 np.around(alpha_high_mass_error_sec,decimals=2),
#                                                 np.around(alpha_low_mass_sec,decimals=2),
#                                                 np.around(alpha_low_mass_error_sec,decimals=2),
#                                                 np.around(Mc_sec,decimals=2),
#                                                 np.around(Mc_error_sec,decimals=2)
#                                                 ))
#         st.plotly_chart(plot_sec, use_container_width=True)

#     with col7:
#         st.subheader("Binary")
#         st.info('$\\alpha_{{A}}={}~\pm~{};~'
#                 '\\alpha_{{B}}={}~\pm~{};~'
#                 'M_{{C}}={}~\pm~{}$'.format(np.around(alpha_high_mass_bin,decimals=2), 
#                                                 np.around(alpha_high_mass_error_bin,decimals=2),
#                                                 np.around(alpha_low_mass_bin,decimals=2),
#                                                 np.around(alpha_low_mass_error_bin,decimals=2),
#                                                 np.around(Mc_bin,decimals=2),
#                                                 np.around(Mc_error_bin,decimals=2)
#                                                 ))
#         st.plotly_chart(plot_bin, use_container_width=True)



































