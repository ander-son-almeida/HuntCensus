# -*- coding: utf-8 -*-
"""
Created on Tue Aug 16 23:50:07 2022

@author: Anderson Almeida
"""




import streamlit as st

st.set_page_config(page_title="Home",layout='centered', page_icon='🔵')

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


st.sidebar.image("images/logo.png", use_column_width=True)

st.write('''
         
         
test test


    ''')
        
st.subheader('The page is still being updated!')

st.write('''

For any questions, information or collaborations, please contact us by email:
    andersonalmeida_sa@outlook.com or hmonteiro@unifei.edu.br
    
    ''')