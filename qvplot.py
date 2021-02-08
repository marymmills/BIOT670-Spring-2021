'''
    BIOT670 Capstone Project - Quad Viewer
    plotting functions
'''
import json
import math

import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_table
import pandas as pd
import numpy as np
import re
import plotly.graph_objects as go
from dash.dependencies import Input, Output

app = dash.Dash(__name__)
styles = {
    'pre': {
        'border': 'thin lightgrey solid',
        'overflowX': 'scroll'
    }
}
xpos = 'DB_CORO_MEAN'
ypos = 'HET_CORO_MEAN'
xneg = 'KOHET_CORO_MEAN'
yneg = 'KOKO_CORO_MEAN'
dfdraft = pd.read_csv('coro_data.tsv', sep='\t')
df = dfdraft.replace(np.nan, 'None', regex=True)

df[xpos] = df[xpos].apply(lambda x: x)
df[ypos] = df[ypos].apply(lambda x: x)

fig = go.Figure()

color_dict={'Integrin':'limegreen',
            'Blood_Coagulation':'firebrick',
            'Cytoskeleton':'orangered',
            'Chemokine_Cytokine_Signaling':'tomato',
            'Chemokine_Cytokine_Signaling, Cytoskeleton':'royalblue',
            'Chemokine_Cytokine_Signaling, Cytoskeleton, Integrin':'seagreen',
            'Chemokine_Cytokine_Signaling, Cytoskeleton, Huntington':'wheat',
            'Chemokine_Cytokine_Signaling, Huntington, Integrin':'yellowgreen',
            'Chemokine_Cytokine_Signaling, Cytoskeleton, Huntington, Integrin':'violet',
            'Huntington, Integrin':'crimson',
            'Parkinson':'lightseagreen',
            'Cytoskeleton, Huntington':'aqua',
            'Chemokine_Cytokine_Signaling, Integrin':'palegreen',
            'Glycolysis, Huntington':'chocolate',
            'ATP_Synthesis':'brown',
            'ATP_Synthesis, Huntington':'burlywood',
            'Glycolysis':'mediumvioletred',
            'Glycolysis, Pyruvate_Metabolism':'cadetblue',
            'Huntington':'goldenrod',
            'Pyruvate_Metabolism':'darkkhaki',
            'Pyruvate_Metabolism, TCA_Cycle':'sienna',
            'TCA_Cycle':'darkred',
            'De_Novo_Purine_Biosynthesis':'mediumpurple',
            'None':'slategray'}

#qvplot.py
fig.add_scatter(x=df[xpos],
                y=df[ypos],
                mode='markers',
                marker=dict(size=10, color=[color_dict[k] for k in df['All_Pathways'].values]),
                text=df['Accession_Number'])
fig.add_scatter(x=df[xpos].apply(lambda x: x * -1),
                y=df[ypos],
                mode='markers',
                marker=dict(size=10, color=[color_dict[k] for k in df['All_Pathways'].values]),
                text=df['Accession_Number'])
fig.add_scatter(x=df[xpos].apply(lambda x: x * -1),
                y=df[ypos].apply(lambda x: x * -1),
                mode='markers',
                marker=dict(size=10, color=[color_dict[k] for k in df['All_Pathways'].values]),
                text=df['Accession_Number'])
fig.add_scatter(x=df[xpos],
                y=df[ypos].apply(lambda x: x * -1),
                mode='markers',
                marker=dict(size=10, color=[color_dict[k] for k in df['All_Pathways'].values]),
                text=df['Accession_Number'])
fig.update_layout(height=500, showlegend=False)
    




app.layout = html.Div([
    dcc.Graph(
        id='basic-interactions',
        figure=fig,
        clear_on_unhover=True
    ),

    html.Div(className='row', children=[
        html.Div([
            dcc.Markdown("""
                **Hover Data**

                Mouse over values in the graph.
            """),
            html.Pre(id='hover-data', style=styles['pre'])
        ], className='three columns')
    ]),
 
    #legend?
    html.Div([
        html.Ul([
            html.Li("Integrin", className='circle', style={'background': 'limegreen','color':'black',
                'list-style':'none','text-indent': '17px'}),
            html.Li("Cytoskeleton", className='circle', style={'background': 'orangered','color':'black',
                'list-style':'none','text-indent': '17px','white-space':'nowrap'}),
            html.Li("Blood Coagulation", className='circle', style={'background': 'firebrick','color':'black',
                'list-style':'none','text-indent': '17px','white-space':'nowrap'}),
            html.Li("Chemokine/Cytokine Signaline", className='circle', style={'background': 'tomato','color':'black',
                'list-style':'none','text-indent': '17px','white-space':'nowrap'}),
            html.Li("None reported", className='circle', style={'background': 'slategrey','color':'black',
                'list-style':'none','text-indent': '17px','white-space':'nowrap'}),            
        ], style={'width':'50%','border-bottom': 'solid 3px', 'border-color':'#00FC87','padding-top': '6px'}
        ),

        #Pathway checklist

    html.Div(
        style={'width':'50%', 'height': '100%', 'float':'left'},
        children=[
            dcc.Checklist(id='pathway_name',
                options=[{'label':str(b),'value':b} for b in sorted(df['All_Pathways'].unique())],
                value=[b for b in sorted(df['All_Pathways'].unique())],
            ),
            ],
        
        ),
             ])


    
])


            
#########
#    html.Div([
#        dash_table.DataTable(
#            id='datatable-interactivity',
#            columns=[{'name':i, 'id':i, 'deletable':False, 'selectable':True} for i in df.columns],
#            data=df.to_dict('records'),
#
#            hidden_columns=['Accession_Number','Identified Proteins (551)','ANOVA Test (P-Value','18129_DB_CORO_1',
#                            '18129_DB_CORO_2','18129_DB_CORO_3','18129_DB_CORO_4','18129_DB_CORO','DB_CORO_MEAN',
#                            '18129_HET_CORO_1','18129_HET_CORO_2','18129_HET_CORO_3','18129_HET_CORO_4','HET_CORO_MEAN',
#                            '18129_KOHET_CORO_1','18129_KOHET_CORO_2','18129_KOHET_CORO_3','18129_KOHET_CORO_4',
#                            'KOHET_CORO_MEAN','18129_KOKO_CORO_1','18129_KOKO_CORO_2','18129_KOKO_CORO_3','KOKO_CORO_MEAN',
#                            'Entry','Molecular Weight'],
#            style_cell={'textAlign':'left'},
#            editable=False,
#            sort_action="native",
#            filter_action='native',
#            sort_mode="multi",
#            column_selectable="multi",
#        ),
#        html.Div(id='datatable-interactivity-container')
#    ])
#    
###########

                    

#unfinished callback&def to update checkboxes
@app.callback(Output('graph', 'figure'),
              Input('pathway_name', 'value'))
def scatter_update(chosen_pathway):
    df_sub = df[(df['All_Pathways'].isin(chosen_pathway))] 
    



@app.callback(
    Output('hover-data', 'children'),
    Input('basic-interactions', 'hoverData'))
def display_hover_data(hoverData):
    name = hoverData['points'][0]['text']
    path = df.loc[df['Accession_Number'] == name, 'All_Pathways'].values[0]
    xpos_val = df.loc[df['Accession_Number'] == name, xpos].values[0]
    ypos_val = df.loc[df['Accession_Number'] == name, ypos].values[0]
    xneg_val = df.loc[df['Accession_Number'] == name, xneg].values[0]
    yneg_val = df.loc[df['Accession_Number'] == name, yneg].values[0]
    blah = f'Protein:\t{name}\nPathways:\t{path}\n{xpos}:\t{xpos_val}\n{ypos}:\t{ypos_val}\n{xneg}:\t{xneg_val}\n{yneg}:\t{yneg_val}'
    return blah





if __name__ == '__main__':
    app.run_server(debug=False)
