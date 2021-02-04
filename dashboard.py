#! /usr/bin/env python

# -*- coding: utf-8 -*-                                                                                                                                                                                                                                                      

from collections import defaultdict
import os,sys,re
import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.express as px
import pandas as pd
from dash.dependencies import Input, Output, State
import dash_table
from dash_table.Format import Format, Padding, Scheme
from dash_table import DataTable, FormatTemplate



## importing socket module
import socket
## getting the hostname by socket.gethostname() method
hostname = socket.gethostname()
## getting the IP address using socket.gethostbyname() method
ip_address = socket.gethostbyname(hostname)
## printing the hostname and ip_address
print(f"Hostname: {hostname}:8080")
print(f"IP Address: {ip_address}:8080")


binFilenameSet = set()
filenameSet = set( ['Collection.xlsx' , 'Anvio_summary.tsv' , 'CheckM.tsv' , 'GTDBtk.tsv' , 'Collection.tsv' ] )
directory = '/env/cns/proj/projet_CSD/scratch/assemblies/Ecro_F_AB1/refinedBins/output'
for root, dirs, files in os.walk(directory):
    for filename in files :
        if filename in filenameSet :
            continue
        else:
            binFilenameSet.add(root+'/'+filename)


binNameSet = set()
data = list()
for filename in binFilenameSet :
    binNameSet.add(os.path.basename(filename.replace('.tsv','')))
    file = open(filename,'r')
    headerList = next(file).rstrip().split('\t')
    for line in file :
        line = line.rstrip()
        data.append(line.split('\t'))
    file.close()

options = []
for binName in binNameSet :
    options.append( {'label': binName, 'value': binName })

output = open('/env/cns/proj/agc/home/rmeheust/scripts/bins.info','w')
output.write('\t'.join(headerList)+'\n')
for elt in data :
    output.write('\t'.join(elt)+'\n')
output.close()

# Create the pandas DataFrame 
pd.options.display.float_format = '{:,.2f}'.format
df = pd.read_csv('/env/cns/proj/agc/home/rmeheust/scripts/bins.info',sep="\t")

#print(df.dtypes)

percentage = FormatTemplate.percentage(1)

columns = [
    dict(id='scaffold', name='Scaffold'),
    {'id' : 'bin' , 'name' : 'Bin Name' , 'hideable' : False , 'presentation' : 'dropdown' , 'editable': True},
    dict(id='refineM_outlier', name='Outliers'),
    dict(id='anvio_length', name='Length', type='numeric' , hideable = True),
    dict(id='anvio_gc', name='GC %', type='numeric', format=percentage , hideable = True),
    dict(id='anvio_coverage', name='Anvio Coverage', type='numeric',format=Format(precision=2, scheme=Scheme.fixed), hideable = True),
    dict(id='Mean coverage', name='RefineM Coverage', type='numeric',format=Format(precision=2, scheme=Scheme.fixed), hideable = True), 
    dict(id='anvio_taxonomy', name='Anvio Taxonomy', hideable = True),
    dict(id='phylum: taxa', name='GTDB (phylum)', hideable = True),
    dict(id='class: taxa', name='GTDB (class)', hideable =True),
    dict(id='order: taxa', name='GTDB (order)', hideable = True),
    dict(id='family: taxa', name='GTDB (family)', hideable = True),
    dict(id='genus: taxa', name='GTDB (genus)', hideable = True),
    dict(id='species: taxa', name='GTDB (species)', hideable = True),
]


#print(df.head())


fig = px.scatter(df, x="GC", y="Mean coverage", color="anvio_taxonomy", facet_col="bin", facet_col_wrap=5, hover_name="scaffold", hover_data=["Length (bp)"  , "family: taxa"  , "genus: taxa" , "species: taxa"])







app = dash.Dash(__name__)

app.layout = html.Div([

    html.H1(children='Hello Dash'),
    html.Div(children='Dash: A web application framework for Python.'),


    html.H2(children='Scatter plots'),
    html.Div([

        dcc.Checklist(
            id='xaxis-column',
            options=options,
            value=list(binNameSet)
        ),
        
        dcc.Dropdown(
            id='legend-dropdown',
            options=[
                {'label': 'ANVIO taxonomy', 'value': 'anvio_taxonomy'},
                {'label': 'GTDB domain taxonomy', 'value': 'domain: taxa'},
                {'label': 'GTDB phylum taxonomy', 'value': 'phylum: taxa'},
                {'label': 'GTDB class taxonomy', 'value': 'class: taxa'},
                {'label': 'GTDB order taxonomy', 'value': 'order: taxa'},
                {'label': 'GTDB family taxonomy', 'value': 'family: taxa'},
                {'label': 'GTDB genus taxonomy', 'value': 'genus: taxa'},
                {'label': 'GTDB species taxonomy', 'value': 'species: taxa'}
            ],
            value='anvio_taxonomy'
        ),

        dcc.Graph(id='boxplot-coverage',figure=fig)

    ]),

    html.H2(children='Datatable'),
    html.Div([
        dash_table.DataTable(
            id='datatable-scaffold',
            data=df.to_dict('records'),
            columns=columns,
            style_cell={
                'textAlign': 'left',
            },        
            style_as_list_view=True,
            style_data_conditional=[
                {
                    'if': {'row_index': 'odd'},
                    'backgroundColor': 'rgb(248, 248, 248)'
                }
            ],
            style_table={'height': '400px', 'overflowY': 'auto'},
            dropdown={
                'bin': {
                    'options': [
                        {'label': i, 'value': i}
                        for i in df['bin'].unique()
                ]
                }
            },
            filter_action="native",
            sort_action="native",
            sort_mode="multi",
            page_size=20,  # we have less data in this example, so setting to 20
            fixed_rows={'headers': True},
            style_header={
                'backgroundColor': 'rgb(230, 230, 230)',
                'fontWeight': 'bold',
                'overflow': 'hidden',
            }
        ),
    ]),

    html.Div(id='table-dropdown-container'),

])

@app.callback(
    Output('boxplot-coverage', 'figure'),
    Output('datatable-scaffold', 'data'),
    Input('xaxis-column', 'value'),
    Input('legend-dropdown', 'value')
)


def update_graph_datatable(checkboxList,legendValue) :
    print('\n')
    print('You have selected: '+str(checkboxList))
    print( 'You have selected "{}"'.format(legendValue) )
    print('\n')

    dff = df[ df['bin'].isin(checkboxList) ]
    fig = px.scatter(dff, x="GC", y="Mean coverage", color=legendValue, facet_col="bin", facet_col_wrap=5,hover_name="scaffold", hover_data=["Length (bp)", "family: taxa" , "genus: taxa" , "species: taxa" ])
    return fig,dff.to_dict('records')


if __name__ == '__main__':
    print('test')

    app.run_server(host=ip_address,debug=True, port = 8080)
