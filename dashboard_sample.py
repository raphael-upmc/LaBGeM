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



def create_conditional_style(df,columns):
    style=[]
    PIXEL_FOR_CHAR = 5
    for elt in columns:
#        print(elt)
        col = elt['id']
        name = elt['name']
#        print(col+'\t'+name)
        name_length = len(name)
        pixel = 70 + round(name_length*PIXEL_FOR_CHAR)
        pixel = str(pixel) + "px"
        style.append({'if': {'column_id': col}, 'minWidth': pixel})
#    print(style)
    return style


## importing socket module
import socket
## getting the hostname by socket.gethostname() method
hostname = socket.gethostname()
## getting the IP address using socket.gethostbyname() method
ip_address = socket.gethostbyname(hostname)
## printing the hostname and ip_address
print("Hostname: "+hostname+":8081")
print(f"IP Address: {ip_address}:8081")



# -------------------------------------------------------------------------------------------------------------------------- #

collection_filename = '/env/cns/proj/projet_CSD/scratch/assemblies/Ecro_F_AB1/refinedBins/output/Collection.tsv'
file = open(collection_filename,'r')
s,project = next(file).strip().split('\t')
print(project)
liste = next(file).strip().split('\t')
sample = liste[1]
liste = next(file).strip().split('\t')
liste = next(file).strip().split('\t')
liste = next(file).strip().split('\t')
next(file)


binName2count = dict()
data = list()
headerList = next(file).rstrip().split('\t')
headerList.append('Name')
for line in file :
    line  = line.rstrip()
    liste = line.split('\t')
    liste[6] = str( float(liste[6])/100.0 )
    liste[7] = str( float(liste[7])/100.0 )
    liste[8] = str( float(liste[8])/100.0 )
    lineage = liste[14]
    print(lineage)
    name = 'TO BE DEFINED'
    print(name)

    if lineage != 'Na' :
        liste1 = lineage.split(';')
        if re.match(r's\_\_',liste1[-1]) :
            if liste1[-1] == 's__' :
                name = liste1[-2].split('__')[1]
            else:
                name = liste1[-1].split('__')[1].replace(' ','_')
        else:
            name = liste1[-1].split('__')[1]
    else:
        name = liste[1].replace(' ','_')
    
        
    binName = name+'__'+project+'__'+sample
    # check name redundancy
    if binName in binName2count :
        nb = str(binName2count[binName])
        binName2count[binName] += 1
        binName = name+'_'+nb+'__'+project+'__'+sample

    else:
        binName2count[binName] = 1

    liste.append(binName)
    data.append( liste )
    
file.close()

output = open('/env/cns/proj/agc/home/rmeheust/scripts/bins.info','w')
output.write('\t'.join(headerList)+'\n')
for elt in data :
    output.write('\t'.join(elt)+'\n')
output.close()

# Create the pandas DataFrame 
pd.options.display.float_format = '{:,.2f}'.format
df = pd.read_csv('/env/cns/proj/agc/home/rmeheust/scripts/bins.info',sep="\t")

# -------------------------------------------------------------------------------------------------------------------------- #






percentage = FormatTemplate.percentage(1)

columns = [
    {'id' : 'Bin' , 'name' : 'id' , 'hideable' : False , 'presentation' : 'dropdown' , 'editable': False},
    {'id' : 'Name' , 'name' : 'Name' , 'hideable' : False , 'presentation' : 'dropdown' , 'editable': True},
    dict(id='Anvio_taxon', name='Anvio taxon', hideable = True),
    dict(id='Anvio_mean_coverage', name='coverage', type='numeric' , format=Format(precision=2, scheme=Scheme.fixed), hideable = True),
    dict(id='Anvio_total_length', name='length', type='numeric' , hideable = True),
    dict(id='Anvio_num_contigs', name='# contigs', type='numeric', hideable = True),
    dict(id='Anvio_N50', name='N50', type='numeric', hideable = True), 
    dict(id='Anvio_GC_content', name='GC%', type='numeric',format=percentage, hideable = True),
    dict(id='Anvio_percent_completion', name='Completeness', type='numeric',format=percentage , hideable = True),
    dict(id='Anvio_percent_redundancy', name='Redundancy', type='numeric',format=percentage , hideable =True), 
    dict(id='Gtdb_classification', name='GTDB classification', hideable = True),
    dict(id='Gtdb_fastani_reference', name='GTDB Reference', hideable = True),
    dict(id='Gtdb_classification_method', name='GTDB classification method', hideable = True),
]


#print(df.head())


app = dash.Dash(__name__)

app.layout = html.Div([

    html.H1(children='Sample '+sample),
    html.Div(children='Dash: A web application framework for Python.'),



    html.H2(children='Bins summary'),
    html.Div(id='datatable-output-container'),

    html.Hr(),

    html.H2(children='Rename bin name'),
    
    html.Div([
        html.Div(
            dcc.Dropdown(
                id='bin-dropdown',
                options = [
                    {'label': i, 'value': i} for i in df['Bin'].unique()
                ],
                placeholder="Select a bin",
                optionHeight=35,
                multi=False,
                searchable=True,
                clearable=True,
                persistence=False,
                persistence_type='memory', # 'memory' (browser tab is refreshed', 'session' (browser tab is closed), 'local' (browser cookies are deleted)
                style={
                    'width': '150px', 
                    'font-size': "100%",
                },
            ),
            style={'display': 'table-cell', 'width':'18%'}),
    
        html.Div(
            dcc.Input(
                id='input-1-state', 
                type='text', 
                value='',
                autoComplete='off',
                debounce=True,
                minLength=10,
                maxLength=50,
                required=True,
                size="20",
                placeholder='suggested name: XXX',
                style={
                    'font-size': "100%",
                },
            ),
        style={'display': 'table-cell', 'width':'50%'}),

        html.Div(
            html.Button(n_clicks=0, children='Submit', id='submit-button-state'),
        style={'display': 'table-cell', 'width':'10%'}),
    ],style={'display': 'table-cell', 'width':'70%'}),
    html.Div(id='dd-output-container')
])



@app.callback(
    dash.dependencies.Output('dd-output-container', 'children'),
    [dash.dependencies.Input('bin-dropdown', 'value')])
def update_output(dd_value):
    if dd_value != None :
        isBin = df["Bin"] == dd_value
        row_index = df.index[isBin].tolist()[0]
        lineage = df.at[row_index,'Gtdb_classification']
        name = ''
        if lineage != 'Na' :
            liste = lineage.split(';')
            if re.match(r's\_\_',liste[-1]) :
                if liste[-1] == 's__' :
                    name = liste[-2].split('__')[1]
                else:
                    name = liste[-1].split('__')[1].replace(' ','_')
            else:
                name = liste[-1].split('__')[1]
        else:
            name = df.at[row_index,'Anvio_taxon'].replace(' ','_')
        return 'Suggested name for '+dd_value+': '+name+'__'+project+'__'+sample
    else:
        return 'Select a bin to get a suggested name'



@app.callback(Output('datatable-output-container', 'children'),
              Input('submit-button-state', 'n_clicks'),
              State('bin-dropdown', 'value'),
              State('input-1-state', 'value'))


def update_datatable(n_clicks, input1, input2):
    if input1 != None :
        isBin = df["Bin"] == input1
        row_index = df.index[isBin].tolist()[0]

        print(isBin)
        print(row_index)
        print(df.head())
        print('input: ')
        print('input 1: '+str(input1))
        print('input 2: '+str(input2))

        # df[row_index]['Name'] = input2
        df.at[row_index,'Name'] = input2
        print(df.head())


    return [
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
            row_selectable='multi',
            selected_rows = [],
            style_cell_conditional = create_conditional_style(df,columns),
            style_table={'minWidth': 95 , 'height': '400px', 'overflowY': 'auto'},
            filter_action="native",
            sort_action="native",
            sort_mode="single",
            page_action = "native",
            page_current=0,
            page_size=20,  # we have less data in this example, so setting to 20
            fixed_rows={'headers': True},
            style_header={
                'backgroundColor': 'rgb(230, 230, 230)',
                'fontWeight': 'bold',
                'overflow': 'hidden',
            }
        )
    ]
    
    
if __name__ == '__main__':
    print('test')

    app.run_server(host=ip_address,debug=True, port = 8081)
