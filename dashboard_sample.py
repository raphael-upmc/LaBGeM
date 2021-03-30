#! /usr/bin/env python

# -*- coding: utf-8 -*-                                                                                                                                                                                                                                                      

from collections import defaultdict
import os,sys,re
import dash
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
import plotly.express as px
import pandas as pd
from dash.dependencies import Input, Output, State
import dash_table
from dash_table.Format import Format, Padding, Scheme
from dash_table import DataTable, FormatTemplate


def checkBinName(binName,sample,project) :
    if len(binName) == 0 :
        return False,'Please, provide a name','danger'

    noFunkyCharacter = r'^[A-Za-z0-9_]*$'
    if not re.match(noFunkyCharacter, binName) :
        return False,'Remove the funky characters','danger'


    liste = binName.split('__')
    if len(liste) == 3 :
        name = liste[0]
        p = liste[1]
        s = liste[2]
        if p != project :
            return False,'Not a project name','danger'
            print('error')
        if s != sample :
            return False,'Not a sample name','danger'
    else:
        return False,'Please format the bin name as follows: NAME__PROJECT__SAMPLE','danger'
    
    return True,'The bin has been renamed','success'


def suggestedName(anvio_lineage,gtdb_lineage,project,sample,binName2count) :
    if gtdb_lineage != 'Na' :
        liste = gtdb_lineage.split(';')
        if re.match(r's\_\_',liste[-1]) :
            if liste[-1] == 's__' :
                name = liste[-2].split('__')[1]
            else:
                name = liste[-1].split('__')[1].replace(' ','_')
        else:
            name = liste[-1].split('__')[1]
        name = name.replace('-','_')
    else:
        name = anvio_lineage.replace(' ','_')
        name = name.replace('-','_')
        
    binName = name+'__'+project+'__'+sample
    # check name redundancy
    if binName in binName2count :
        nb = str(binName2count[binName])
        binName2count[binName] += 1
        binName = name+'_'+nb+'__'+project+'__'+sample

    else:
        binName2count[binName] = 1
    return binName

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

port = str('8085')
print("Hostname: "+hostname+":"+port)
print(f"IP Address: {ip_address}:{port}")



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
    gtdb_lineage = liste[14]
    anvio_lineage = liste[1]
    print(gtdb_lineage)
    binName = suggestedName(anvio_lineage,gtdb_lineage,project,sample,binName2count)

    result,msg,color = checkBinName(binName,sample,project)
    if result :
        liste.append(binName)
    else:
        liste.append('TO BE DEFINED')
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


app = dash.Dash(__name__, external_stylesheets=[dbc.themes.YETI]) # https://bootswatch.com/default/

app.layout = html.Div([

    dbc.Row(dbc.Col(html.H1(children='Sample '+sample),
                    width={'size' : '3' , 'offset' : '1', 'order': '1'}
                    )
            ),

    dbc.Row(dbc.Col(html.Div(children='Dash: A web application framework for Python.'))),


    dbc.Row(dbc.Col(html.Hr())),

    dbc.Row(dbc.Col(html.H2(children='Bins summary'))),

    dbc.Row(dbc.Col(html.Div(id='datatable-output-container'))),

    dbc.Row(dbc.Col(html.Hr())),

    dbc.Row(dbc.Col(html.H2(children='Rename bin name'))),

    dbc.Row([
        dbc.Col(
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
                # style={
                #     'width': '300px', 
                #     'font-size': "100%",
                # },
            ),width={'size' : '2' , 'order': '1' , 'offset' : '0'} # 'offset' : '0'
        ),
        dbc.Col(
            dcc.Input(
                id='input-1-state', 
                type='text', 
                value='',
                autoComplete='off',
                debounce=True,
                minLength=10,
                maxLength=50,
                required=True,
                size="50",
                placeholder='suggested name: XXX',
                # style={
                #     'font-size': "100%",
                # },
            ),width={'size' : '3' , 'order': '2'}
        ),
        dbc.Col(
            html.Button(n_clicks=0, children='Submit', id='submit-button-state')
            ,width={'size' : '2' , 'order': '3'}
        )
    ]),

    dbc.Row(dbc.Col(html.Div(id='dd-output-container'))),
    dbc.Row(dbc.Col(html.Div(id='dd-submit-container'),width={'size' : '3' , 'offset' : '0'})),
    dbc.Row(dbc.Col(html.Button(n_clicks=0, children='Check and Save', id='save-button-state'),width={'size' : '3'}))

], style={'marginBottom': 50, 'marginTop': 50 , 'marginLeft' : 25})



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
        return 'Select the bin you want to rename'







@app.callback([Output('datatable-output-container', 'children'),
              Output('dd-submit-container', 'children')],
              Input('submit-button-state', 'n_clicks'),
              State('bin-dropdown', 'value'),
              State('input-1-state', 'value'))


def update_datatable(n_clicks, selectedBin, binName):
    if selectedBin != None :
        isBin = df["Bin"] == selectedBin
        row_index = df.index[isBin].tolist()[0]

        print(isBin)
        print(row_index)
        print(df.head())
        print('input: ')
        print('input 1: '+str(selectedBin))
        print('input 2: '+str(binName))

        result,msg,colorAlert = checkBinName(binName,sample,project) # have to check for redundancy

        if result :
            df.at[row_index,'Name'] = binName
        print(df.head())

        alertButton = dbc.Alert(msg, color=colorAlert)

    else:
        alertButton = ''

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
            style_table={'minWidth': 95 , 'overflowY': 'auto'}, # 'height': '250px',
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
        ),
        alertButton
    ]
    
    
if __name__ == '__main__':
    print('test')

    app.run_server(host=ip_address,debug=True, port = int(port))
