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
import json
import numpy as np
from dash.exceptions import PreventUpdate
import sqlite3
from sqlite3 import Error
import os.path, time


def create_conditional_style(columns):
    style=[]
    PIXEL_FOR_CHAR = 6
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


def create_connection(db_file):
    """ create a database connection to the SQLite database
        specified by db_file
    :param db_file: database file
    :return: Connection object or None
    """
    conn = None
    try:
        conn = sqlite3.connect(db_file , check_same_thread=True)
        return conn
    except Error as e:
        print(e)

    return conn


## importing socket module
import socket
## getting the hostname by socket.gethostname() method
hostname = socket.gethostname()
## getting the IP address using socket.gethostbyname() method
ip_address = socket.gethostbyname(hostname)
## printing the hostname and ip_address

port = str('8080')
print("Hostname: "+hostname+":"+port)
print(f"IP Address: {ip_address}:{port}")



################
# loading data #
################

# tab 1 #

sample = 'Esil_F_AA1'
refinedBin_output_directory = '/env/cns/proj/projet_CSD/scratch/assemblies/'+sample+'/refinedBins/output'
bins_summary_filename = '/env/cns/proj/projet_CSD/scratch/assemblies/'+sample+'/refinedBins/output/Bins_summary.tsv'
refinedBin_directory = '/env/cns/proj/projet_CSD/scratch/assemblies/'+sample+'/refinedBins'


config_filename = '/env/cns/proj/projet_CSD/scratch/assemblies/'+sample+'/assembly/info.json'
with open(config_filename) as f:
    data = json.load(f)

project = data['project']
sample = data['sample']


db_filename = refinedBin_directory+'/'+'REFINEDBINS.db'

percentage = FormatTemplate.percentage(1)

columns_summary = [
    {'id' : 'anvio_id' , 'name' : 'Anvio_id' , 'hideable' : False , 'presentation' : 'dropdown' , 'editable': False},
    {'id' : 'name' , 'name' : 'Name' , 'hideable' : False , 'presentation' : 'dropdown' , 'editable': True},
    dict(id='anvio_taxonomy', name='Anvio taxon', hideable = True),
    dict(id='anvio_coverage', name='coverage', type='numeric' , format=Format(precision=2, scheme=Scheme.fixed), hideable = True),
    dict(id='anvio_length', name='length', type='numeric' , hideable = True),
    dict(id='anvio_contig_nb', name='# contigs', type='numeric', hideable = True),
    dict(id='anvio_N50', name='N50', type='numeric', hideable = True), 
    dict(id='anvio_gc', name='GC%', type='numeric',format=percentage, hideable = True),
    dict(id='anvio_completeness', name='Completeness', type='numeric',format=percentage , hideable = True),
    dict(id='anvio_contamination', name='Redundancy', type='numeric',format=percentage , hideable =True), 
    dict(id='gtdb_taxonomy', name='GTDB classification', hideable = True)
#    dict(id='gtdb_fastani_reference', name='GTDB Reference', hideable = True),
#    dict(id='gtdb_classification_method', name='GTDB classification method', hideable = True),
]


columns_bin = [
    dict(id='scaffold_id', name='Scaffold'),
    {'id' : 'anvio_id' , 'name' : 'Anvio_id' , 'hideable' : False , 'presentation' : 'dropdown' , 'editable': False},
    {'id' : 'anvio_updated_id' , 'name' : 'New Bin Name' , 'hideable' : False , 'presentation' : 'dropdown' , 'editable': True},
    dict(id='refineM_outlier', name='Outliers'),
    dict(id='refineM_length', name='Length', type='numeric' ,format=Format(precision=0, scheme=Scheme.fixed) , hideable = True),
    dict(id='refineM_gc', name='GC %', type='numeric', format=percentage , hideable = True),
    dict(id='anvio_coverage', name='Anvio Coverage', type='numeric',format=Format(precision=1, scheme=Scheme.fixed), hideable = True),
    dict(id='refineM_coverage', name='RefineM Coverage', type='numeric',format=Format(precision=1, scheme=Scheme.fixed), hideable = True), 
    dict(id='anvio_taxonomy', name='Anvio Taxonomy', hideable = True),
    dict(id='refineM_domain', name='GTDB (domain)', hideable = True),
    dict(id='refineM_phylum', name='GTDB (phylum)', hideable = True),
    dict(id='refineM_class', name='GTDB (class)', hideable =True),
    dict(id='refineM_order', name='GTDB (order)', hideable = True),
    dict(id='refineM_family', name='GTDB (family)', hideable = True),
    dict(id='refineM_genus', name='GTDB (genus)', hideable = True),
    dict(id='refineM_species', name='GTDB (species)', hideable = True),
]



##########
# layout #
##########

bins_tab_content = html.Div([
    dbc.Card(
        dbc.CardBody( 
            [
                dbc.Row(dbc.Col(html.Hr())),
                dbc.Row([
                    dbc.Col(html.H2(children='Bins distribution')),
                    dbc.Col(
                        html.Div(id='starting_msg_tab2_id'),width={'size' : '6' , 'order' : '2'}
                    )
                ]),
                dbc.Row(dbc.Col(html.Hr())),
                dbc.Row([
                    dcc.ConfirmDialog(
                        id='confirm-dialog',
                        displayed=False,
                        message='Please choose Bins!',
                    )
                ]),

                dbc.Row(dbc.Col(html.Hr())),

                dbc.Row([
                    dbc.Col(
                        dcc.Dropdown(
                            id='binList_tab2_id',
                            options=[],
                            value= [],
                            multi=True
                        ),width={'size' : '7' , 'order': '1' , 'offset' : '0'} # 'offset' : '0'  
                    ),
                    dbc.Col(
                        dcc.Dropdown(
                            id='legend_dropdown_tab2_id',
                            options=[
                                {'label': 'ANVIO taxonomy', 'value': 'anvio_taxonomy'},
                                {'label': 'GTDB domain taxonomy', 'value': 'refineM_domain'},
                                {'label': 'GTDB phylum taxonomy', 'value': 'refineM_phylum'},
                                {'label': 'GTDB class taxonomy', 'value': 'refineM_class'},
                                {'label': 'GTDB order taxonomy', 'value': 'refineM_order'},
                                {'label': 'GTDB family taxonomy', 'value': 'refineM_family'},
                                {'label': 'GTDB genus taxonomy', 'value': 'refineM_genus'},
                                {'label': 'GTDB species taxonomy', 'value': 'refineM_species'}
                            ],
                            value='anvio_taxonomy'
                        ),width={'size' : '3' , 'order': '2' , 'offset' : '0'} # 'offset' : '0'  
                    ),
                    dbc.Col(
                        dbc.Button(color="primary", outline=True, className="mr-1",n_clicks=0, children='Update scatterplot', id='update_scatter_tab2_id'),width={'order': '3' , 'size' : '2'}
                    )
                ]),

                dbc.Row(dbc.Col(dcc.Graph(id='scatter_tab2_id',style={'width': '100%', 'height': '75vh'}))),

                dbc.Row(dbc.Col(html.Hr())),

                dbc.Row(dbc.Col(
                    dash_table.DataTable(
                        css=[{"selector": ".Select-menu-outer", "rule": "display: block !important"}], # Add this line for dropdown
                        columns=columns_bin,
                        id = 'datatable_scaffold_tab2_id',
                        style_cell={
                            'textAlign': 'left',
                        },        
                        style_as_list_view=True,
                        style_cell_conditional = create_conditional_style(columns_bin),
                        style_data_conditional=[
                            {
                                'if': {'row_index': 'odd'},
                                'backgroundColor': 'rgb(248, 248, 248)'
                            }
                        ],
                        style_table={'height': '800px', 'overflowY': 'auto'},
                        # dropdown={
                        #     'anvio_updated_id': {
                        #         'options': [ {'label': 'test1', 'value': 'test1'}  ]
                        #     }
                        # },
                        filter_action="native",
                        sort_action="native",
                        sort_mode="single",
                        page_size=20,  # we have less data in this example, so setting to 20
                        fixed_rows={'headers': True},
                        style_header={
                            'backgroundColor': 'rgb(230, 230, 230)',
                            'fontWeight': 'bold',
                            'overflow': 'hidden',
                        }
                    )
                )
            ),
            dbc.Row(dbc.Col(html.Hr()))
            ]
        ),
        className="mt-3",style={'width':'150%'}
    )
])


summary_tab_content = html.Div([
    dbc.Card(
        dbc.CardBody( 
            [
                dbc.Row([
                    dbc.Col(html.H2(children='Bins summary')),
                    dbc.Col(
                        html.Div(id='starting_msg_tab1_id'),width={'size' : '6' , 'order' : '2'}
                    )
                ]),

                dbc.Row(dbc.Col(html.Hr())),
                
                dbc.Row(dbc.Col(html.Hr())),
                dbc.Row(dbc.Col(
                    dash_table.DataTable(
                        id='bins_datatable_tab1_id',
                        columns=columns_summary,
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
                        style_cell_conditional = create_conditional_style(columns_summary),
                        style_table={'minWidth': '100%' , 'overflowY': 'auto' , }, # 'height': '250px',
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
                )),


                dbc.Row(dbc.Col(html.Hr())),


                dbc.Row(dbc.Col(html.H2(children='Update bin name'))), # update bin name                
                dbc.Row([
                    dbc.Col(
                        dcc.Dropdown(
                            id='bin_dropdown_tab1_id',
                            options = [],
                            placeholder="Select a bin",
                            optionHeight=35,
                            multi=False,
                            searchable=True,
                            clearable=True,
                            persistence=False,
                            #persistence_type='session', # 'memory' (browser tab is refreshed', 'session' (browser tab is closed), 'local' (browser cookies are deleted)
                        ),width={'size' : '3' , 'order': '1' , 'offset' : '0'} # 'offset' : '0'
                    ),
                    dbc.Col(
                        dcc.Input(
                            id='input_tab1_id', 
                            type='text', 
                            value='',
                            autoComplete='off',
                            debounce=True,
                            minLength=10,
                            maxLength=50,
                            required=True,
                            #bs_size="lg",
                            size="60"
                        ),width={'size' : '4' , 'order': '2'}
                    ),
                    dbc.Col(
                        dbc.Button(color="primary", outline=False, className="mr-1",n_clicks=0, children='Update', id='submit_bin_name_tab1_id')
                        ,width={'order': '3'}
                    ),
                    dbc.Col(
                        html.Div(id='submit_new_bin_name_msg_tab1_id'),width={'size' : '3' , 'order' : '4'}
                    )
                ]),
                dbc.Row(dbc.Col(html.Div(id='suggested_bin_name_tab1_id'))),

                dbc.Row(dbc.Col(html.Hr())),

                dbc.Row(dbc.Col(html.H2(children='Add a new bin'))), # update bin name                
                dbc.Row([
                    dbc.Col(
                        dcc.Input(
                            id='input_add_bin_tab1_id', 
                            type='text', 
                            value='',
                            autoComplete='off',
                            debounce=True,
                            minLength=10,
                            maxLength=50,
                            required=True,
                            #bs_size="lg",
                            size="60"
                        ),width={'size' : '4' , 'order': '2'}
                    ),
                    dbc.Col(
                        dbc.Button(color="primary", outline=False, className="mr-1",n_clicks=0, children='Add', id='add_new_bin_name_tab1_id')
                        ,width={'order': '3'}
                    ),
                    dbc.Col(
                        html.Div(id='add_new_bin_name_msg_tab1_id'),width={'size' : '3' , 'order' : '4'}
                    )
                ]),
                dbc.Row(dbc.Col(html.Hr())),



            ]
        ),style={"width": "150%" , "className":"mt-3"},
    )
])



app = dash.Dash(__name__, external_stylesheets=[dbc.themes.YETI]) # https://bootswatch.com/default/

app.layout = dbc.Container(
    [

    dbc.Row(dbc.Col(html.H1(children = 'Sample '+sample),width={'size' : '6' , 'offset' : '0', 'order': '1'})),
    dbc.Row(dbc.Col(html.Div(children='Dash: A web application framework for Python.'))),

    dbc.Tabs(
        [
            dbc.Tab(summary_tab_content,label="Summary", tab_id="summary"),
            dbc.Tab(bins_tab_content,label="Bins", tab_id="bins")
        ],
        id="tabs",
        active_tab="summary",
    ),

    ], style={'marginBottom': 50, 'marginTop': 25 , 'marginLeft' : 0 , 'marginRight' : 0}
)





##################
# callbacks tab1 #
##################

print('Callback...')


@app.callback([Output('bins_datatable_tab1_id', 'data'),
               Output('bin_dropdown_tab1_id','options'),
               Output('submit_new_bin_name_msg_tab1_id','children')],
              [Input('submit_bin_name_tab1_id', 'n_clicks')], #, prevent_initial_call=True)
              [State('bin_dropdown_tab1_id', 'value'),
              State('input_tab1_id', 'value')]
          )

def populate_datatable(n_clicks_update,anvio_id,name):
    print('\n\npopulate datatable (TAB1)')
    print('connecting to the sqlite db done...')
    conn = sqlite3.connect(db_filename , check_same_thread=True)
    update_msg = ''

    input_triggered = dash.callback_context.triggered[0]["prop_id"].split(".")[0]
    print(input_triggered)
    if input_triggered == "submit_bin_name_tab1_id": 
        print('update '+db_filename+' with a new bin name '+name+' ('+anvio_id+')')
        sql_update_query = """Update bins set name = ? where anvio_id = ?"""
        data = (name , anvio_id)
        cursor = conn.cursor()
        cursor.execute(sql_update_query, data)
        cursor.close()
        print(conn.commit())
        print("Record Updated successfully")
        update_msg = "Record Updated successfully"


    print('displaying the datatable....')
    query = 'SELECT anvio_id , name , anvio_length , anvio_gc , anvio_contig_nb , anvio_N50 , anvio_completeness , anvio_contamination, anvio_coverage , anvio_taxonomy , gtdb_taxonomy FROM bins'
    print(query)
    df = pd.read_sql_query(query, con=conn)
    print(df.head())
    conn.close()
    print('close connection')
    return df.to_dict('records') , [ {'label': i, 'value': i} for i in df['anvio_id'].unique() ] , update_msg



@app.callback(
    [Output('suggested_bin_name_tab1_id', 'children'),
     Output('input_tab1_id', 'value')],
    Input('bin_dropdown_tab1_id', 'value'))

def suggestedBinName(dd_value):
    print('\n\nSuggested Bin Name (TAB1)')
    if dd_value != None :
        conn = sqlite3.connect(db_filename , check_same_thread=True)
        query = 'SELECT anvio_id , anvio_taxonomy , gtdb_taxonomy FROM bins'
        df = pd.read_sql_query(query, con=conn)
        conn.close()
        isBin = df["anvio_id"] == dd_value
        row_index = df.index[isBin].tolist()[0]
        lineage = df.at[row_index,'gtdb_taxonomy']
        name = ''
        if lineage != 'NULL' :
            liste = lineage.split(';')
            if re.match(r's\_\_',liste[-1]) :
                if liste[-1] == 's__' :
                    name = liste[-2].split('__')[1]
                else:
                    name = liste[-1].split('__')[1].replace(' ','_')
            else:
                name = liste[-1].split('__')[1]
        else:
            name = df.at[row_index,'anvio_taxonomy'].replace(' ','_')
        return 'Suggested name for '+dd_value+': '+name+'__'+project+'__'+sample,''
    else:
        return 'Select the bin you want to rename',''



##################
# callbacks tab2 #
##################

@app.callback([Output('binList_tab2_id', 'options'),
               Output('binList_tab2_id', 'value'),
               Output('starting_msg_tab2_id', 'children')],
              [Input('add_new_bin_name_tab1_id', 'n_clicks')] #, prevent_initial_call=True)
          )

def initBinDropdown(binList):
    print('\n\npopulate bin dropdown (TAB2)')
    print('connecting to the sqlite db done...')
    conn = sqlite3.connect(db_filename , check_same_thread=True)

    input_triggered = dash.callback_context.triggered[0]["prop_id"].split(".")[0]
    if input_triggered == "add_new_bin_name_tab1_id" :
        msg = 'Updating the bins list'
    else:
        msg = 'Initializing the bins list'

    query = 'SELECT anvio_id FROM bins'
    print(query)
    df = pd.read_sql_query(query, con=conn)
    print(df.head())
    conn.close()
    print('close connection')
    options_binList = [ {'label': i, 'value': i} for i in df['anvio_id'].unique() ]
    options_binList.append({'label': 'Unbinned' , 'value': 'Unbinned'})
    print(options_binList)
    return  options_binList , df['anvio_id'].unique() , msg 




@app.callback([Output('datatable_scaffold_tab2_id', 'data'),
               Output(component_id = 'datatable_scaffold_tab2_id', component_property = 'dropdown')],
               [Input('binList_tab2_id', 'value'), #, prevent_initial_call=True)
               Input('update_scatter_tab2_id', 'n_clicks')], #, prevent_initial_call=True)
           )

def populate_datatable(binList , n_clicks_update):
    print('\n\npopulate scaffold datatable (TAB2)')
    print('connecting to the sqlite db done...')

    conn = sqlite3.connect(db_filename , check_same_thread=True)

    input_triggered = dash.callback_context.triggered[0]["prop_id"].split(".")[0]
    print(input_triggered)
    if input_triggered == "update_scatter_tab2_id": 
        print('update '+db_filename+' with a new scaffold assignements')
        # sql_update_query = """Update bins set name = ? where anvio_id = ?"""
        # data = (name , anvio_id)
        # cursor = conn.cursor()
        # cursor.execute(sql_update_query, data)
        # cursor.close()
        # print(conn.commit())
        print("Record Updated successfully")


    print('displaying the datatable....')
    query = 'SELECT scaffold_id , anvio_id , anvio_updated_id , refineM_outlier , refineM_length , refineM_gc , anvio_coverage , refineM_coverage , anvio_taxonomy , refineM_domain , refineM_phylum , refineM_class , refineM_order , refineM_family , refineM_genus , refineM_species FROM scaffolds WHERE anvio_updated_id IN ( \''+'\' , \''.join(binList)+'\' )'
    print(query)
    df = pd.read_sql_query(query, con=conn)
    print(df.head())
    conn.close()
    print('close connection')
    dropdown={
        'anvio_updated_id': {
            'options': [ 
                {'label': i, 'value': i} for i in df['anvio_id'].unique()
            ]
        }
    }
    dropdown['anvio_updated_id']['options'].append( { 'label' : 'Unbinned' , 'value' : 'Unbinned' } )
    return df.to_dict('records') , dropdown


@app.callback(
    [Output('scatter_tab2_id', 'figure')],
    [Input('legend_dropdown_tab2_id', 'value'),
    Input('datatable_scaffold_tab2_id', 'data')]
)

def display_graph(legendValue,dataset) :
    print('\n\ndisplay_graph (TAB2)')
    print(legendValue)
    print(len(pd.DataFrame(dataset)))
    fig = px.scatter(pd.DataFrame(dataset), x="refineM_gc", y="refineM_coverage", color=legendValue, facet_col="anvio_updated_id", facet_col_wrap=3,hover_name="scaffold_id", hover_data=["refineM_length", "refineM_class", "refineM_order", "refineM_family" , "refineM_genus" , "refineM_species" ])#,height=800)

    return [fig]


if __name__ == '__main__':
    app.run_server(host=ip_address,debug=True, port = int(port))

