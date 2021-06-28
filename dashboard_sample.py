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




def checkUniqness(df,row_index,binName) :
    #df.at[row_index,'Name'] = binName
    name2count = defaultdict(int)
    for name in df['Name'] :
        name2count[name] += 1
    
    if binName in name2count :
        return False,'This bin name already exists','danger'
    else:
        return True,'The bin has been renamed','success'


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
            if liste[-1] == 's__' : # if gtdb species name is empty, look for genus, family, order, class, phylum, domain name...
                name = 'Unknown'
                print(reversed(liste))
                for taxa in reversed(liste) :
                    if taxa.split('__')[1] != '' :
                        name = taxa.split('__')[1]
                        break
                    else:
                        continue
            else: # if has a gtdb species name then concatenate genus and species name
                name = liste[-1].split('__')[1].replace(' ','_')
        else: # no species s__ field...
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


def load_bins_summary_from_saved_file(bins_summary_filename,bin2name_filename):
    binName2count = dict()
    data = list()

    bin2name = dict()
    file = open(bin2name_filename)
    header = next(file)
    for line in file :
        line = line.rstrip()
        binName,name = line.split('\t')
        bin2name[ binName ] = name
    file.close()

    file = open(bins_summary_filename,'r')
    headerList = next(file).rstrip().split('\t')
    headerList.append('Name')
    for line in file :
        line  = line.rstrip()
        liste = line.split('\t')
        liste[2] = float(liste[2])
        liste[3] = int(liste[3])
        liste[4] = int(liste[4])
        liste[5] = int(liste[5])
        liste[6] = float(liste[6])/100.0
        liste[7] = float(liste[7])/100.0
        liste[8] = float(liste[8])/100.0
        gtdb_lineage = liste[14]
        anvio_lineage = liste[1]
        
        liste.append(bin2name[ liste[0] ])
        data.append( liste )
    file.close()

    df_summary = pd.DataFrame( data, columns = headerList )
    
    return df_summary, binName2count,'Bin names have been loaded from a previous saved work'

def load_bins_summary_from_scratch(bins_summary_filename) :
    binName2count = dict()
    data = list()

    file = open(bins_summary_filename,'r')
    headerList = next(file).rstrip().split('\t')
    headerList.append('Name')
    for line in file :
        line  = line.rstrip()
        liste = line.split('\t')
        liste[2] = float(liste[2])
        liste[3] = int(liste[3])
        liste[4] = int(liste[4])
        liste[5] = int(liste[5])
        liste[6] = float(liste[6])/100.0
        liste[7] = float(liste[7])/100.0
        liste[8] = float(liste[8])/100.0
        gtdb_lineage = liste[14]
        anvio_lineage = liste[1]
        binName = suggestedName(anvio_lineage,gtdb_lineage,project,sample,binName2count)
        result,msg,color = checkBinName(binName,sample,project)
        if result :
            liste.append(binName)
        else:
            liste.append('TO BE DEFINED')
        data.append( liste )
    file.close()

    df_summary = pd.DataFrame( data, columns = headerList )
    
    return df_summary, binName2count,'Suggested bin names, please check the names before saving them'


def load_scaffold_from_scratch(directory):
    binFilenameSet = set()
    filenameSet = set( [ 'Collection.xlsx' , 'Anvio_summary.tsv' , 'CheckM.tsv' , 'GTDBtk.tsv' , 'Bins_summary.tsv' , 'Sample_summary.tsv' , 'partial_scaffolds.tsv' ] )
    for root, dirs, files in os.walk(directory):
        for filename in files :
            if filename in filenameSet :
                continue
            else:
                binFilenameSet.add(root+'/'+filename)


    binNameSet = set()
    data_bin = list()
    for filename in binFilenameSet :
        lengthSet = set()
        binNameSet.add(os.path.basename(filename.replace('.tsv','')))
        file = open(filename,'r')
        headerList = next(file).rstrip().split('\t')
        for line in file :
            line = line.rstrip()
            liste = line.split('\t')
            sublist = []
            for i in [0,1,2,3,4,6,7,8,9,10,13,18,23,28,33,38,43]:
                #            print(i)
                if i >= len(liste) :
                    sublist.append(np.nan)
                else:
                    if i == 3 or i == 8 : # int
                        if liste[i] == 'Na' :
                            sublist.append(np.nan)
                        else:
                            sublist.append(int(liste[i]))
                    elif i == 6 or i == 4 or i == 9 or i == 10 : # float
                        if liste[i] == 'Na' :
                            sublist.append(np.nan)
                        else:
                            sublist.append(float(liste[i]))
                    else:
                        sublist.append(liste[i])
            data_bin.append(sublist)
            lengthSet.add(len(line.split('\t')))
        file.close()
    if len(lengthSet) != 1:
        print('ERROR DURING REFININGBINS, CONTACT RAPHAEL ('+filename+') '+str(lengthSet))

    headerSublist = []
    for i in [0,1,2,3,4,6,7,8,9,10,13,18,23,28,33,38,43]:
        headerSublist.append(headerList[i])

    df_bin = pd.DataFrame( data_bin, columns = headerSublist )
    #print(df_bin.head())
    return df_bin


def load_scaffold_from_saved_file(scaffold2bin_filename,directory):
    
    scaffold2bin = dict()
    file = open(scaffold2bin_filename)
    header = next(file)
    for line in file :
        line = line.rstrip()
        scaffold,binName = line.split('\t')
        scaffold2bin[scaffold] = binName
    file.close()

    binFilenameSet = set()
    filenameSet = set( [ 'Collection.xlsx' , 'Anvio_summary.tsv' , 'CheckM.tsv' , 'GTDBtk.tsv' , 'Bins_summary.tsv' , 'Sample_summary.tsv' , 'partial_scaffolds.tsv' ] )
    for root, dirs, files in os.walk(directory):
        for filename in files :
            if filename in filenameSet :
                continue
            else:
                binFilenameSet.add(root+'/'+filename)

    binNameSet = set()
    data_bin = list()
    for filename in binFilenameSet :
        lengthSet = set()
        binNameSet.add(os.path.basename(filename.replace('.tsv','')))
        file = open(filename,'r')
        headerList = next(file).rstrip().split('\t')
        for line in file :
            line = line.rstrip()
            liste = line.split('\t')
            sublist = []
            for i in [0,1,2,3,4,6,7,8,9,10,13,18,23,28,33,38,43]:
                #            print(i)
                if i >= len(liste) :
                    sublist.append(np.nan)
                else:
                    if i == 3 or i == 8 : # int
                        if liste[i] == 'Na' :
                            sublist.append(np.nan)
                        else:
                            sublist.append(int(liste[i]))
                    elif i == 6 or i == 4 or i == 9 or i == 10 : # float
                        if liste[i] == 'Na' :
                            sublist.append(np.nan)
                        else:
                            sublist.append(float(liste[i]))
                    else:
                        sublist.append(liste[i])
            sublist[1] = scaffold2bin[sublist[0]]
            data_bin.append(sublist)
            lengthSet.add(len(line.split('\t')))
        file.close()
    if len(lengthSet) != 1:
        print('ERROR DURING REFININGBINS, CONTACT RAPHAEL ('+filename+') '+str(lengthSet))

    headerSublist = []
    for i in [0,1,2,3,4,6,7,8,9,10,13,18,23,28,33,38,43]:
        headerSublist.append(headerList[i])

    df_bin = pd.DataFrame( data_bin, columns = headerSublist )
    #print(df_bin.head())
    return df_bin



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

sample = sys.argv[1]
refinedBin_output_directory = '/env/cns/proj/projet_CSD/scratch/assemblies/'+sample+'/refinedBins/output'
bins_summary_filename = '/env/cns/proj/projet_CSD/scratch/assemblies/'+sample+'/refinedBins/output/Bins_summary.tsv'
refinedBin_directory = '/env/cns/proj/projet_CSD/scratch/assemblies/'+sample+'/refinedBins'


config_filename = '/env/cns/proj/projet_CSD/scratch/assemblies/'+sample+'/assembly/info.json'
with open(config_filename) as f:
    data = json.load(f)

project = data['project']
sample = data['sample']


bin2name_filename = refinedBin_directory+'/'+"bin2name.tsv"
if not os.path.exists(bin2name_filename):
    df_summary,binName2count,starting_msg = load_bins_summary_from_scratch(bins_summary_filename)
else:
    df_summary,binName2count,starting_msg = load_bins_summary_from_saved_file(bins_summary_filename,bin2name_filename)


print(df_summary.dtypes)

percentage = FormatTemplate.percentage(1)

columns_summary = [
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



# tab 2 

scaffold2bin_filename = refinedBin_directory+'/'+"scaffold2bin.tsv"
print('loading scaffold data...')
df_bin = load_scaffold_from_scratch(refinedBin_output_directory)

print(df_bin.dtypes)
print(df_bin.head())


options_checklist = []
for binName in df_summary['Bin'].unique() :
    options_checklist.append( {'label': binName, 'value': binName })
options_checklist.append( {'label': 'Unbinned', 'value': 'Unbinned' })


percentage = FormatTemplate.percentage(1)
columns_bin = [
    dict(id='scaffold', name='Scaffold'),
    {'id' : 'bin' , 'name' : 'Bin Name' , 'hideable' : False , 'presentation' : 'dropdown' , 'editable': False},
    {'id' : 'new_bin' , 'name' : 'New Bin Name' , 'hideable' : False , 'presentation' : 'dropdown' , 'editable': True},
    dict(id='refineM_outlier', name='Outliers'),
    dict(id='refineM_length', name='Length', type='numeric' ,format=Format(precision=0, scheme=Scheme.fixed) , hideable = True),
    dict(id='refineM_gc', name='GC %', type='numeric', format=percentage , hideable = True),
    dict(id='anvio_coverage', name='Anvio Coverage', type='numeric',format=Format(precision=1, scheme=Scheme.fixed), hideable = True),
    dict(id='refineM_coverage', name='RefineM Coverage', type='numeric',format=Format(precision=1, scheme=Scheme.fixed), hideable = True), 
    dict(id='anvio_taxonomy', name='Anvio Taxonomy', hideable = True),
    dict(id='phylum: taxa', name='GTDB (phylum)', hideable = True),
    dict(id='class: taxa', name='GTDB (class)', hideable =True),
    dict(id='order: taxa', name='GTDB (order)', hideable = True),
    dict(id='family: taxa', name='GTDB (family)', hideable = True),
    dict(id='genus: taxa', name='GTDB (genus)', hideable = True),
    dict(id='species: taxa', name='GTDB (species)', hideable = True),
]



##########
# layout #
##########

summary_tab_content = dbc.Card(
    dbc.CardBody( 
        [
            dbc.Row([
                dbc.Col(html.H2(children='Bins summary')),
                dbc.Col(
                    dbc.Alert(starting_msg,color="success",duration=10000)
                )
            ]),
            dbc.Row(dbc.Col(html.Hr())),
            dbc.Row([
                dbc.Col(

                    dbc.ButtonGroup(
                        [dbc.Button(children="Load saved names",id='load_saved_names_id',n_clicks=0,outline=False,color="primary"), dbc.Button("Restart from scratch",id='restart_from_scratch_id',n_clicks=0,outline=False,color="primary"),dbc.Button(color="primary", outline=False,n_clicks=0, children='Save the bin names', id='save_to_csv_id')], 
                        size="lg",
                        className="mr-1",
                    ),width={'size' : '5' , 'order': '1' , 'offset' : '0'} # 'offset' : '0'
                ),
                dbc.Col(
                    html.Div(id='button_group_msg_id'),width={'size' : '6' , 'order' : '2'}
                )
            ]),

            dbc.Row(dbc.Col(html.Hr())),
            dbc.Row(dbc.Col(
                dash_table.DataTable(
                    id='datatable_bin_summary',
                    data=df_summary.to_dict('records'),
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
                    row_selectable='multi',
                    selected_rows = [],
                    style_cell_conditional = create_conditional_style(df_summary,columns_summary),
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
            dbc.Row(dbc.Col(html.H2(children='Rename bin name'))),
            dbc.Row([
                dbc.Col(
                    dcc.Dropdown(
                        id='bin-dropdown',
                        options = [
                            {'label': i, 'value': i} for i in df_summary['Bin'].unique()
                        ],
                        placeholder="Select a bin",
                        optionHeight=35,
                        multi=False,
                        searchable=True,
                        clearable=True,
                        persistence=False,
                        persistence_type='memory', # 'memory' (browser tab is refreshed', 'session' (browser tab is closed), 'local' (browser cookies are deleted)
                    ),width={'size' : '3' , 'order': '1' , 'offset' : '0'} # 'offset' : '0'
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
                        #bs_size="lg",
                        size="60"
                    ),width={'size' : '4' , 'order': '2'}
                ),
                dbc.Col(
                    dbc.Button(color="primary", outline=False, className="mr-1",n_clicks=0, children='Submit', id='submit_new_bin_name')
                    ,width={'order': '3'}
                ),
                dbc.Col(
                    html.Div(id='submit_new_bin_name_msg'),width={'size' : '3' , 'order' : '4'}
                )
            ]),
            dbc.Row(dbc.Col(html.Div(id='suggested_bin_name'))),
            dbc.Row(dbc.Col(html.Hr()))
        ]
    ),
    style={"width": "150%" , "className":"mt-3"},
)


bins_tab_content = html.Div([
    dcc.Store(id='memory_output_tab2_id',storage_type='session'),    
    dbc.Card(
        dbc.CardBody( 
            [
                dbc.Row(dbc.Col(html.Hr())),
                dbc.Row(dbc.Col(html.H2(children='Bins distribution'))),
                dbc.Row(dbc.Col(html.Hr())),
                dbc.Row([
                    dcc.ConfirmDialog(
                        id='confirm-dialog',
                        displayed=False,
                        message='Please choose Bins!',
                    )
                ]),

                dbc.Row([
                    dbc.Col(
                        
                        dbc.ButtonGroup(
                            [dbc.Button(children="Load",id='load_saved_names_tab2_id',n_clicks=0,outline=False,color="primary"), dbc.Button("Restart from scratch",id='restart_from_scratch_tab2_id',n_clicks=0,outline=False,color="primary"),dbc.Button(color="primary", outline=False,n_clicks=0, children='Save', id='save_to_csv_tab2_id')], 
                            size="lg",
                            className="mr-1",
                        ),width={'size' : '5' , 'order': '1' , 'offset' : '0'} # 'offset' : '0'
                    ),
                    dbc.Col(
                        html.Div(id='button_group_msg_tab2_id'),width={'size' : '6' , 'order' : '2'}
                    )
                ]),
                dbc.Row(dbc.Col(html.Hr())),

                dbc.Row([
                    dbc.Col(
                        dcc.Dropdown(
                            id='binList_tab2_id',
                            options=options_checklist,
                            value= df_summary['Bin'].unique(),
                            multi=True
                        ),width={'size' : '7' , 'order': '1' , 'offset' : '0'} # 'offset' : '0'  
                    ),
                    dbc.Col(
                        dcc.Dropdown(
                            id='legend_dropdown_tab2_id',
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
                        data=df_bin.to_dict('records'),
                        css=[{"selector": ".Select-menu-outer", "rule": "display: block !important"}], # Add this line for dropdown
                        columns=columns_bin,
                        id = 'datatable_scaffold_tab2_id',
                        style_cell={
                            'textAlign': 'left',
                        },        
                        style_as_list_view=True,
                        style_cell_conditional = create_conditional_style(df_bin,columns_bin),
                        style_data_conditional=[
                            {
                                'if': {'row_index': 'odd'},
                                'backgroundColor': 'rgb(248, 248, 248)'
                            }
                        ],
                        style_table={'height': '800px', 'overflowY': 'auto'},
                        dropdown={
                            'new_bin': {
                                'options': [
                                    {'label': i, 'value': i} for i in df_bin['bin'].unique()
                                ]
                            }
                        },
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


app = dash.Dash(__name__, external_stylesheets=[dbc.themes.YETI]) # https://bootswatch.com/default/

app.layout = dbc.Container(
    [

    dbc.Row(dbc.Col(html.H1(children = 'Sample '+sample),width={'size' : '6' , 'offset' : '0', 'order': '1'})),
    dbc.Row(dbc.Col(html.Div(children='Dash: A web application framework for Python.'))),

    dbc.Tabs(
        [
            dbc.Tab(summary_tab_content,label="Summary", tab_id="summary"),
            dbc.Tab(bins_tab_content,label="Bins", tab_id="bins"),
        ],
        id="tabs",
        active_tab="summary",
    ),

    ], style={'marginBottom': 50, 'marginTop': 25 , 'marginLeft' : 0 , 'marginRight' : 0}
)





##################
# callbacks tab1 #
##################



@app.callback(
    [Output('suggested_bin_name', 'children'),
     Output('input-1-state', 'value')],
    Input('bin-dropdown', 'value'))

def update_output(dd_value):
    if dd_value != None :
        isBin = df_summary["Bin"] == dd_value
        row_index = df_summary.index[isBin].tolist()[0]
        lineage = df_summary.at[row_index,'Gtdb_classification']
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
            name = df_summary.at[row_index,'Anvio_taxon'].replace(' ','_')
        return 'Suggested name for '+dd_value+': '+name+'__'+project+'__'+sample,''
    else:
        return 'Select the bin you want to rename',''


@app.callback([Output('datatable_bin_summary', 'data'),
               Output('button_group_msg_id', 'children'),
               Output('submit_new_bin_name_msg', 'children')],
              [Input('save_to_csv_id', 'n_clicks'),
               Input('load_saved_names_id', 'n_clicks'),
               Input('restart_from_scratch_id', 'n_clicks'),
               Input('submit_new_bin_name', 'n_clicks')],
              [State('bin-dropdown', 'value'),
               State('datatable_bin_summary', 'data'),
               State('input-1-state', 'value')])


def update_datatable(n_clicks_save,n_clicks_load,n_clicks_restart,n_clicks_rename_bin, selectedBin, dataset, binName):

    input_triggered = dash.callback_context.triggered[0]["prop_id"].split(".")[0]

    if input_triggered == "save_to_csv_id":
        print("The data has been saved to your folder.")
        df = pd.DataFrame(dataset)
        df[["Bin", "Name"]].to_csv(bin2name_filename,sep='\t',index=False)
        return [dash.no_update,dbc.Alert("The data has been saved in "+bin2name_filename+".",color='success',duration=3000),'']

    elif input_triggered == "load_saved_names_id": # load_bins_summary_from_saved_file
        print("Previous bin names have been erased so you can start from scratch.")
        df_summary_tmp,binName2count,starting_msg = load_bins_summary_from_saved_file(bins_summary_filename,bin2name_filename)
        return [df_summary_tmp.to_dict('records'),dbc.Alert("The data has been loaded from previous work.",color='success',duration=3000),'']

    elif input_triggered == "restart_from_scratch_id":
        print("The data has been loaded from scratch.")
        df_summary_tmp,binName2count,starting_msg = load_bins_summary_from_scratch(bins_summary_filename)
        print(df_summary_tmp.head())
        return [df_summary_tmp.to_dict('records'),dbc.Alert("The data has been loaded from scratch.",color='success',duration=3000),'']

    elif input_triggered == "submit_new_bin_name":

        if selectedBin != None :
            df = pd.DataFrame(dataset)
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
                # check for uniqness
                result2,msg2,colorAlert2 = checkUniqness(df,row_index,binName)
                if result2 :
                    df.at[row_index,'Name'] = binName
                    #print(df_summary.head())
                alertButton = dbc.Alert(msg2, color=colorAlert2,duration=3000)
            else:
                alertButton = dbc.Alert(msg, color=colorAlert,duration=3000)
        else:
            alertButton = ''
        return [df.to_dict('records'),alertButton,'']
    else:
        print('initial call')
        raise PreventUpdate



    
##################
# callbacks tab2 #
##################


# https://dash.plotly.com/datatable/interactivity
# https://community.plotly.com/t/dash-table-datatable-filtering-and-sorting-doesnt-seem-to-work/16362

@app.callback(
    [Output('memory_output_tab2_id', 'data')],
    [Input('update_scatter_tab2_id','n_clicks'),
     Input('load_saved_names_tab2_id','n_clicks'),
     Input('restart_from_scratch_tab2_id','n_clicks'),
     Input('save_to_csv_tab2_id','n_clicks')],
    [State('memory_output_tab2_id', 'data'),
     State('datatable_scaffold_tab2_id', 'data')]
)

def update_store(n_clicks_update,n_clicks_load,n_clicks_restart, n_clicks_save,data_scaffold2bin,dataset):
    print('\n\nupdate dcc.store (TAB2)')
    input_triggered = dash.callback_context.triggered[0]["prop_id"].split(".")[0]
    print(input_triggered)
    if input_triggered == "update_scatter_tab2_id": # updating df_bin
        print('updating the dcc.store component...')
        df_scaffold2bin = pd.DataFrame(data_scaffold2bin)
        cpt = 0
        for row in dataset :
            cpt += 1
            isScaffold = df_scaffold2bin["scaffold"] == row['scaffold']
            row_index = df_scaffold2bin.index[isScaffold].tolist()[0]
            if df_scaffold2bin.at[row_index,'new_bin'] != row['new_bin'] :
                print(row['scaffold']+'\t'+df_scaffold2bin.at[row_index,'new_bin']+'\t==>\t'+row['new_bin'])    
                df_scaffold2bin.at[row_index,'new_bin'] = row['new_bin']

    elif input_triggered == "load_saved_names_tab2_id": # updating df_bin
        print('loading the saved scaffold2bin...')
        df_scaffold2bin = pd.read_csv(scaffold2bin_filename,sep='\t')
    elif input_triggered == "restart_from_scratch_tab2_id": # updating df_bin
        print('starting from scratch')
        df_scaffold2bin = pd.DataFrame( data = { 'scaffold' : df_bin['scaffold'] , 'new_bin' : df_bin['bin'] } )

    elif input_triggered == "save_to_csv_tab2_id": # updating df_bin
        print('saving scaffold2bin...')
        df_scaffold2bin = pd.DataFrame(data_scaffold2bin)
        df_scaffold2bin.to_csv(scaffold2bin_filename,sep='\t',index=False)
    else:
        print('initial call')
        if os.path.exists(scaffold2bin_filename) :
            print(scaffold2bin_filename+' exists')
            df_scaffold2bin = pd.read_csv(scaffold2bin_filename,sep='\t')
            
        else:
            print(scaffold2bin_filename+' does not exist')
            df_scaffold2bin = pd.DataFrame( data = { 'scaffold' : df_bin['scaffold'] , 'new_bin' : df_bin['bin'] } )
            df_scaffold2bin.to_csv(scaffold2bin_filename,sep='\t',index=False)

    print('dcc end\n\n\n')
    return [df_scaffold2bin.to_dict('records')]



@app.callback(
    [Output('datatable_scaffold_tab2_id', 'data'),
     Output(component_id='confirm-dialog', component_property='displayed')],
    [Input('binList_tab2_id', 'value'),
     Input('memory_output_tab2_id', 'data')]
)

def update_datatable(binList,dataset) :  # https://kanoki.org/2019/04/06/pandas-map-dictionary-values-with-dataframe-columns/
    print('\n\nupdate scaffold datatable (TAB2)')
    print(binList)
    df_tmp = pd.DataFrame(dataset)

    if len(binList)==0:
        return dash.no_update,True
    else:
        return pd.merge( df_bin[ df_bin['bin'].isin(binList) ],df_tmp[ df_tmp['new_bin'].isin(binList) ] , on = 'scaffold' ).to_dict('records') , False


@app.callback(
    [Output('scatter_tab2_id', 'figure')],
    [Input('legend_dropdown_tab2_id', 'value'),
    Input('datatable_scaffold_tab2_id', 'data')]
)

def display_graph(legendValue,dataset) :
    print('\n\ndisplay_graph (TAB2)')
    print(legendValue)
    print(len(pd.DataFrame(dataset)))
    #print(dataset['type'])
    #print(dataset['namespace'].keys())
    fig = px.scatter(pd.DataFrame(dataset), x="refineM_gc", y="refineM_coverage", color=legendValue, facet_col="new_bin", facet_col_wrap=3,hover_name="scaffold", hover_data=["refineM_length", "class: taxa", "order: taxa", "family: taxa" , "genus: taxa" , "species: taxa" ])#,height=800)

    return [fig]


if __name__ == '__main__':
    app.run_server(host=ip_address,debug=True, port = int(port))
