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



# -------------------------------------------------------------------------------------------------------------------------- #
# tab 1 #

sample = sys.argv[1]
directory = '/env/cns/proj/projet_CSD/scratch/assemblies/'+sample+'/refinedBins/output'
bins_summary_filename = '/env/cns/proj/projet_CSD/scratch/assemblies/'+sample+'/refinedBins/output/Bins_summary.tsv'

config_filename = '/env/cns/proj/projet_CSD/scratch/assemblies/'+sample+'/assembly/info.json'
with open(config_filename) as f:
    data = json.load(f)

project = data['project']
sample = data['sample']


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
    #print(gtdb_lineage)
    binName = suggestedName(anvio_lineage,gtdb_lineage,project,sample,binName2count)

    result,msg,color = checkBinName(binName,sample,project)
    if result :
        liste.append(binName)
    else:
        liste.append('TO BE DEFINED')
    data.append( liste )
file.close()


# Create the pandas DataFrame 
pd.options.display.float_format = '{:,.2f}'.format
df_summary = pd.DataFrame( data, columns = headerList )
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



# -------------------------------------------------------------------------------------------------------------------------- #
# tab 2 

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

options_checklist = []
for binName in binNameSet :
    options_checklist.append( {'label': binName, 'value': binName })

print('data_bin')
print(len(data_bin))
print(len(data_bin[0]))
print(data_bin[0])
print()

headerSublist = []
for i in [0,1,2,3,4,6,7,8,9,10,13,18,23,28,33,38,43]:
    headerSublist.append(headerList[i])

df_bin = pd.DataFrame( data_bin, columns = headerSublist )
#print(headerSublist)


#print(df_bin.dtypes)


# Create the pandas DataFrame 
#pd.options.display.float_format = '{:,.2f}'.format




#print(df.dtypes)

percentage = FormatTemplate.percentage(1)


columns_bin = [
    dict(id='scaffold', name='Scaffold'),
    {'id' : 'bin' , 'name' : 'Bin Name' , 'hideable' : False , 'presentation' : 'dropdown' , 'editable': True},
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
            dbc.Row(dbc.Col(html.H2(children='Bins summary'))),
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
                        size="50"
                        #bs_size="lg"
                        #placeholder='suggested name: XXX',
                        # style={
                        #     'font-size': "100%",
                        # },
                    ),width={'size' : '5' , 'order': '2'}
                ),
                dbc.Col(
                    html.Button(n_clicks=0, children='Submit', id='submit_new_bin_name')
                    ,width={'size' : '2' , 'order': '3'}
                ),
                dbc.Col(
                    html.Div(id='submit_new_bin_name_msg'),width={'size' : '3' , 'order' : '4'}
                )
            ]),
            dbc.Row(dbc.Col(html.Div(id='suggested_bin_name'))),
            dbc.Row(dbc.Col(html.Hr())),
            dbc.Row([
                dbc.Col(
                    html.Button(n_clicks=0, children='Check and Save', id='save_to_csv')
                    ,width={'size' : '2'}
                ),
                dbc.Col(
                    html.Div(id='save_to_csv_msg')
                    ,width={'size' : '4' , 'order' : '2'}
                )
            ])
        ]
    ),
    style={"width": "150%" , "className":"mt-3"},
)


bins_tab_content = html.Div(
    dbc.Card(
        dbc.CardBody( 
            [
                dbc.Row(dbc.Col(html.Hr())),
                dbc.Row(dbc.Col(html.H2(children='Bins distribution'))),
                dbc.Row(dbc.Col(html.Hr())),
                dbc.Row([
                    dbc.Col(
                        dcc.Checklist(
                            id='xaxis_column',
                            options=options_checklist,
                            value=list(binNameSet)
                        )#,width={'size' : '8' , 'order': '1' , 'offset' : '0'} # 'offset' : '0'  
                    )
                ]),

                dbc.Row([
                    dbc.Col(
                        dcc.Dropdown(
                            id='legend_dropdown',
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
                        ),width={'size' : '3' , 'order': '1' , 'offset' : '0'} # 'offset' : '0'  
                    ),
                    dbc.Col(
                        html.Button(n_clicks=0, children='Update scatterplot', id='update_scatter'),width={'order': '1' , 'size' : '3'}
                    )
                ]),
                dbc.Row(dbc.Col(dcc.Graph(id='scatter_id',style={'width': '100%', 'height': '75vh'}))),
                dbc.Row(dbc.Col(html.Hr())),
                dbc.Row(dbc.Col(
                    dash_table.DataTable(
                        css=[{"selector": ".Select-menu-outer", "rule": "display: block !important"}], # Add this line for dropdown
                        columns=columns_bin,
                        id = 'datatable_scaffold_id',
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
                            'bin': {
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
)


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
               Output('submit_new_bin_name_msg', 'children')],
              Input('submit_new_bin_name', 'n_clicks'),
              State('bin-dropdown', 'value'),
              State('input-1-state', 'value'))


def update_datatable(n_clicks, selectedBin, binName):
    if selectedBin != None :
        isBin = df_summary["Bin"] == selectedBin
        row_index = df_summary.index[isBin].tolist()[0]

        print(isBin)
        print(row_index)
        print(df_summary.head())
        print('input: ')
        print('input 1: '+str(selectedBin))
        print('input 2: '+str(binName))

        result,msg,colorAlert = checkBinName(binName,sample,project) # have to check for redundancy

        if result :
            # check for uniqness
            result2,msg2,colorAlert2 = checkUniqness(df_summary,row_index,binName)
            if result2 :
                df_summary.at[row_index,'Name'] = binName
                print(df_summary.head())
            alertButton = dbc.Alert(msg2, color=colorAlert2,duration=3000)
        else:
            alertButton = dbc.Alert(msg, color=colorAlert,duration=3000)

    else:
        alertButton = ''


    return df_summary.to_dict('records'),alertButton


@app.callback(
    [Output('save_to_csv_msg', 'children')],
    [Input('save_to_csv', 'n_clicks')],
    [State('datatable_bin_summary', 'data')]
)

def df_to_csv(n_clicks, dataset ):    
    input_triggered = dash.callback_context.triggered[0]["prop_id"].split(".")[0]
    print(input_triggered)

    if input_triggered == "save_to_csv":
        df = pd.DataFrame(dataset)
        df.to_csv("Your_Bin_Summary_Data.csv")

        return [dbc.Alert("The data has been saved to your folder.", color="success",duration=3000)]
    else:
        return ['']

    
##################
# callbacks tab2 #
##################


# https://dash.plotly.com/datatable/interactivity
# https://community.plotly.com/t/dash-table-datatable-filtering-and-sorting-doesnt-seem-to-work/16362

@app.callback(
    [Output('datatable_scaffold_id', 'data')],
    [Input('xaxis_column', 'value'),
     Input('update_scatter','n_clicks')],
    [State('datatable_scaffold_id', 'data')]
)

def display_datatable(checkboxList,n_clicks,dataset) :
    print('update datatable')
    print(checkboxList)
    print(columns_bin)
    print(len(df_bin[ df_bin['bin'].isin(checkboxList) ]))

    input_triggered = dash.callback_context.triggered[0]["prop_id"].split(".")[0]
    if input_triggered == "update_scatter": # updating df_bin
        for row in dataset :
            isScaffold = df_bin["scaffold"] == row['scaffold']
            row_index = df_bin.index[isScaffold].tolist()[0]
            if df_bin.at[row_index,'bin'] != row['bin'] :
                print(row['scaffold']+'\t'+df_bin.at[row_index,'bin']+'\t==>\t'+row['bin'])
                df_bin.at[row_index,'bin'] = row['bin']

    return [df_bin[ df_bin['bin'].isin(checkboxList) ].to_dict('records')]


@app.callback(
    [Output('scatter_id', 'figure')],
    [Input('legend_dropdown', 'value'),
    Input('datatable_scaffold_id', 'data')]
)

def display_graph(legendValue,dataset) :
    print('display_graph')
    print(legendValue)
    print(len(pd.DataFrame(dataset)))
    #print(dataset['type'])
    #print(dataset['namespace'].keys())
    fig = px.scatter(pd.DataFrame(dataset), x="refineM_gc", y="refineM_coverage", color=legendValue, facet_col="bin", facet_col_wrap=3,hover_name="scaffold", hover_data=["refineM_length", "class: taxa", "order: taxa", "family: taxa" , "genus: taxa" , "species: taxa" ])#,height=800)

    return [fig]

    
if __name__ == '__main__':
    print('test')
    app.run_server(host=ip_address,debug=True, port = int(port))
