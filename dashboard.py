#! /usr/bin/env python

# -*- coding: utf-8 -*-                                                                                                                                                                                                                                                      


import os,sys,re
import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.express as px
import pandas as pd
from dash.dependencies import Input, Output


## importing socket module
import socket
## getting the hostname by socket.gethostname() method
hostname = socket.gethostname()
## getting the IP address using socket.gethostbyname() method
ip_address = socket.gethostbyname(hostname)
## printing the hostname and ip_address
print(f"Hostname: {hostname}")
print(f"IP Address: {ip_address}")






binFilenameSet = set()
filenameSet = set( ['Collection.xlsx' , 'Anvio_summary.tsv' , 'CheckM.tsv' , 'GTDBtk.tsv' , 'Collection.tsv' ] )
directory = '/env/cns/proj/projet_CSD/scratch/assemblies/Efas_M_CS1/refinedBins/output'
for root, dirs, files in os.walk(directory):
    for filename in files :
        if filename in filenameSet :
            continue
        else:
            binFilenameSet.add(root+'/'+filename)
            print(root+'/'+filename)


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
df = pd.read_csv('/env/cns/proj/agc/home/rmeheust/scripts/bins.info',sep="\t")

fig = px.scatter(df, x="bin", y="Mean coverage", color="anvio_taxonomy")



app = dash.Dash(__name__)

app.layout = html.Div([

    html.H1(children='Hello Dash'),
    html.Div(children='Dash: A web application framework for Python.'),
    dcc.Graph(id='boxplot-coverage',figure=fig),

    html.Div([
        dcc.Checklist(
            id='xaxis-column',
            options=options,
            value=list(binNameSet)
        )
    ]),
    html.Div(id='my-output'),

])

@app.callback(
    Output('boxplot-coverage', 'figure'),
#    Output(component_id='my-output', component_property='children'),
    Input('xaxis-column', 'value')
)

def update_graph(checkboxList) :
    print(checkboxList)
    dff = df[ df['bin'].isin(checkboxList) ]
    print('dff: '+str(set(dff['bin'])))
    fig = px.scatter(dff, x="bin", y="Mean coverage", color="anvio_taxonomy")
    return fig




if __name__ == '__main__':
    print('test')

    app.run_server(host=ip_address,debug=True, port = 8080)
