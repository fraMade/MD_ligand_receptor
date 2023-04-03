from plotly.subplots import make_subplots
import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
import numpy as np
import json



# get subplot
def get_go_express(fig_ex):
    fig_traces = []
    for trace in range(len(fig_ex["data"])):
        fig_traces.append(fig_ex["data"][trace])
    return fig_traces[0]


'''
Given the permanences data as a dataframe, the function return a plot to visualize the informations.
If no data is available it will return a "No recorded interactions" caption.
'''
def permanence_to_bar(permanenceData):
    
    ## check if the dataframe is empty
    if permanenceData.empty:
        fig = go.Figure()
        fig.update_layout(
            xaxis =  { "visible": False },
            yaxis = { "visible": False },
            annotations = [
                {   
                    "text": "No recorded interactions ",
                    "xref": "paper",
                    "yref": "paper",
                    "showarrow": False,
                    "font": {
                        "size": 28
                    }
                }
            ]
        )
        return fig ## return a plot with "No recorded interactions " caption
    
    return  px.bar(permanenceData, x='BondType', y='#Timestamps', color = 'BondType', title = 'Ligand-Receptor interactions per bond type')  


'''
Given a dataframe containing the interactions done for each atom, the function return a stacked-bar plot showing the count of interactions for each atom, divided by bond-type. The function also filters out the atoms that are not active for at least a specified threshold percentage.
If no data is available it returns a plot with "No recorded interactions" caption.
'''
def atoms_permanence_to_barplot(atoms_permanences, threshold):
    
    if  atoms_permanences.empty:
        fig = go.Figure()
        fig.update_layout(
            xaxis =  { "visible": False },
            yaxis = { "visible": False },
            annotations = [
                {   
                    "text": "No recorded interactions",
                    "xref": "paper",
                    "yref": "paper",
                    "showarrow": False,
                    "font": {
                        "size": 28
                    }
                }
            ]
        )
        return fig
    
    temp = atoms_permanences[atoms_permanences["percentage"] >= threshold].iloc[:,:-1]
    temp.index.name = 'Ligand atoms'
    temp = pd.melt(temp.reset_index(), id_vars='Ligand atoms', value_vars=['hbonds',
                                                                           'pi-stacks',
                                                                           'hydrophobic-interactions',
                                                                           'water-bridge'], 
                   var_name='Bond types', value_name='#timeframes')
    
    fig = px.bar(temp, x='Ligand atoms', y='#timeframes', barmode='stack', color = 'Bond types', title = 'atoms permanence with ' + str(threshold)+"% threshold")
    fig.update_layout(
    xaxis = go.layout.XAxis(
        tickangle = 90)
    )
        
        
    fig.update_layout(
    xaxis = dict(
        tickmode = 'array',
        tickvals = temp['Ligand atoms'],
        ticktext = temp['Ligand atoms']
    )
)

    return fig




'''
The function return a multiplot composed by a scatter plot, to visualize the timeframes in which the atoms are involved, and a bar plot to visualize the permanence percentage.
'''
def plt_res_activity(df, timestamp_num, what='', threshold = 0):
    if not df.empty:
        prcs = df.groupby(['Receptor']).pipe(lambda grp: round(grp.count() / (timestamp_num+1) *100, 2)).reset_index().rename({"Timestamp": 'Percentage'}, axis=1).sort_values('Receptor', ascending = False)
        filteredprcs= prcs[prcs['Percentage']>threshold]
        df = df[df['Receptor'].isin(filteredprcs['Receptor'])]

    if  df.empty  :
        fig = go.Figure()
        fig.update_layout(
            xaxis =  { "visible": False },
            yaxis = { "visible": False },
            annotations = [
                {   
                    "text": "No recorded interactions for "+what,
                    "xref": "paper",
                    "yref": "paper",
                    "showarrow": False,
                    "font": {
                        "size": 28
                    }
                }
            ]
        )
        return fig
    fig = make_subplots(
    rows=1, cols=2, column_widths=[0.6, 0.4],subplot_titles=(what,'Permanence'), shared_yaxes=True, horizontal_spacing = 0.05)
    
    fig1 = go.Figure(data=go.Scattergl(
    x = df["Timestamp"],
    y = df["Receptor"],
    mode='markers',
    marker_size=10,
    marker_symbol="line-ns"
    ))
    # fig1.update_traces(marker_symbol="line-ns", selector=dict(type='scattergl'))
    # fig1 = px.scatter(df, y="Receptor", x="Timestamp", symbol_sequence = ['line-ns'])
    fig.add_trace(get_go_express(fig1), row=1, col=1)
    
    
    fig2 = px.bar(filteredprcs, y="Receptor", x="Percentage", text_auto=True)
  
    fig.add_trace(get_go_express(fig2), row=1, col=2)
    fig['layout']['title'].update(text= what +  " interactions")
    fig['layout']['xaxis2'].update(range=[0,100])
    fig2.update_xaxes(range=[0,100] , dtick=5, autorange=False)

    
    return fig
    
    
    
'''
The function return an heatmap plot showing the percentage of interactions between any pair of atoms.
'''
def pltres_lig_heatmap(df, prf, timestamp_number, min_val=0):
    
    df = df.applymap(lambda x: round((float(x/(timestamp_number+1))*100)))
    df = df[(df>=min_val).any(axis=1)]
    df = df.T
    df = df[(df>=min_val).any(axis=1)]
    if df.empty:
        fig = go.Figure()
        fig.update_layout(
            xaxis =  { "visible": False },
            yaxis = { "visible": False },
            annotations = [
                {   
                    "text": "No recorded interactions for selected bond type ",
                    "xref": "paper",
                    "yref": "paper",
                    "showarrow": False,
                    "font": {
                        "size": 28
                    }
                }
            ]
        )
        return fig
   
    zz = df.to_numpy().tolist()
    zz_test = [list(map(str, zzz)) for zzz in zz]
    fig = px.imshow(df, aspect="auto", title = "Ligand-"+prf+" interactions",range_color =[0,100])
    fig.update_traces(text = zz_test, texttemplate = "%{text}")
    fig.update_layout(coloraxis_colorbar=dict(
    title="Permanence",
    ticks="outside", ticksuffix=" %",
))
    return fig


                        
                
                        