from dash import dcc, html
layout = html.Div([
    #first row
    html.Div(children=[
         #first column
        html.Div([
        dcc.Upload(
            id='upload-data-json',
            children=html.Div([
                'Drag and Drop JSON FILE or ',
                html.A('Select Files')
                ]),
            style={
            'width': '80%',
            'height': '60px',
            'lineHeight': '60px',
            'borderWidth': '1px',
            'borderStyle': 'dashed',
            'borderRadius': '5px',
            'textAlign': 'center',
            'margin': 'auto'
        },
        ),
                ]),
        html.Div([
         dcc.Upload(
             
            id='upload-data-csv',
            children=html.Div([
                'Drag and Drop LIG CSV FILE or ',
                html.A('Select Files')
                ]),
            style={
            'width': '80%',
            'height': '60px',
            'lineHeight': '60px',
            'borderWidth': '1px',
            'borderStyle': 'dashed',
            'borderRadius': '5px',
            'textAlign': 'center',
            'margin': 'auto'
        },
        ),
                ]),
        html.Div([   
         dcc.Upload(
             
            id='upload-rcpt-csv',
            children=html.Div([
                'Drag and Drop RCPT CSV FILE or ',
                html.A('Select Files')
                ]),
            style={
            'width': '80%',
            'height': '60px',
            'lineHeight': '60px',
            'borderWidth': '1px',
            'borderStyle': 'dashed',
            'borderRadius': '5px',
            'textAlign': 'center',
            'margin': 'auto'
        },
        ),
             
                ]),    
            ],style={'display': 'inline-block','width':'100%', 'margin': 'auto'}),


        dcc.Store(id ='output-csv-data-upload'),
        dcc.Store(id ='output-json-data-upload'),
        dcc.Store(id = 'output-csv-rcpt-upload'),
        
       
        
         html.Div(children=[
             html.P( "Filter by treshold (or select range in histogram):", className = "control_label"),
						dcc.Slider(id = "ts_slider",
							            min = 0,
							            max = 100,
                                        
							            value = 1,
							            className = "dcc_control"),
         ],style={"margin-top":"2%"}), #style={'display': 'inline-block', 'vertical-align': 'top', 'margin-left': '3vw', 'margin-top': '3vw'}),
            
         
         html.Div([
				html.Div(   [dcc.Graph(id = "bondp")],
                            id = "bondpContainer",
                            className = "pretty_container six columns",
				),
                html.Div(   [dcc.Graph(id = "prembond")],
                            id = "prembondContainer",
                            className = "pretty_container six columns",
				),
			],
			#className = "row flex-display",
		),
            
    
         html.Div(children=[
             html.Button("Download CSV", id="bondp_csv",style={'margin-left':'5%'}),
             dcc.Download(id="download-bondp-csv"),
             
             html.Button("Download CSV", id="prembond_csv",style={'margin-left':'41%'}),
             dcc.Download(id="download-prembond-csv"),
             
            ],className = "row flex-display",),
    
          #second row
         
         html.Div(children=[
             #second column
             html.Div(children=[
                html.P("Filter by Bond type:", className = "control_label"),
                            dcc.Dropdown(   id = "bondtype",
                                            options = [ {"label": "Hbonds", "value": "hbonds"},
                                                        {"label": "Pi-stacks", "value": "pi-stacks"},
                                                        {"label":"Hydrophobic", "value": "hydrophobic-interactions"},
                                                        {"label":"Water-bridge", "value": "water-bridge"}
                                                         ],
                                            multi = True,
                                            value = ["hbonds","pi-stacks","hydrophobic-interactions","water-bridge"],
                                            className = "dcc_control"),
          html.Div([
				html.Div(   [dcc.Graph(id = "nucleomap")],
                            id = "nucelomapContainer",
                            className = "pretty_container six columns",
				),
                html.Div(   [dcc.Graph(id = "aminomap")],
                            id = "aminomapContainer",
                            className = "pretty_container six columns",
				),
			],
			className = "row flex-display",
		),
                 
             ]),
          ]),
            dcc.Store(id='nuc_csv'),
            dcc.Store(id='ammin_csv'),
            dcc.Store(id= "atom_res_data"),
            dcc.Store(id= "atom_data"),
            dcc.Store(id="dna-lig-activity"),
            dcc.Store(id="prot-lig-activity"),
            dcc.Store(id="mapatoms"),
            dcc.Store(id="timestamp_num"),
            dcc.Store(id="atoms_permanences_per_bondtype"),
            dcc.Store(id="ligands_intercations_per_bondtype"),
   
    
        html.Div(children=[
             html.Button("Download CSV", id="nucleomap_csv",style={'margin-left':'5%'}),
             dcc.Download(id="download-nucleomap-csv"),
            
              html.Button("Download CSV", id="aminomap_csv",style={'margin-left':'41%'}),
             dcc.Download(id="download-aminomap-csv"),
             
            ], style= {}),
    
        

           
            #third row 
            html.Div(children=[
    
             html.P( "Filter interaction activity by treshold :", className = "control_label"),
						dcc.Slider(id = "inter_slider",
							            min = 0,
							            max = 100,
                                        
							            value = 1,
							            className = "dcc_control"),
					],),
    
         html.Div([
				html.Div(   [dcc.Graph(id = "resactivity")],
                            id = "resactivityContainer",
                            className = "pretty_container twelve columns",
				),
			],
			className = "row flex-display",
		),
    
        html.Div(children=[
             html.Button("Download CSV", id="resactivity_csv",style={'margin-left':'5%'}),
             dcc.Download(id="download-resactivity-csv"),
        ]),
    
    
       html.Div([
				html.Div(   [dcc.Graph(id = "resactivityamino")],
                            id = "resactivityaminoContainer",
                            className = "pretty_container twelve columns",
				),
			],
			className = "row flex-display",
		),
    
     html.Div(children=[
             html.Button("Download CSV", id="resactivityamino_csv",style={'margin-left':'5%'}),
             dcc.Download(id="download-resactivityamino-csv"),
        ]),
    
     #second row
         html.Div(children=[
    
                            html.P("Select the atom:", className = "control_label"),
                            dcc.Dropdown(   id = "atom",
                                            options = [''],
                                            #multi = False,
                                            #value = atoms[0],
                                            className = "dcc_control"),
                            dcc.Dropdown(   id = "recs",
                                            options =  [''],
                                            #multi = False,
                                            #value = recs[0],
                                            className = "dcc_control"),

             ],),#style={'display': 'inline-block', 'vertical-align': 'top', 'margin-left': '3vw', 'margin-top': '3vw'}),

        
     html.Div([
				html.Div(   [dcc.Graph(id = "atompermbondtype")],
                            id = "atompermbondtypeContainer",
                            className = "pretty_container twelve columns",
				),
			],
			className = "row flex-display",
		),
    
      html.Div(children=[
             html.Button("Download CSV", id="atom_csv",style={'margin-left':'5%'}),
             dcc.Download(id="download-atom_activity-csv"),
        ]),
    
        html.Div([
				html.Div(   [dcc.Graph(id = "atompermbondtype2")],
                            id = "atompermbondtype2Container",
                            className = "pretty_container twelve columns",
				),
			],
			className = "row flex-display",
		),
        html.Div(children=[
             html.Button("Download CSV", id="atom_res_csv",style={'margin-left':'5%'}),
             dcc.Download(id="download-atom_res_activity-csv"),
        ]),
        
    ])