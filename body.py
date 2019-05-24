"""
body for the website
@author: Scott Campit
"""

# Create content and store as "cards"
_body = dbc.Container(
    dbc.Row(
        [
            dbc.Col([
                html.Div([
                    html.H1("MetOncoFit")
                    ], style={'marginTop':30, 'marginBottom':30}),
                html.Div([
                    html.H4("A Machine Learning Algorithm to Identify Common Biochemical and Topological Attributes of Metabolic Genes Recurrently Dysregulated in Tumors"),
                    html.P(
                    """
                    MetOncoFit is a data-driven approach that uses biochemical and metabolic attributes to predict tumor differential expression, copy number variation, and patient survival.
                    """
                    )
                ], style={'marginTop':30, 'marginBottom':30}),
                html.Div([
                    html.H4("Installing MetOncoFit"),
                    html.P(
                    """
                    We recommend installing MetOncoFit inside a virtual environment.

                    To install the package, use PyPI:
                    """
                    )
                ], style={'marginTop':30, 'marginBottom':30}),
                html.Div([
                    html.P('''
                    > pip install metoncofit
                    ''')],
                    style={'marginTop':10, 'marginBottom':10, 'font-family':'courier', 'background-color':'#393939', 'color':'#EC7F37', 'text-indent':'1em'}
                ),
                html.Div([
                    html.H4("Contributing to the MetOncoFit project"),
                    html.P(
                    """
                    Contributions are welcome! Please read the contributions guide to get started.
                    """),
                    html.P(
                    """
                    """
                    ),
                    html.P(
                    """
                    To support the MetOncoFit project, you can cite our publication:
                    """
                    ),
                    html.P(
                    """
                    """
                    ),
                    html.P(
                    """
                    Oruganty, K., Campit, S.E., Mamde, S., & Chandrasekaran, S. Common biochemical and topological attributes of metabolic genes recurrently dysregulated in tumors.
                    """
                    )
                ], style={'marginTop':30, 'marginBottom':30}),
                html.Div([
                    html.H4("MetOncoFit interactive explorer"),
                ])
            ])
        ]
    )
)
