"""
"""

from flask import Flask
from flask.helpers import get_root_path

from config import BaseConfig

def create_app():
    server = Flask(__name__)
    server.config.from_object(BaseConfig)

    register_dashapps(server)
    register_extensions(server)
    register_blueprints(server)

    return server

def register_dashapps(app):
    from app.dashapp.layout import layout
    from app.dashapp.callbacks import register_callbacks

    # Metatags
    meta_viewport = {
        'name':'viewport',
        'content':'width=device-width, intiial-scale=1, shrink-to-fit=no'
    }

    dashapp = dash.Dash(__name__,
    server = app,
    url_base_pathname=
    )
