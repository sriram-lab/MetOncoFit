#! /usr/bin/python3.6

import logging
import sys
from run_app import app as application

logging.basicConfig(stream=sys.stderr)
sys.path.insert(0, '/home/metoncofit-website/metoncofit/')

application.secret_key = "ncrc-520-1225"
