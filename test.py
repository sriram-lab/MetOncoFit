# -*- encoding: utf-8 -*-
import flask
metoncofit = flask.Flask(__name__)

@metoncofit.route('/')

def hello():
    return "hello"

if __name__ == "__main__":
    metoncofit.run()
