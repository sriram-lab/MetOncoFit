#!/usr/bin/env python
# -*- encoding: utf-8 -*-
"""
MetOncoFit URL Configuration

"""
from django.contrib import admin
from django.urls import path
from django.conf.urls import patterns, include, url
from django.views.generic import TemplateView

urlpatterns = patterns('',
    url(r'^$', TemplateView.as_view(template_name="index.html"))
)
