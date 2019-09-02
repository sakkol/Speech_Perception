# -*- coding: utf-8 -*-
"""
Created on Mon Sep  2 13:33:16 2019

@author: newuser
"""

def explicit():
    from google.cloud import storage

    # Explicitly use service account credentials by specifying the private key
    # file.
    storage_client = storage.Client.from_service_account_json('My First Project-dc1f01a6b01a.json')

    # Make an authenticated API request
    buckets = list(storage_client.list_buckets())
    print(buckets)