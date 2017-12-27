#!/usr/bin/python
# -*- coding: utf-8 -*-

import psycopg2
import sys


con = None

try:
     
    con = psycopg2.connect(database='testdb', user='janbodnar') 
    cur = con.cursor()
    cur.execute('SELECT version()')          
    ver = cur.fetchone()
    print(ver)    
    


    
finally:
    
    if con:
        con.close()
