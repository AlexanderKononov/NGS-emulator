import psycopg2
import sys

con = None

try:

    con = psycopg2.connect(database='testdb', user='janbodnar')
    cur = con.cursor()
    cur.execute('SELECT version()')
    ver = cur.fetchone()
    print(ver)


#except psycopg2.DatabaseError:
 #   print('Error')
    # sys.exit(1)


finally:

    if con:
        con.close()