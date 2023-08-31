# This package contains functions to work with sqlite3 database.

from projdirs import datadir
import sqlite3
import pandas as pd
import numpy as np
import sys,os

if sys.platform == 'win32':
    connector = '\\'
elif sys.platform == 'linux':
    connector = '/'
else: connector = 'thequickbrownfoxjumpsoverthelazydog'

def list_tables(db, db_dir = datadir + 'arbitrage%s' %connector):
    """ Finds the list if tables in a sql database file.
    The name of the database file is set by db.
    db is the name of the database file with extension '.db' """

    connection = sqlite3.connect(db_dir + db)
    cursor = connection.cursor()
    cursor.execute("SELECT name from sqlite_master where type='table'")
    List=([row for row in cursor.fetchall()[0]])
    connection.close
    return(List) 


def list_columns(db, table, data_type=True, 
                 db_dir = datadir + 'arbitrage%s' %connector):
    """This function lists the column names of a table
    in a database file.
    db: database file 'file_name.db'
    table: 'table_name'
    datatype: whether you want to see the type of data in
    each column ot not.
    datatype=False returns a list of column names"""
    data_file = db_dir + db
    if not os.path.isfile(data_file):
          print('Data file not found!')
          return()
    else:
          connection = sqlite3.connect(data_file)
          cursor = connection.cursor()
          cursor.execute('PRAGMA table_info(' + table+ ')')
          if data_type:
                List = [row[1:3] for row in cursor.fetchall()]
          else:
                List = [row[1] for row in cursor.fetchall()]
          connection.close
          return(List)

def list_unique_index(db, db_dir = datadir + 'arbitrage%s' %connector):
    """This function lists all the unique indices of a database.
    db: database file e.g. 'file_name.db' """

    connection = sqlite3.connect(db_dir + db)
    cursor = connection.cursor()
    cursor.execute(
        "select type, name, tbl_name, sql from sqlite_master where type='index'")
    List=([row for row in cursor.fetchall()])
    connection.close
    return(List)

def get_table(db, table, db_dir = datadir + 'arbitrage%s' %connector):
    """This function loads all columns from a table in a database file.
    table: database table name 'table_name'
    db: database name ('file_name.db')
    
    THE PV PLANT CLASS DEPENDS ON THIS FUNCTION. IF CHANGES ARE MADE THEN
    THE PV PLANT CLASS SHOULD BE TESTED FOR SUCCESSFUL INSTANTIATION.
    
    Note from Jeff: I don't know how to work with data types when importing 
    SQL data so I've done data conversions after importing data from the 
    SQL database. There is probably a better way to do this based on the SQL
    data types that I don't know about.

    Example:
    cols = ['state','year','sh','obj','rte', 'loss', 'cap']
    table = 'storage_value'
    db = 'storage_value_arbitrage.db'"""
    
    data_file = db_dir + db
    if not os.path.isfile(data_file):
          print('Data file not found!')
          return()
    else:
        columns = np.array(list_columns(db, table, db_dir = db_dir, data_type = True))
        connection = sqlite3.connect(data_file)
        cursor = connection.cursor()
        cursor.execute("SELECT * FROM " + table)
        data = np.array(cursor.fetchall(), dtype = object)
        data[:,columns[:,1]=='REAL'] = data[:,columns[:,1]=='REAL'].astype(float)
        data[:,columns[:,1]=='INTEGER'] = data[:,columns[:,1]=='INTEGER'].astype(int)
        df = pd.DataFrame(data[:,1:], index = data[:,0],
                          columns = columns[1:,0])
        connection.close
        return(df)

def get_data(cols, table, db, db_dir = datadir + 'arbitrage%s' %connector):
    """This function loads the specified columns from a table in a database file.
    cols: a list of requested database column names [col1, col2, ...]
    table: database table name 'table_name'
    db: database name ('file_name.db')

    Example:
    cols = ['state','year','sh','obj','rte', 'loss', 'cap'] if not known enter 'N/A'
    table = 'storage_value'; if not known enter 'N/A'
    db = 'storage_value_arbitrage.db'
    sm.get_data(cols, table, db)"""
    data_file = db_dir + db
    if table=='N/A':
          table = list_tables(db)[0]
    if cols=='N/A':
          cols = [col[0] for col in list_columns(db, table)]
    if not os.path.isfile(data_file):
          print('Data file not found!')
          return()
    else:
          connection = sqlite3.connect(data_file)
          cursor = connection.cursor()
          cursor.execute("SELECT " + "["+"], [".join(cols)+"]" + " FROM " + table)
          data = pd.DataFrame(cursor.fetchall())
          connection.close
          data.columns = cols
          return(data)

def create_table(table, db, cols, create_unique_idx=False, idx_cols=[],
                 db_dir = datadir + 'arbitrage%s' %connector):
    """ This function creates a new table in a new/existing database file.
    table: 'table_name'
    db: database file name 'file_name.db'
    cols: can contain the data format e.g. cols = ['state CHAR(3)', 'year INT']
    'create_unique_idx': whether to create a unique index on the table or not
    'unique_idx': what columns to use for the unique index; 
    if create_unique_idx==False, then no entry is required for unique_idx
    """
    connection = sqlite3.connect(db_dir + db)
    cursor = connection.cursor()
    cursor.execute(
        "SELECT name FROM sqlite_master WHERE type='table' AND name="+
        str(['table_test']).strip('[]'));
    sql_command = "CREATE TABLE IF NOT EXISTS %s (%s)" %(table, ",".join(cols))
    cursor.execute(sql_command)
    connection.commit()
    connection.close
    if create_unique_idx:
        create_unique_index(db,table, idx_cols)
    return('Database table was created!')
    

def create_unique_index(db,table, cols, idx='idx',
                        db_dir = datadir + 'arbitrage%s' %connector):
    """This function creates a unique index in a table if the unique index doesn't
    already exist.
    table: 'table_name'
    db: database file name '*.db'
    cols: a list of column names in the db e.g. ['col1', 'col2', ...]
    idx: index name, 'ID' """
    
    connection = sqlite3.connect(db_dir + db)
    cursor = connection.cursor()
    sql_command = "CREATE UNIQUE INDEX IF NOT EXISTS %s ON %s (%s)"%(idx, table, ",".join(cols))
    cursor.execute(sql_command)
    connection.commit()
    connection.close
    return('Unique index was created or already existed!')

def replace_into_db(df, db, table, cols, create_unique_idx=False, idx_cols=[],
                    db_dir = datadir + 'arbitrage%s' %connector):
    """
    This function writes data from 
    cols: a list of column names in the db e.g. ['col1', 'col2', ...]
    """
    
    sqlite3.register_adapter(np.int64, lambda val: int(val))
    sqlite3.register_adapter(np.int32, lambda val: int(val))
    create_table(table, db, cols, create_unique_idx, idx_cols)
    connection = sqlite3.connect(db_dir + db)
    cursor = connection.cursor()
    sql_command = "REPLACE INTO %s (%s) VALUES (%s)"%(table, ",".join(cols),
                                ",".join(np.repeat('?', len(cols)).tolist()))
    res = df.to_records(index=False)
    cursor.executemany(sql_command, res)
    connection.commit()
    connection.close
    return('Data was recorded into the database!')

def add_column_into_db(db,table,df,update_uidx,idx,
                       db_dir = datadir + 'arbitrage%s' %connector):
    """
    db: the name of the database to be edited
    table: the name of the table in the database to be edited
    df: the dataframe data that replaces the existing table.
    update_uidx: True or False to create a unique index on the new table
    idx: the list of column names for the unique index
    
    """
    sqlite3.register_adapter(np.int64, lambda val: int(val))
    sqlite3.register_adapter(np.int32, lambda val: int(val))
    connection = sqlite3.connect(db_dir + db)
    cursor = connection.cursor()
    cursor.execute("DROP TABLE "+table)
    connection.commit()
    connection.close
    replace_into_db(df, db, table, df.columns.tolist(),
                    create_unique_idx=update_uidx,
                    idx_cols=idx)


def delete_table(db,table, db_dir = datadir + 'arbitrage%s' %connector):
    """
    db: the name of the database to be edited
    table: the name of the table in the database to be edited
    df: the dataframe data that replaces the existing table.
    update_uidx: True or False to create a unique index on the new table
    idx: the list of column names for the unique index
    
    """
    sqlite3.register_adapter(np.int64, lambda val: int(val))
    sqlite3.register_adapter(np.int32, lambda val: int(val))
    connection = sqlite3.connect(db_dir + db)
    cursor = connection.cursor()
    cursor.execute("DROP TABLE "+table)
    connection.commit()
    connection.close

