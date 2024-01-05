import psycopg2  # PostgreSQL
#import pubchempy as pcp  # PubChem baza podataka
#import csv

pgUser = "postgres"
pgPassword = "root"



def SearchDescriptor(descriptor):
    
    if ", " in descriptor:
        descriptor = descriptor.strip().split(", ")
    elif "," in descriptor:
        descriptor = descriptor.strip().split(",")
    elif " " in descriptor:
        descriptor = descriptor.strip().split(" ")
    
    
    
    conn = psycopg2.connect(
        database="kemoinformatika",
        user=pgUser,
        password=pgPassword,
        host='127.0.0.1',
        port='5432'
        )
    cursor = conn.cursor()

    if isinstance(descriptor, str):
        sql = """select * 
        from molekule as m
        inner join svojstva as s
        on m.id_molekule = s.id_molekule
        where s.mirisno_svojstvo like '{desc}'""".format(desc = descriptor)
        cursor.execute(sql)
        x = cursor.fetchall()
    
    else:
        sql = """select * 
        from molekule as m
        inner join svojstva as s
        on m.id_molekule = s.id_molekule
        where """
        for i in descriptor:
            sulfix = "s.mirisno_svojstvo like '{desc}' or ".format(desc = i)
            sql = sql + sulfix
        cursor.execute(sql[:-3])
        x = cursor.fetchall()

    conn.close()
    return x



def SearchSmiles(smiles):
    
    conn = psycopg2.connect(
        database="kemoinformatika",
        user=pgUser,
        password=pgPassword,
        host='127.0.0.1',
        port='5432'
        )
    cursor = conn.cursor()
    sql = """ select * 
    from molekule
    where smiles_kod like '{sm}'
    """.format(sm = smiles)
    cursor.execute(sql)
    x = cursor.fetchone()
    sql = """ select mirisno_svojstvo
    from svojstva
    where id_molekule = {idm}
    """.format(idm = x[0])
    cursor.execute(sql)
    y = cursor.fetchall()
    return x,y
    

def SearchName(name):
        
    conn = psycopg2.connect(
        database="kemoinformatika",
        user=pgUser,
        password=pgPassword,
        host='127.0.0.1',
        port='5432'
        )
    cursor = conn.cursor()
    sql = """ select * 
    from molekule
    where ime_molekule like '{nm}'
    """.format(nm = name)
    cursor.execute(sql)
    x = cursor.fetchone()
    sql = """ select mirisno_svojstvo
    from svojstva
    where id_molekule = {idm}
    """.format(idm = x[0])
    cursor.execute(sql)
    y = cursor.fetchall()
    return x,y
    

#print(SearchDescriptor("fishy nutty"))
#print(SearchSmiles("CC(=O)C1NCCS1"))
#print(SearchName("1-(1,3-thiazolidin-2-yl)ethanone"))