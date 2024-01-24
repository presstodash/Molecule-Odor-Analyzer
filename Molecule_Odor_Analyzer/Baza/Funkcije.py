import psycopg2  # PostgreSQL
import pubchempy as pcp  # PubChem baza podataka
#import csv
from rdkit import Chem
from rdkit.Chem import Draw
from Molecule_Odor_Analyzer.Baza.Baza_account import getUser, getPassword
#from Baza_account import getPassword, getUser

pgUser = getUser()
pgPassword = getPassword()



def get_compound_details(smiles_code):
    compound = pcp.get_compounds(smiles_code, 'smiles')[0]  # Fetch the first compound from the search results
    compound_id = compound.cid  # Compound ID
    name = compound.iupac_name  # Compound name
    formula = compound.molecular_formula  # Compound formula
    weight = compound.molecular_weight  # Molecular weight
    smiles = compound.canonical_smiles
    return [compound_id, name, smiles, formula, weight]


def get_compound_details_name(name):
    compound = pcp.get_compounds(name, 'name')[0]  # Fetch the first compound from the search results
    compound_id = compound.cid  # Compound ID
    name = compound.iupac_name  # Compound name
    formula = compound.molecular_formula  # Compound formula
    weight = compound.molecular_weight  # Molecular weight
    smiles = compound.canonical_smiles
    return [compound_id, name, smiles, formula, weight]


def SearchDescriptor(descriptor):
    
    if ", " in descriptor:
        descriptor = descriptor.strip().split(", ")
    elif "," in descriptor:
        descriptor = descriptor.strip().split(",")
    elif " " in descriptor:
        descriptor = descriptor.strip().split(" ")
    
    # tu za svaki element treba spellchecker
    
    conn = psycopg2.connect(
        database="kemoinformatika",
        user=pgUser,
        password=pgPassword,
        host='127.0.0.1',
        port='5432'
        )
    cursor = conn.cursor()

    if isinstance(descriptor, str):
        sql = """select m.* 
        from molekule as m
        inner join svojstva as s
        on m.id_molekule = s.id_molekule
        where s.mirisno_svojstvo like '{desc}'""".format(desc = descriptor)
        cursor.execute(sql)
        x = cursor.fetchall()
    
    else:
        sql = """SELECT m.*
        FROM Molekule m
        JOIN (
            SELECT ID_Molekule
            FROM Svojstva
            WHERE Mirisno_Svojstvo IN ('{desc}'""".format(desc = descriptor[0])
        
        for i in range(1, len(descriptor)):
            sulfix = ", '{desc}'".format(desc = descriptor[i])
            sql = sql + sulfix
        sql = sql +""")  
                GROUP BY ID_Molekule
                HAVING COUNT(DISTINCT Mirisno_Svojstvo) >= {num}
            ) s ON m.ID_Molekule = s.ID_Molekule; """.format(num = len(descriptor))
        cursor.execute(sql)
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
    sm = 1
    x = cursor.fetchone() # ako postoji molekula, netreba se spajat na pubchem
    #print(x)
    if x != None:
        sql = """ select mirisno_svojstvo
        from svojstva
        where id_molekule = {idm}
        """.format(idm = x[0])
        cursor.execute(sql)
        y = cursor.fetchall()
    else:
        try:
            x = get_compound_details(smiles) # testirati koji error dolazi kad se ne da smiles kod
            # y je popis deskriptora koje AI daje
            y = [] # za sad
        except Exception as e:
            #print("Online fail")
            sm = None
            y=[]
        
    if sm is not None:
        mol = Chem.MolFromSmiles(x[2])# provjerit dal u ovom trenutku postoji smiles kod, i nastavit
        if mol is not None:
            img = Draw.MolToImage(mol, size=(300, 300))
            img = img.convert('RGBA')

            imgdata = img.getdata()
            newimgdata = [(r, g, b, 0) if ((r > 100) and (g > 100) and (b > 100)) else (r, g, b, a) for r, g, b, a in imgdata]
            img.putdata(newimgdata)
            
    else:
        img = None
    
    return x,y, img
    

def SearchName(name):
    
    name = name.lower()
    
    # tu ubaci spellchecker
        
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
    sm = 1
    cursor.execute(sql)
    x = cursor.fetchone()
    if x != None:
        sql = """ select mirisno_svojstvo
        from svojstva
        where id_molekule = {idm}
        """.format(idm = x[0])
        cursor.execute(sql)
        y = cursor.fetchall()
    else:
        try:
            x = get_compound_details_name(name) # testirati koji error dolazi kad se ne da smiles kod
            # y je popis deskriptora koje AI daje
            y = [] # za sad
        except Exception as e:
            #print("Online fail")
            sm = None
            y=[]
        
    if sm is not None:
        mol = Chem.MolFromSmiles(x[2])# provjerit dal u ovom trenutku postoji smiles kod, i nastavit
        if mol is not None:
            img = Draw.MolToImage(mol, size=(300, 300))
            img = img.convert('RGBA')

            imgdata = img.getdata()
            newimgdata = [(r, g, b, 0) if ((r > 100) and (g > 100) and (b > 100)) else (r, g, b, a) for r, g, b, a in imgdata]
            img.putdata(newimgdata)
    else:
        img = None
    
    
        
    return x,y,img
   

def getScent(molId, num):
    conn = psycopg2.connect(
        database="kemoinformatika",
        user=pgUser,
        password=pgPassword,
        host='127.0.0.1',
        port='5432'
        )
    cursor = conn.cursor()
    sql = """ select * from slicnostmiris
    where id_molekule1 = {id1} or id_molekule2 ={id2}
    order by slicnostm desc limit {n}
    """.format(id1 = molId, id2 = molId, n = num)
    
    cursor.execute(sql)
    x = cursor.fetchall()
    #print(x)
    ids = []
    for i in x:
        if i[0] != molId:
            ids.append(i[0])
        elif i[1] != molId:
            ids.append(i[1])
    mol = []
    for i in ids:
        sql="""Select m.ime_molekule, m.smiles_kod, s.slicnostm from molekule as m
        inner join slicnostmiris as s
        on m.id_molekule = s.id_molekule1 or m.id_molekule = s.id_molekule2
        where id_molekule = {idm} and (s.id_molekule1 = {id2} or s.id_molekule2 = {id3})""".format(idm = i, id2 = molId, id3 = molId)
        cursor.execute(sql)
        m = cursor.fetchone()
        mol.append(m)
    return mol




def getTanimoto(molId, num):
    conn = psycopg2.connect(
        database="kemoinformatika",
        user=pgUser,
        password=pgPassword,
        host='127.0.0.1',
        port='5432'
        )
    cursor = conn.cursor()
    sql = """ select * from slicnostsvojstva
    where id_molekule1 = {id1} or id_molekule2 ={id2}
    order by slicnosts desc limit {n}
    """.format(id1 = molId, id2 = molId, n = num)
    
    cursor.execute(sql)
    x = cursor.fetchall()
    #print(x)
    ids = []
    for i in x:
        if i[0] != molId:
            ids.append(i[0])
        elif i[1] != molId:
            ids.append(i[1])
    mol = []
    for i in ids:
        sql="""Select m.ime_molekule, m.smiles_kod, s.slicnosts from molekule as m
        inner join slicnostsvojstva as s
        on m.id_molekule = s.id_molekule1 or m.id_molekule = s.id_molekule2
        where id_molekule = {idm} and (s.id_molekule1 = {id2} or s.id_molekule2 = {id3})""".format(idm = i, id2 = molId, id3 = molId)
        cursor.execute(sql)
        m = cursor.fetchone()
        mol.append(m)
    return mol
        

#print(SearchDescriptor("fishy, nutty, woody"))
#print(SearchSmiles("CC1CCCC(N1)C"))
#print(SearchName("aspirin"))
#print(getScent(8082, 10))
#print(getTanimoto(8082, 10))
