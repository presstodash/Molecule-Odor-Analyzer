import psycopg2  # PostgreSQL
import pubchempy as pcp  # PubChem baza podataka
import csv

file_path = 'curated_GS_LF_merged_4983.csv'

firstLine = True

pgUser = "postgres"
pgPassword = "root"


x = 1



def get_compound_details(smiles_code):
    compound = pcp.get_compounds(smiles_code, 'smiles')[0]  # Fetch the first compound from the search results
    compound_id = compound.cid  # Compound ID
    name = compound.iupac_name  # Compound name
    formula = compound.molecular_formula  # Compound formula
    weight = compound.molecular_weight  # Molecular weight
    return compound_id, name, formula, weight





with open(file_path, 'r', newline='') as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        if firstLine:
            firstLine = False
        else:
            try:
                
                conn = psycopg2.connect(
                    database="kemoinformatika",
                    user=pgUser,
                    password=pgPassword,
                    host='127.0.0.1',
                    port='5432'
                    )
                cursor = conn.cursor()
                smiles = row[0]  # Prvi stupac
                #print(smiles)
                sql = "SELECT * FROM Molekule WHERE smiles_kod like '{sm}'".format(sm = smiles)
                cursor.execute(sql)
                mol = cursor.fetchone()
                
                if mol == None:
                    
                    scents = row[1].split(";")  # Drugi stupac
                    compound_id, name, formula, weight = get_compound_details(smiles)
                    
                    if name == None:
                        name = ""
                    if compound_id is None:
                        continue
                
                    name = name.replace("'","''")
                    sql = "INSERT INTO Molekule(ID_Molekule, Ime_Molekule, SMILES_Kod, Kemijska_Formula, Molekulska_Masa) VALUES ({cid}, '{ime}', '{smi}', '{kf}', {mm})".format(cid=compound_id, ime=name, smi=smiles, kf=formula, mm=weight)
                    #print(sql)
                    cursor.execute(sql)
                    conn.commit()
                
                    for i in scents:
                        sql = "INSERT INTO SVOJSTVA(ID_Molekule, Mirisno_Svojstvo) VALUES({cid}, '{sv}')".format(cid = compound_id, sv = i)
                        cursor.execute(sql)
                        conn.commit()
                        #print(compound_id)
                    
                    conn.close()
            except Exception as e:
                print(x)
                print(e)
                print(smiles)
                x = x + 1

