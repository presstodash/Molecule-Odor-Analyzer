import psycopg2  # PostgreSQL
import pubchempy as pcp  # PubChem baza podataka
import csv
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.DataStructs import TanimotoSimilarity
from Molecule_Odor_Analyzer.Baza.Baza_account import getUser, getPassword

pgUser = getUser()
pgPassword = getPassword()

file_path = 'Molecule_Odor_Analyzer\Baza\curated_GS_LF_merged_4983.csv'

firstLine = True



x = 1



def get_compound_details(smiles_code):
    compound = pcp.get_compounds(smiles_code, 'smiles')[0]  # Fetch the first compound from the search results
    compound_id = compound.cid  # Compound ID
    name = compound.iupac_name  # Compound name
    formula = compound.molecular_formula  # Compound formula
    weight = compound.molecular_weight  # Molecular weight
    smiles = compound.canonical_smiles
    return compound_id, name, formula, weight, smiles





with open(file_path, 'r', newline='') as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        if firstLine:
            firstLine = False
        else:
            try:
                print(x)
                x +=1
                
                conn = psycopg2.connect(
                    database="kemoinformatika",
                    user=pgUser,
                    password=pgPassword,
                    host='127.0.0.1',
                    port='5432'
                    )
                cursor = conn.cursor()
                smiles = row[0].upper()  # Prvi stupac
                #print(smiles)
                sql = "SELECT * FROM Molekule WHERE smiles_kod like '{sm}'".format(sm = smiles)
                cursor.execute(sql)
                mol = cursor.fetchone()
                
                if mol == None:
                    
                    scents = row[1].split(";")  # Drugi stupac
                    compound_id, name, formula, weight, smiles = get_compound_details(smiles)
                    
                    if name == None:
                        name = ""
                    if compound_id is None:
                        continue
                
                    name = name.replace("'","''")
                    sql = "INSERT INTO Molekule(ID_Molekule, Ime_Molekule, SMILES_Kod, Kemijska_Formula, Molekulska_Masa) VALUES ({cid}, '{ime}', '{smi}', '{kf}', {mm})".format(cid=compound_id, ime=name.lower(), smi=smiles, kf=formula, mm=weight)
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
                print(e)
                print(smiles)
                
                
                

conn = psycopg2.connect(
    database="kemoinformatika",
    user=pgUser,
    password=pgPassword,
    host='127.0.0.1',
    port='5432'
    )
cursor = conn.cursor()

sql = "SELECT ID_Molekule from Molekule"
cursor.execute(sql)
ids = cursor.fetchall()
conn.close()




#print(ids[0][0])


for i in range(len(ids)):
    for j in range(i + 1, len(ids)):
        try:
            conn = psycopg2.connect(
                database="kemoinformatika",
                user=pgUser,
                password=pgPassword,
                host='127.0.0.1',
                port='5432'
                )
            cursor = conn.cursor()
            sql = "SELECT SMILES_Kod from Molekule where ID_Molekule = {idm}".format(idm=ids[i][0])
            cursor.execute(sql)
            sm1 = cursor.fetchone()[0]
            sql = "SELECT SMILES_Kod from Molekule where ID_Molekule = {idm}".format(idm=ids[j][0])
            cursor.execute(sql)
            sm2 = cursor.fetchone()[0]
            
            mol1 = Chem.MolFromSmiles(sm1)
            mol2 = Chem.MolFromSmiles(sm2)
            
            if mol1 is not None and mol2 is not None:
                fp1 = AllChem.RDKFingerprint(mol1)
                fp2 = AllChem.RDKFingerprint(mol2)
                
                similarity = TanimotoSimilarity(fp1, fp2)
                
                sql = "INSERT INTO SlicnostSvojstva(ID_Molekule1, ID_Molekule2, SlicnostS) VALUES({cid1}, {cid2}, {sl})".format(cid1=ids[i][0], cid2=ids[j][0], sl=similarity)
                cursor.execute(sql)
                conn.commit()
                
                sql = "SELECT Mirisno_Svojstvo from Svojstva where ID_Molekule = {idm}".format(idm=ids[i][0])
                cursor.execute(sql)
                mol1 = cursor.fetchall()
                sql = "SELECT Mirisno_Svojstvo from Svojstva where ID_Molekule = {idm}".format(idm=ids[j][0])
                cursor.execute(sql)
                mol2 = cursor.fetchall()
                
                list1 = [k[0] for k in mol1]
                list2 = [k[0] for k in mol2]

                total = len(list1)
                match = 0
                for i in list2:
                    if i in list1:
                        match += 1
                    else :
                        total += 1
                
                if match > 0:
                    sql = "INSERT INTO SlicnostMiris(ID_Molekule1, ID_Molekule2, SlicnostM) VALUES({cid1}, {cid2}, {sl})".format(
                        cid1=ids[i][0], cid2=ids[j][0], sl=match / total)
                    cursor.execute(sql)
                    conn.commit()

            conn.close()
        except Exception as e:
            print(e)
