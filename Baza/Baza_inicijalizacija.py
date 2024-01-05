import psycopg2  # PostgreSQL
#import pubchempy  # PubChem baza podataka

pgUser = "postgres"
pgPassword = "root"

# povezivanje s postgres bazom
conn = psycopg2.connect(
    database="postgres",
    user=pgUser,
    password=pgPassword,
    host='127.0.0.1',
    port='5432'
)
conn.autocommit = True

# kreiranje kurskora
cursor = conn.cursor()

cursor.execute("DROP DATABASE IF EXISTS Kemoinformatika")

# Kreiranje baze podataka
sql = "CREATE DATABASE Kemoinformatika"

cursor.execute(sql)
print("Baza podataka uspješno napravljena.")

# Odspajanje od postgres baze
conn.close()


# Spajanje na novu bazu podataka
conn = psycopg2.connect(
    database="kemoinformatika",
    user=pgUser,
    password=pgPassword,
    host='127.0.0.1',
    port='5432'
)
cursor = conn.cursor()

cursor.execute("DROP TABLE IF EXISTS Molekule")

sql = """
CREATE TABLE Molekule(
    ID_Molekule INTEGER NOT NULL PRIMARY KEY,
    Ime_Molekule VARCHAR(150),
    SMILES_Kod VARCHAR(150),
    Kemijska_Formula VARCHAR(50),
    Molekulska_Masa DECIMAL
)
"""

cursor.execute(sql)
print("Tablica Molekule uspješno napravljena.")
conn.commit()

cursor.execute("DROP TABLE IF EXISTS Svojstva")

sql = """
CREATE TABLE Svojstva(
    ID_Molekule INTEGER NOT NULL,
    Mirisno_Svojstvo VARCHAR(20),
    PRIMARY KEY (ID_Molekule, Mirisno_Svojstvo),
    FOREIGN KEY (ID_Molekule) REFERENCES Molekule(ID_Molekule)
)
"""

cursor.execute(sql)
print("Tablica Svojstva uspješno napravljena.")
conn.commit()

conn.close()
