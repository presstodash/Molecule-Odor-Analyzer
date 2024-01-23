import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

#training
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.metrics import roc_auc_score


#saving
import joblib

# Preprocessing

csv_path = 'curated_GS_LF_merged_4983.csv'  
df = pd.read_csv(csv_path)


smiles_data = df['nonStereoSMILES'].values  

#krecemo od treceg stupca, prvi je smiles ,drugi sadrzi tekst zapis deskriptora, zelimo brojeve
characteristics_data = df.iloc[:, 2:].values  

mol_objects = [Chem.MolFromSmiles(smiles) for smiles in smiles_data]

fingerprint_size = 2048
maccs_fingerprints = [AllChem.GetMACCSKeysFingerprint(mol) for mol in mol_objects]

boolean_fingerprint_matrix = np.array([[bool(int(b)) for b in fp.ToBitString()] for fp in maccs_fingerprints], dtype=bool)

boolean_characteristics_data = np.array([[bool(b) for b in row] for row in characteristics_data], dtype=bool)

# ----------------------------
# Model training
# ----------------------------

X_train, X_val, y_train, y_val = train_test_split(boolean_fingerprint_matrix, boolean_characteristics_data, test_size=0.2, random_state=42)

model = RandomForestClassifier(n_estimators=128, random_state=42)

model.fit(X_train, y_train)

predictions = model.predict(X_val)

accuracy = accuracy_score(y_val, predictions)
print(f"Validation Accuracy: {accuracy}")


predictions_probability = model.predict_proba(X_val)[:, 1]
auroc = roc_auc_score(y_val, predictions_probability)
print(f"AUROC Score: {auroc}")

# ---------------------------
# Saving the model
# ---------------------------

joblib.dump(model, 'random_forest_modelmaccs.joblib')
