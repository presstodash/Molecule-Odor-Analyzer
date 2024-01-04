import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

#training
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.model_selection import GridSearchCV

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
ecfp_fingerprints = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=fingerprint_size) for mol in mol_objects]

boolean_fingerprint_matrix = np.array([[bool(int(b)) for b in fp.ToBitString()] for fp in ecfp_fingerprints], dtype=bool)

boolean_characteristics_data = np.array([[bool(b) for b in row] for row in characteristics_data], dtype=bool)

# ----------------------------
# Model training
# ----------------------------

X_train, X_val, y_train, y_val = train_test_split(boolean_fingerprint_matrix, boolean_characteristics_data, test_size=0.2, random_state=42)

param_grid = {
    'n_estimators': [50, 100, 200],
    'max_depth': [None, 10, 20],
    'min_samples_split': [2, 5, 10],
    'min_samples_leaf': [1, 2, 4]
}

rf = RandomForestClassifier(random_state=42)

grid_search = GridSearchCV(estimator=rf, param_grid=param_grid, cv=5, scoring='accuracy')
grid_search.fit(X_train, y_train)

print("Best Hyperparameters:", grid_search.best_params_)

best_model = grid_search.best_estimator_

best_predictions = best_model.predict(X_val)
best_accuracy = accuracy_score(y_val, best_predictions)
print(f"Best Validation Accuracy: {best_accuracy}")

# ---------------------------
# Saving the model
# ---------------------------

joblib.dump(model, 'random_forest_model2.joblib')
