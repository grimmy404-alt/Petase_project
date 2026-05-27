import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, classification_report
import pickle

# Load features

X = np.load("data/features.npy")
y = np.load("data/labels.npy")

# Split data

X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42
)

print(f"Training set: {len(X_train)}")
print(f"Test set: {len(X_test)}")

# Train model

print("\nTraining Random Forest...")
model = RandomForestClassifier(n_estimators=100, random_state=42)
model.fit(X_train, y_train)

# Evaluate

y_pdn = model.predict(X_test)
accuracy = accuracy_score(y_test, y_pdn)

print(f"\nAccuracy: {accuracy:.2%}")
print("\nClassification Report:")
print(classification_report(y_test, y_pdn, 
                          target_names=["Non-PETase", "PETase"]))

# Save model

with open("models/petase_classifier.pkl", "wb") as f:
    pickle.dump(model, f)

print("\nModel saved to models/petase_classifier.pkl")