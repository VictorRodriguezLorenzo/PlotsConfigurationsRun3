# EvaluateDNN.py

import numpy as np
from tensorflow.keras.models import load_model
import os
import joblib
os.environ["TF_USE_LEGACY_KERAS"] = "1"

# Load the model
model = load_model('/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/Full2023v12/DNNmodels/Models/model_dnn_model.h5', compile = False)

# Load scaler
scaler = joblib.load(
    "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/Full2023v12/DNNmodels/Models/scaler_dnn_model.pkl"
)

feature_names = scaler.feature_names_in_

def load_neural_network(inputs):
    try:
        df = pd.DataFrame([inputs], columns=feature_names)
        inputs_scaled = scaler.transform(df)
        result = model.predict(inputs_scaled, verbose=0)
        print(float(result[0][0]))
        return [float(result[0][0])]
    except Exception as e:
        return [-99.9]

if __name__ == "__main__":
    # Including some test code here for local testing
    test_input = [
    1, 2, 3, 4, 5, 6, 7,
    ]
    result = load_neural_network(test_input)
    print("Result:", result)
