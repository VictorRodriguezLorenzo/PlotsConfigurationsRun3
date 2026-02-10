# EvaluateDNN.py

import numpy as np
from tensorflow.keras.models import load_model
import os
os.environ["TF_USE_LEGACY_KERAS"] = "1"

# Load the model
model = load_model('/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/Full2022v12/DNNmodels/Models/model_dnn_model_DNN_.h5', compile = False)

def load_neural_network(inputs):
    try:
        result = model.predict(np.array(inputs).reshape(1, -1), verbose = 0)
        return [result[0][0]]
    except Exception as e:
        return [-99.9]

if __name__ == "__main__":
    # Including some test code here for local testing
    test_input = [
    1, 2, 3, 4, 5, 6, 7,
    ]
    result = load_neural_network(test_input)
    print("Result:", result)
