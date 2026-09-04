# EvaluateDNN_ttDM.py

import os
import numpy as np

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

import tensorflow as tf
from tensorflow.keras.models import load_model
from tensorflow.keras.layers import Dense
from tensorflow.keras import regularizers
import joblib


# ---------------------------------------------------------
# Custom layer used by the ttDM DNN
# ---------------------------------------------------------

class AffineConditioning(tf.keras.layers.Layer):
    def __init__(self, units, l2_weight=1e-5, l2_bias=1e-6, **kwargs):
        super().__init__(**kwargs)

        self.units = units
        self.l2_weight = l2_weight
        self.l2_bias = l2_bias

        self.scale = Dense(
            units,
            activation=None,
            kernel_regularizer=regularizers.l2(l2_weight),
            bias_regularizer=regularizers.l2(l2_bias),
        )

        self.bias = Dense(
            units,
            activation=None,
            kernel_regularizer=regularizers.l2(l2_weight),
            bias_regularizer=regularizers.l2(l2_bias),
        )

    def build(self, input_shape):
        hidden_shape, mass_shape = input_shape

        self.scale.build(mass_shape)
        self.bias.build(mass_shape)

        super().build(input_shape)

    def call(self, inputs):
        hidden, mass = inputs
        return hidden * self.scale(mass) + self.bias(mass)

    def get_config(self):
        config = super().get_config()
        config.update({
            "units": self.units,
            "l2_weight": self.l2_weight,
            "l2_bias": self.l2_bias,
        })
        return config

MODEL_DIR = "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/Full2022EEv12/DNNmodels/Models"


# ---------------------------------------------------------
# Load model
# ---------------------------------------------------------

model = load_model(
    f"{MODEL_DIR}/model_DNN_ttDM_s.keras",
    custom_objects={"AffineConditioning": AffineConditioning},
    compile=False,
)


# ---------------------------------------------------------
# Load scaler
# ---------------------------------------------------------

scaler = joblib.load(f"{MODEL_DIR}/scaler_model_DNN_ttDM_s.pkl")

feature_names = scaler.feature_names_in_
assert feature_names[-1] == 'mPhi'

_scaler_mean = getattr(scaler, "mean_", None)
_scaler_scale = getattr(scaler, "scale_", None)

if _scaler_mean is not None and _scaler_scale is not None:
    _scaler_mean = np.asarray(_scaler_mean, dtype=np.float32)
    _scaler_inv_scale = 1.0 / np.asarray(_scaler_scale, dtype=np.float32)


@tf.function(reduce_retracing=True)
def _model_predict(inputs):
    return model(inputs, training=False)


def _scale_inputs(inputs):
    if _scaler_mean is not None and _scaler_scale is not None:
        inputs -= _scaler_mean
        inputs *= _scaler_inv_scale
        return inputs

    return scaler.transform(inputs)


def _predict(inputs):
    inputs_scaled = _scale_inputs(np.asarray(inputs, dtype=np.float32))
    return _model_predict(inputs_scaled).numpy().reshape(-1)


def load_neural_network_all_mPhi(inputs, masses):
    try:
        base_inputs = np.asarray(inputs, dtype=np.float32).reshape(-1)
        mass_points = np.asarray(masses, dtype=np.float32).reshape(-1)

        inputs_with_mass = np.empty(
            (len(mass_points), len(base_inputs) + 1),
            dtype=np.float32
        )

        inputs_with_mass[:, :-1] = base_inputs
        inputs_with_mass[:, -1] = mass_points

        result = _predict(inputs_with_mass)

        return result.tolist()

    except Exception:
        return [-99.9] * len(masses)


if __name__ == "__main__":
    test_input = np.arange(len(feature_names) - 1, dtype=np.float32)
    test_masses = [
        10, 50, 100, 150, 200, 250, 300,
        350, 400, 500, 600, 700, 800, 1000
    ]

    result = load_neural_network_all_mPhi(test_input, test_masses)

    print("Result:", result)

