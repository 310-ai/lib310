import numpy as np
import pandas as pd


class ModelOutput(object):
    def __init__(self, outputs, obs: pd.DataFrame):
        self.raw_outputs = outputs
        self.obs = obs

        self.n_samples = len(outputs)

        self.embeddings = np.array([output[0][0] for output in outputs]).reshape(self.n_samples, -1)

    def get_hidden_states(self, layer: Optional[int] = None):
        if layer is None:
            layer = -1

        return self.embeddings
