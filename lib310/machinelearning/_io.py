from typing import Optional

import numpy as np

from lib310.data import ProteinDataTable


class GoAnnotationOutput(object):
    def __init__(self, outputs, obs: ProteinDataTable):
        self.raw_outputs = outputs
        self.obs = obs

        self.n_samples = len(outputs)

        self.embeddings = np.array([output[0][0] for output in outputs]).reshape(self.n_samples, -1)

    def get_hidden_states(self, layer: Optional[int] = None):
        if layer is None:
            layer = -1

        return self.embeddings


class CausalLMOutput(object):
    def __init__(self, outputs):
        self.raw_outputs = outputs

        self.generated_proteins = [output['generated_text'].strip() for output in outputs]
