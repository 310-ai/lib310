import os.path
from typing import Optional

import torch.cuda

from ._base import BaseMLTask
from transformers import AutoModelForCausalLM, AutoTokenizer, TextGenerationPipeline

from ._io import CausalLMOutput

SUPPORTED_BERT_MODELS = ['ProtGPT2']
MODELS_MAP = {
    'ProtGPT2': 'nferruz/ProtGPT2'
}


class AutoRegressiveLanguageModeling(BaseMLTask):
    def __init__(self, model, tokenizer, use_gpu: Optional[bool] = False):
        super().__init__(model, tokenizer)
        self.use_gpu = use_gpu and torch.cuda.is_available()
        self.pipeline = None

    @classmethod
    def from_pretrained(cls, model: str, version: str, use_gpu: Optional[bool] = False):
        model_name = model

        if model_name in SUPPORTED_BERT_MODELS:
            if os.path.exists(version):
                model = AutoModelForCausalLM.from_pretrained(os.path.abspath(version))
            else:
                model = AutoModelForCausalLM.from_pretrained(MODELS_MAP[model_name])
            tokenizer = AutoTokenizer.from_pretrained(MODELS_MAP[model_name])
        else:
            raise Exception('Invalid model or version passed!')

        return cls(model, tokenizer, use_gpu=use_gpu)

    def fit(self, *args, **kwargs):
        pass

    def run(self, num_samples: int = 32, batch_size: int = 32, **kwargs):
        if self.pipeline is None:
            self.pipeline = TextGenerationPipeline(model=self.model, tokenizer=self.tokenizer, batch_size=batch_size,
                                                   device=0 if self.use_gpu else -1, task='LM')
        outputs = self.pipeline('', max_length=kwargs.pop('max_length', 512), do_sample=True,
                                top_k=kwargs.pop('top_k', 950),
                                repetition_penalty=kwargs.pop('repitition_penalty', 1.2),
                                num_return_sequences=num_samples,
                                eos_token_id=kwargs.pop('eos_token_id', 0), **kwargs)

        return CausalLMOutput(outputs)
