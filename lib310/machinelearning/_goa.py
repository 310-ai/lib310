import os

from typing import Optional

from tqdm import tqdm
import torch
from transformers.pipelines.base import Dataset

from ._base import BaseMLTask
from transformers import AutoTokenizer, BertForSequenceClassification

from ._io import GoAnnotationOutput
from transformers import FeatureExtractionPipeline

SUPPORTED_BERT_MODELS = ['prot_bert']
MODELS_MAP = {
    'prot_bert': 'Rostlab/prot_bert'
}


class GoADataset(Dataset):
    def __init__(self, sequences):
        self.seqs = sequences

    def __getitem__(self, index):
        return self.seqs[index]

    def __len__(self):
        return len(self.seqs)


class GoAnnotation(BaseMLTask):
    def __init__(self, model, tokenizer, num_labels, use_gpu):
        super().__init__(model, tokenizer)
        self.num_labels = num_labels
        self.use_gpu = use_gpu and torch.cuda.is_available()
        self.pipeline = None

    @classmethod
    def from_pretrained(cls, model: str, version: str, num_labels: int, use_gpu: Optional[bool] = True):
        model_name = model

        if model_name in SUPPORTED_BERT_MODELS:
            if os.path.exists(version):
                model = BertForSequenceClassification.from_pretrained(os.path.abspath(version), num_labels=num_labels)
            else:
                model = BertForSequenceClassification.from_pretrained(model_name, num_labels=num_labels)
            tokenizer = AutoTokenizer.from_pretrained(MODELS_MAP[model_name])
        else:
            raise Exception('Invalid model or version passed!')

        return cls(model, tokenizer, num_labels, use_gpu=use_gpu)

    def fit(self, *args, **kwargs):
        raise NotImplementedError

    def run(self, dataset, labels: Optional[list] = None,
            return_last_hidden_state: Optional[bool] = True, **kwargs):

        sequences = list(dataset['Sequence'].values)
        sequences = [' '.join(seq) for seq in sequences]

        goa_dataset = GoADataset(sequences)

        # def database(seqs):
        #     encoded_sequence = self.tokenizer(seqs, padding=True, return_tensors='pt')
        #     encoded_sequence = {key: value.to(self.model.device) for key, value in encoded_sequence.items()}
        #     return encoded_sequence
        #
        # model_outputs = []
        # batches = [sequences[i: i + batch_size] for i in range(0, len(sequences), batch_size)]
        # for batch in tqdm(batches):
        #     encoded = database(batch)
        #     with torch.no_grad():
        #         model_outputs.append(self.model(**encoded, output_hidden_states=return_last_hidden_state))
        if self.pipeline is None:
            self.pipeline = FeatureExtractionPipeline(model=self.model.bert, tokenizer=self.tokenizer,
                                                      task='GoA', device=0 if self.use_gpu else -1)
        for output in tqdm(self.pipeline(goa_dataset, batch_size=kwargs.get('batch_size', 32), total=len(goa_dataset))):
            pass

        return GoAnnotationOutput(outputs=output, obs=dataset.copy())
