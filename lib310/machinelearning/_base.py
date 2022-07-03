from abc import ABC, abstractmethod


class BaseMLTask(ABC):
    def __init__(self, model, tokenizer, **kwargs):
        self.model = model
        self.tokenizer = tokenizer

    @abstractmethod
    def fit(self, *args, **kwargs):
        raise NotImplementedError

    @abstractmethod
    def run(self, *args, **kwargs):
        raise NotImplementedError
