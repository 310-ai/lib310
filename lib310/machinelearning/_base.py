from abc import ABC, abstractmethod


class BaseMLTask(ABC):
    def __init__(self, *args, **kwargs):
        pass

    @abstractmethod
    def fit(self, *args, **kwargs):
        raise NotImplementedError

    @abstractmethod
    def run(self, *args, **kwargs):
        raise NotImplementedError

    def __get_embeddings(self, *args, **kwargs):
        raise NotImplementedError

    @abstractmethod
    def save(self, *args, **kwargs):
        raise NotImplementedError
