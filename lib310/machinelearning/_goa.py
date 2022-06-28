from ._base import BaseMLTask


class GoAnnotation(BaseMLTask):
    def __init__(self):
        super().__init__()

    def fit(self, *args, **kwargs):
        raise NotImplementedError

    def run(self, *args, **kwargs):
        pass

    def __get_embeddings(self, *args, **kwargs):
        pass

    def save(self, *args, **kwargs):
        pass
