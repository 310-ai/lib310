class VolumeLimitException(Exception):

    def __init__(self, message="The query uses volume more than our limit"):
        super().__init__(message)