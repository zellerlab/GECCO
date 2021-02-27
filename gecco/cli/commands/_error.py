
class InvalidArgument(ValueError):
    pass

class CommandExit(Exception):
    def __init__(self, code):
        self.code = code
