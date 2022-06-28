import logging
from typing import Union

from rich.console import Console
from rich.logging import RichHandler

b10_logger = logging.getLogger("lib310")


class B10Config:
    def __init__(self, verbosity: int = logging.INFO):
        self.verbosity = verbosity

    def reset_logging_handler(self):
        """
        Resets "lib310" log handler to a basic RichHandler().
        This is useful if piping outputs to a file.
        """
        b10_logger.removeHandler(b10_logger.handlers[0])
        ch = RichHandler(level=self._verbosity, show_path=False, show_time=False)
        formatter = logging.Formatter("%(message)s")
        ch.setFormatter(formatter)
        b10_logger.addHandler(ch)

    @property
    def verbosity(self):
        return self._verbosity

    @verbosity.setter
    def verbosity(self, level: Union[int, str]):
        self._verbosity = level
        b10_logger.setLevel(level)
        if len(b10_logger.handlers) == 0:
            console = Console(force_terminal=True)
            if console.is_jupyter is True:
                console.is_jupyter = False
            ch = RichHandler(
                level=level, show_path=False, console=console, show_time=False
            )
            formatter = logging.Formatter("%(message)s")
            ch.setFormatter(formatter)
            b10_logger.addHandler(ch)
        else:
            b10_logger.setLevel(level)


settings = B10Config()
