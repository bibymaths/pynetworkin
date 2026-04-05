from loguru import logger as _logger
from rich.console import Console
from rich.theme import Theme

# Remove default loguru handler
_logger.remove()

THEME = Theme({
    "info": "bold cyan",
    "success": "bold green",
    "warning": "bold yellow",
    "error": "bold red",
    "header": "bold white on dark_blue",
    "muted": "dim",
    "score": "bold magenta",
})

console = Console(theme=THEME, highlight=True)

# Route loguru -> rich
_logger.add(
    lambda msg: console.print(str(msg), end=""),
    format="{time:HH:mm:ss} | {level: <8} | {message}",
    level="INFO",
)

# ---- Wrapper class ----

class Logger:
    def info(self, msg, *args):
        console.print(f"[info]{msg.format(*args)}[/]")

    def success(self, msg, *args):
        console.print(f"[success]{msg.format(*args)}[/]")

    def warning(self, msg, *args):
        console.print(f"[warning]{msg.format(*args)}[/]")

    def error(self, msg, *args):
        console.print(f"[error]{msg.format(*args)}[/]")

    def header(self, msg, *args):
        console.print(f"[header] {msg.format(*args)} [/]\n")

    def muted(self, msg, *args):
        console.print(f"[muted]{msg.format(*args)}[/]")

    def score(self, msg, *args):
        console.print(f"[score]{msg.format(*args)}[/]")


# This is what your code imports
logger = Logger()