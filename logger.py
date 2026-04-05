"""logger.py — Centralised logging and console output for pynetworkin."""

from loguru import logger
from rich.console import Console
from rich.theme import Theme

logger.remove()

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

logger.add(
    lambda msg: console.print(str(msg), end=""),
    format="{time:HH:mm:ss} | {level: <8} | {name}:{line} - {message}",
    colorize=False,
    level="INFO",
)