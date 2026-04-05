"""logger.py – Centralised logging and console output for pynetworkin."""

from loguru import logger
from rich.console import Console
from rich.theme import Theme

# Remove the default loguru handler to avoid duplicate output.
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
    lambda msg: console.print(msg, end=""),
    format="<green>{time:HH:mm:ss}</green> | <level>{level: <8}</level> | <cyan>{name}</cyan>:<cyan>{line}</cyan> — <level>{message}</level>",
    colorize=True,
    level="INFO",
)
