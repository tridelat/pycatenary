import warnings
from functools import wraps


def deprecated(message: str):
    """Decorator to mark functions as deprecated."""

    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            warnings.warn(message, DeprecationWarning, stacklevel=2)
            return func(*args, **kwargs)

        return wrapper

    return decorator
