import sys
import logging
from functools import wraps

def graceful_exit(main_func):
    @wraps(main_func)
    def wrapper(*args, **kwargs):
        try:
            return main_func(*args, **kwargs)
        except KeyboardInterrupt:
            sys.exit(130)
        except Exception as e:
            logger = logging.getLogger('pipeline')
            if logger.handlers:
                logger.critical("A critical, unhandled exception occurred. The programme will now terminate.")
                logger.exception(e)
            else:
                sys.exit(1)

    return wrapper