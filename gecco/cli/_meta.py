import functools
import logging
import warnings

import verboselogs


class classproperty(object):
    """A class property decorator.
    """

    def __init__(self, f):
        self.f = f

    def __get__(self, obj, owner):
        return self.f(owner)


class BraceAdapter(logging.LoggerAdapter, verboselogs.VerboseLogger):
    """An logging adapter for `VerboseLogger` to use new-style formatting.
    """

    class Message(object):
        def __init__(self, fmt, args):
            self.fmt = str(fmt)
            self.args = args

        def __str__(self):
            return self.fmt.format(*self.args)

    def __init__(self, logger, extra=None):
        super(BraceAdapter, self).__init__(logger, extra or {})

    def log(self, level, msg, *args, **kwargs):
        if self.isEnabledFor(level):
            msg, kwargs = self.process(msg, kwargs)
            self.logger._log(level, self.Message(msg, args), (), **kwargs)

    def notice(self, msg, *args, **kwargs):
        if self.isEnabledFor(verboselogs.NOTICE):
            msg, kwargs = self.process(msg, kwargs)
            self.logger.notice(self.Message(msg, args), **kwargs)

    def spam(self, msg, *args, **kwargs):
        if self.isEnabledFor(verboselogs.SPAM):
            msg, kwargs = self.process(msg, kwargs)
            self.logger.spam(self.Message(msg, args), **kwargs)

    def verbose(self, msg, *args, **kwargs):
        if self.isEnabledFor(verboselogs.VERBOSE):
            msg, kwargs = self.process(msg, kwargs)
            self.logger.verbose(self.Message(msg, args), **kwargs)

    def success(self, msg, *args, **kwargs):
        if self.isEnabledFor(verboselogs.SUCCESS):
            msg, kwargs = self.process(msg, kwargs)
            self.logger.success(self.Message(msg, args), **kwargs)


def wrap_warnings(logger: logging.Logger):
    """Have the function patch `warnings.showwarning` with the given logger.

    Arguments:
        logger (~logging.logger): the logger to wrap warnings with when
            the decorated function is called.

    Returns:
        `function`: a decorator function.

    """

    def new_showwarning(message, category, filename, lineno, file=None, line=None):
        logger.warning(message)

    def decorator(func):
        @functools.wraps(func)
        def new_func(*args, **kwargs):
            old_showwarning = warnings.showwarning
            warnings.showwarning = new_showwarning
            try:
                return func(*args, **kwargs)
            finally:
                warnings.showwarning = old_showwarning

        return new_func

    return decorator
