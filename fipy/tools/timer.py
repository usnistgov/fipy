import time

class Timer(object):
    """Context manager that measures time elapsed in context

    Defaults to nanosecond precision (although probably only microsecond or
    even millisecond accuracy).

    >>> with Timer() as timer:
    ...     pass
    >>> print("elapsed: {elapsed} ns".format(elapsed=timer.elapsed)) # doctest: +ELLIPSIS
    elapsed: ... ns

    Parameters
    ----------
    timer : callable, optional
        Function that returns a time

            ``timer() -> int or float``

        The difference between successive calls to ``timer()`` should give
        the elapsed time at the desired resolution and type.
        (default: :func:`~time.perf_counter_ns`)
    """
    def __init__(self, timer=None):
        self.timer = timer or getattr(time, "perf_counter_ns", self.clock_ns)

        self.start_time = self.stop_time = 0
        self.running = False

    @staticmethod
    def clock_ns():
        """Substitute "nanosecond" timer for Python 2.7
        """
        return int(time.clock() * 1e9)

    @property
    def elapsed(self):
        """Time measured so far
        """
        if self.running:
            return self.timer() - self.start_time
        else:
            return self.stop_time - self.start_time

    def __enter__(self):
        self.running = True
        self.start_time = self.stop_time = self.timer()

        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.stop_time = self.timer()
        self.running = False

        return False

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
