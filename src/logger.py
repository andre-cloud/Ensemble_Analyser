import logging
import os

LOG_FORMAT = "%(message)s"

DEBUG = bool(os.getenv("DEBUG"))


def ordinal(n):
    return "%d-%s" % (n, "tsnrhtdd"[(n // 10 % 10 != 1) * (n % 10 < 4) * n % 10 :: 4])


def create_log(output):
    """
    Creating an logger instance.

    output | str : output filename

    return : logger instance
    """

    logging.basicConfig(
        filename=output,
        level=logging.DEBUG if DEBUG else logging.INFO,
        format=LOG_FORMAT,
        filemode="w",
    )

    log = logging.getLogger()
    # sys.stdout = StreamToLogger(log,logging.INFO)
    # sys.stderr = StreamToLogger(log,logging.ERROR)

    log.warning(f"DEBUG mode: {DEBUG}")

    return logging.getLogger()


if __name__ == "__main__":  # pragma: no cover:
    log = create_log("tests/output_test.out")
