import os
import sys
import logging


class StreamFormatter(logging.Formatter):

    error_fmt = "%(levelname)s - %(module)s, line %(lineno)d:\n%(message)s"
    warning_fmt = "%(levelname)s:\n%(message)s"
    debug_fmt = "%(asctime)s\n%(message)s"
    info_fmt = "%(message)s"

    def __init__(self):
        super().__init__(fmt="%(levelname)s: %(message)s",
                         datefmt="%d/%m %H:%M:%S")

    def format(self, record):

        # Save the original format configured by the user
        # when the logger formatter was instantiated
        format_orig = self._style._fmt

        # Replace the original format with one customized by logging level
        if record.levelno == logging.DEBUG:
            self._style._fmt = StreamFormatter.debug_fmt

        elif record.levelno == logging.INFO:
            self._style._fmt = StreamFormatter.info_fmt

        elif record.levelno == logging.ERROR:
            self._style._fmt = StreamFormatter.error_fmt

        elif record.levelno == logging.WARNING:
            self._style._fmt = StreamFormatter.warning_fmt

        # Call the original formatter class to do the grunt work
        result = logging.Formatter.format(self, record)

        # Restore the original format configured by the user
        self._style._fmt = format_orig

        return result


class FileFormatter(logging.Formatter):
    """Log file output rules (removes all colour characters)."""
    fmt = logging.Formatter(fmt="[%(asctime)s] %(levelname)s: %(message)s",
                            datefmt="%Y-%m-%d %H:%M:%S")

    default_fmt = logging.Formatter(fmt="[%(asctime)s] INFO: %(message)s",
                                    datefmt="%Y-%m-%d %H:%M:%S")
    debug_fmt = logging.Formatter(fmt="[%(asctime)s] DEBUG: %(message)s",
                                  datefmt="%Y-%m-%d %H:%M:%S")
    info_fmt = logging.Formatter(fmt="[%(asctime)s] INFO: %(message)s",
                                 datefmt="%Y-%m-%d %H:%M:%S")
    warn_fmt = logging.Formatter(fmt="[%(asctime)s] WARNING: %(message)s",
                                 datefmt="%Y-%m-%d %H:%M:%S")
    err_fmt = logging.Formatter(fmt="[%(asctime)s] ERROR: %(message)s",
                                datefmt="%Y-%m-%d %H:%M:%S")

    def format(self, record):
        if record.levelno >= logging.ERROR:
            return self.err_fmt.format(record)
        elif record.levelno >= logging.WARNING:
            return self.warn_fmt.format(record)
        elif record.levelno >= logging.INFO:
            return self.info_fmt.format(record)
        elif record.levelno >= logging.DEBUG:
            return self.debug_fmt.format(record)
        else:
            return self.default_fmt.format(record)


def logger_name() -> str:
    return __name__.split('.')[0]


def reactivate_log() -> None:
    ts_logger = logging.getLogger(logger_name())
    for h in ts_logger.handlers:
        if h.__class__.__name__ == "StreamHandler":
            h.propagate = True
    ts_logger.disabled = False
    ts_logger.is_silent = False
    return


def prep_logging(log_file=None, verbosity=False, silent=False, stream=sys.stderr) -> None:
    """
    Allows for multiple file handlers to be added to the root logger, but only a single stream handler.
    The new file handlers must be removed outside of this function explicitly

    :param log_file: Path to a file to write the TreeSAPP log
    :param verbosity: Whether debug-level information should be written (True) or not (False)
    :param stream: Which stream, sys.stdout or sys.stderr, should the console logger write to?
    :param silent: Suppresses writing to the output stream if True
    :return: None
    """
    if verbosity:
        stream_level = logging.DEBUG
    else:
        stream_level = logging.INFO

    # Detect whether handlers are already present and return if true
    ts_logger = logging.getLogger(logger_name())
    if silent:
        ts_logger.disabled = True
    ts_logger.setLevel(logging.DEBUG)
    if len(ts_logger.handlers):
        return

    # Set the console handler normally writing to stdout/stderr
    ts_stream_logger = logging.StreamHandler(stream=stream)
    ts_stream_logger.setLevel(stream_level)
    ts_stream_logger.terminator = ''
    ts_stream_logger.setFormatter(StreamFormatter())
    ts_logger.addHandler(ts_stream_logger)

    ts_logger.is_silent = False
    if silent:
        ts_logger.is_silent = True
        ts_stream_logger.setLevel(logging.ERROR)
        ts_logger.disabled = True

    if log_file:
        if not os.path.isabs(log_file):
            log_file = os.path.join(os.getcwd(), os.path.dirname(log_file), os.path.basename(log_file))
        log_dir = os.path.dirname(log_file)
        try:
            if log_dir and not os.path.isdir(log_dir):
                os.mkdir(log_dir)
        except (IOError, OSError):
            sys.stderr.write("ERROR: Unable to make directory '" + log_dir + "'.\n")
            sys.exit(3)
        ts_file_logger = logging.FileHandler(log_file, 'w')
        ts_file_logger.setFormatter(FileFormatter())
        ts_file_logger.setLevel(logging.DEBUG)
        ts_logger.addHandler(ts_file_logger)
        ts_logger.propagate = False

    return
