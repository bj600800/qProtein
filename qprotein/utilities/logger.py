"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/03/02

# Description: logging class with stream handler for the project.
# ------------------------------------------------------------------------------
"""

import logging


class LevelFilter:
    def __init__(self, handler_level):
        self.handler_level = handler_level

    def filter(self, record):
        return record.levelno <= self.handler_level


def setup_log(name):
    # define format for handlers
    debug_format = '%(asctime)s [%(levelname)s]' \
                     ' - FileName: %(filename)s' \
                     ' - FuncName: %(funcName)s' \
                     ' - Line number: %(lineno)d' \
                     ' - [%(message)s]'

    info_format = '%(asctime)s [%(levelname)s]: %(message)s'

    warning_format = '%(asctime)s [%(levelname)s]' \
                     ' - FileName: %(filename)s' \
                     ' - FuncName: %(funcName)s' \
                     ' - Line number: %(lineno)d' \
                     ' - [%(message)s]'

    # define handler for debug level
    debug_handler = logging.StreamHandler()
    debug_handler.setLevel(logging.DEBUG)
    debug_handler.setFormatter(logging.Formatter(debug_format))
    debug_handler.addFilter(LevelFilter(logging.DEBUG))

    # define handler for info level
    info_handler = logging.StreamHandler()
    info_handler.setLevel(logging.INFO)
    info_handler.setFormatter(logging.Formatter(info_format))
    info_handler.addFilter(LevelFilter(logging.INFO))

    # define handler for warning level
    warning_handler = logging.StreamHandler()
    warning_handler.setLevel(logging.WARNING)
    warning_handler.setFormatter(logging.Formatter(warning_format))
    warning_handler.addFilter(LevelFilter(logging.WARNING))

    # define handler for error and critical
    error_handler = logging.StreamHandler()
    error_handler.setLevel(logging.ERROR)
    error_handler.setFormatter(logging.Formatter(warning_format))

    # initiation logger
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    logger.addHandler(debug_handler)
    logger.addHandler(info_handler)
    logger.addHandler(warning_handler)
    logger.addHandler(error_handler)

    return logger
