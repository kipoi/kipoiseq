import abc

from kipoiseq import Interval
from kipoiseq.transforms.functional import rc_dna

import logging

log = logging.getLogger(__name__)

__all__ = [
    "BaseExtractor"
]


class BaseExtractor(object):
    __metaclass__ = abc.ABCMeta

    _use_strand: bool

    # main method
    @abc.abstractmethod
    def extract(self, interval: Interval, *args, **kwargs) -> str:
        raise NotImplementedError

    @property
    def use_strand(self):
        return self._use_strand

    # closing files
    def __del__(self):
        return self.close()

    def close(self):
        # implemented by the subclass
        pass
