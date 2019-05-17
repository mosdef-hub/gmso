import logging
import topology.utils.log
from topology.tests.base_test import BaseTest


class TestLogging(BaseTest):
    def test_open_write_close_log(self):
        topology.utils.log.start_logging()
        log = logging.getLogger("TopLog")
        msg = "Test debug message"
        log.debug(msg)
        topology.utils.log.end_logging()
        logfile = open('TopDebug.out', 'r').readlines()
        assert msg in logfile[0]
