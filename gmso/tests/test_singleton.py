from gmso.utils.singleton import Singleton
from gmso.tests.base_test import BaseTest


class TestSingleton(BaseTest):
    def test_unique_id(self):
        singletons = [id(Singleton()) for i in range(10)]
        assert len(set(singletons)) == 1
