class Singleton(object):
    _inst = None

    def __new__(cls):
        if cls._inst is None:
            cls._inst = super(Singleton, cls).__new__(cls)
        return cls._inst
