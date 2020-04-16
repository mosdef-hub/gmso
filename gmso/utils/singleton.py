class Singleton(object):
    """A helper super class to create singletons.

    A singleton class is a class for which no more than one
    instance can exist. Any class subclassing from Singleton
    will have only one instance.
    """
    _inst = None

    def __new__(cls):
        if cls._inst is None:
            cls._inst = super(Singleton, cls).__new__(cls)
        return cls._inst
