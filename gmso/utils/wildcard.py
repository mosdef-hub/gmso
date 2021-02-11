from typing import Iterator
from itertools import product

import numpy as np

from gmso.utils._constants import FF_TOKENS_SEPARATOR

__all__ = [
    'WildCardTokenizer'
]


class WildCardTokenizer:
    def __init__(self, token: str) -> None:
        self.token = token
        self.tokens_chain = []
        self._initialize()

    def _initialize(self) -> None:
        self.tokens_chain.append(self.token)
        max_len = len(self.token)
        self.tokens_chain = list(f'{self.token[0:max_len-j]}'
                                 if j == 0 else f'{self.token[0:max_len-j]}*'
                                 for j in range(0, max_len+1))

    @staticmethod
    def group_tokens(*tokenizers: Iterator['WildCardTokenizer']) -> str:
        tokens_chains_gen = (tokenizer.tokens_chain for tokenizer in tokenizers)
        step = np.prod([len(tokenizer.tokens_chain) for tokenizer in tokenizers[1:]])
        possible_tokens = list(product(*tokens_chains_gen))
        j = 0
        yielded = set()
        for k in range(0, len(tokenizers[0].tokens_chain)):
            print('here', k)
            while True:
                try:
                    idx = j*step + k
                    j += 1
                    if idx not in yielded:
                        yield FF_TOKENS_SEPARATOR.join(possible_tokens[j*step+k])
                    else:
                        continue
                except IndexError:
                    j = 0
                    break

