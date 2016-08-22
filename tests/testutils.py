from random import choice


def seq_generator(length):
    """Random sequence generator."""
    return ''.join(choice('CGTA') for _ in range(length))
