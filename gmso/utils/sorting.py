"""Sorting utilities for alphanumeric strings."""
import re


def _atoi(text):
    """Convert a string to an int."""
    return int(text) if text.isdigit() else text


def natural_sort(text):
    """Given an alphanumeric string, sort using the natural sort algorithm."""
    return [_atoi(a) for a in re.split(r"(\d+)", text)]
