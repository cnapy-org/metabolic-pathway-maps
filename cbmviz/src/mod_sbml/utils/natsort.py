# ---------------------------------------------------------
# natsort.py: Natural string sorting.
# ---------------------------------------------------------

# By Seo Sanghyeon.  Some changes by Connelly Barnes.

import re
from functools import cmp_to_key


def try_int(s):
    """Convert to integer if possible."""
    try:
        return int(s)
    except:
        return s


def natsort_key(s):
    """Used internally to get a tuple by which s is sorted."""
    # numbers sort before text (as in Python 2) and are never compared with text
    return [(0, part, '') if isinstance(part, int) else (1, 0, part)
            for part in map(try_int, re.findall(r'(\d+|\D+)', s))]


def natcmp(a, b):
    """Natural string comparison, case sensitive."""
    key_a, key_b = natsort_key(a), natsort_key(b)
    return (key_a > key_b) - (key_a < key_b)


def natcasecmp(a, b):
    """Natural string comparison, ignores case."""
    return natcmp(a.lower(), b.lower())


def natsorted(seq, cmp=natcmp):
    """Returns a list, sorted by natural string sort."""
    return sorted(seq, key=cmp_to_key(cmp))
