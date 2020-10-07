"""

Set of functions to access CHD/UKBB sample metadata

"""

import hail as hl


def get_sample_metadata() -> hl.Table:
    """
    Return HT with sample metadata

    :return: HailTable
    """
    ht_path = 'update_path'  # TODO
    return hl.read_table(path=ht_path, min_partitions=500)

