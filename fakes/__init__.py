from .env import check_dependencies

SYSTEM_DEPENDENCIES = {
    'sextractor (sex)': (
        ['sex', '--version'],
        lambda v: v.split()[2],
        '2.18.0'
    ),
    'psfex': (
        ['psfex', '--version'],
        lambda v: v.split()[2],
        '3.21.1'
    )
}

check_dependencies(SYSTEM_DEPENDENCIES)

from .fake import *
