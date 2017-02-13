"""
CBMPy: CBConfig module
======================
PySCeS Constraint Based Modelling (http://cbmpy.sourceforge.net)
Copyright (C) 2009-2016 Brett G. Olivier, VU University Amsterdam, Amsterdam, The Netherlands

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>

Author: Brett G. Olivier
Contact email: bgoli@users.sourceforge.net
Last edit: $Author: bgoli $ ($Id: CBConfig.py 561 2017-01-24 13:36:09Z bgoli $)

"""
## gets rid of "invalid variable name" info
# pylint: disable=C0103
## gets rid of "line to long" info
# pylint: disable=C0301
## use with caution: gets rid of module xxx has no member errors (run once enabled)
# pylint: disable=E1101


# preparing for Python 3 port
#from __future__ import division, print_function
#from __future__ import absolute_import
#from __future__ import unicode_literals

# release

__CBCONFIG__ = {'VERSION_MAJOR' : 0,
                'VERSION_MINOR' : 1,
                'VERSION_MICRO' : 0,
                'VERSION_STATUS' : 'alpha',
                'VERSION' : '',
                'DEBUG' : False
               }

def current_version():
    """
    Return the current CBMPy version as a string

    """
    return '{}.{}.{}'.format(__CBCONFIG__['VERSION_MAJOR'], __CBCONFIG__['VERSION_MINOR'], __CBCONFIG__['VERSION_MICRO'])

def current_version_tuple():
    """
    Return the current CBMPy version as a tuple (x, y, z)

    """
    return (__CBCONFIG__['VERSION_MAJOR'], __CBCONFIG__['VERSION_MINOR'], __CBCONFIG__['VERSION_MICRO'])

__CBCONFIG__['VERSION'] = current_version()


