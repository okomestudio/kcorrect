# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Python implementation of Yanny idlutils IDL procedures.

Based on yanny_read.pro in idlutil v4_10_7.  The idea is to do what
those IDL routines do to get kcorrect to run, so the codes are not
very well written and in fact very ugly.

TODO:

  THIS MODULE OFFER A SET OF FUNCTIONALITY TO MAKE KCORRECT WORK; NEED
  MUCH MORE ATTENTION BEFORE GENERIC USE.

  -- Make things work strictly by following Yanny structure format
  -- Better documentation
  -- Comment code to make understandable

REVISION HISTORY:

  Jul 31, 2012 - A bit of cleaning up for astropy integration.
  Dec 28, 2011 - Some minor clean up.
  Feb 1, 2006 - First alpha version.
"""
from __future__ import print_function, division


def _combine_broken_lines(strs):
    """Combines lines broken by backslash.

    Given a list of strings, this function concatenates one entry to
    the next if the last non-whitespace character of the first string
    is a backslash (\).  A modified list of strings will be returned.

    Parameters
    ----------
    strs : list of str
    """
    if not strs:
        return []
    ret = [strs.pop(0)]
    while strs:
        if ret[-1].endswith('\\'):
            ret[-1] = ''.join([ret[-1], strs.pop(0)])
        else:
            ret.append(strs.pop(0))
    return ret


class InvalidDefStr(Exception): pass
class InvalidDataStr(Exception): pass


class YannyStruct(dict):
    """Yanny structure

    This structure is used like a dict object.
    """

    def __init__(self, defstr):
        self.name = ''
        self.trans_func = []
        self.is_list = []
        self.element_names = []
        self.defstr = defstr
        self.define_from_Yanny_struct()
        super(dict, self).__init__()

    def define_from_Yanny_struct(self, defstr=''):
        """Construct structure from given string.

        Define a structure from the given string defining a Yanny
        structure format.
        """
        if not defstr:
            defstr = self.defstr
        defstr = defstr.strip()
        defstr = defstr.rstrip(';')

        j1 = defstr.find('{')
        s1, s2 = defstr[0:j1].split()
        if not (s1 == 'typedef' and s2 == 'struct'):
            raise InvalidDefStr('typedef is not struct.')
        j2 = defstr.rfind('}')
        self.name = defstr[j2+1:].strip()
        deftypes = defstr[j1+1:j2].strip().rstrip(';')

        # Define type conversion mapping to be used
        tfuncs = {'char': 'str', 'short': 'int', 'int': 'int',
                  'long': 'int', 'float': 'float', 'double': 'float'}

        for each in deftypes.split(';'):
            try:
                etype, ename = each.split()
            except ValueError:
                raise InvalidDefStr('Invalid data type definition.')

            # Python does not need a length for char array
            if etype.startswith('char['):
                etype = 'char'
            # Set an unknown type to 'char'
            if etype not in tfuncs.keys():
                etype = 'char'
            # Provide a conversion function
            self.trans_func.append(tfuncs[etype])

            # Remove length info from char name definition
            if etype == 'char':
                j1 = ename.rfind('[')
                if j1 > 0:
                    ename = ename[:j1]

            # Is this element an array? (Only 1-dimensional arrays
            # supported here.)
            j1 = ename.find('[')
            if j1 > 0:
                j2 = ename.find(']')
                if j2 < 0:
                    raise InvalidDefStr('Invalid data type definition.')
                length = int(ename[j1+1:j2])
                ename = ename[:j1]
                self.is_list.append(length)
            else:
                self.is_list.append(0)

            self.element_names.append(ename)
            self[ename] = []

    def parse_data_line(self, s):
        """Read in a set of parameters from a string."""
        tokens = s.split('"')
        if not len(tokens) % 2:
            raise InvalidDataStr('Open double quote block found.')

        elements = []
        for i in xrange(0, len(tokens), 2):
            words = tokens[i].replace('{', '').replace('}', '').split()
            for j in xrange(len(words)):
                words[j] = words[j].strip()
            if i == 0:
                elements.extend(words)
            else:
                elements.append(tokens[i-1])
                elements.extend(words)

        if elements[0].upper() != self.name:
            msg = 'Invalid data string for {0}.'.format(self.name)
            raise InvalidDataStr(msg)

        elements = elements[1:]
        elements = iter(elements)
        elem_no = 0
        while True:
            try:
                elem = elements.next()
            except StopIteration:
                break
            if not self.is_list[elem_no]:
                exec ('val = {0}("{1}")'.format(self.trans_func[elem_no], elem))
            else:
                val = []
                exec ('val.append({0}("{1}"))'
                      .format(self.trans_func[elem_no], elem))
                for i in xrange(self.is_list[elem_no] - 1):
                    elem = elements.next()
                    exec ('val.append({0}("{1}"))'
                          .format(self.trans_func[elem_no], elem))

            self[self.element_names[elem_no]].append(val)
            elem_no += 1


class YannyFileRead(object):
    """Read structures from a Yanny format file.

    Notes
    -----
    The enum feature has not been implemented (currently only reads in
    the portion of file pertaining to enum structure).
    """

    def __init__(self, file):
        self.file = file
        self.texts = {'enum': [], 'struct': [],
                      'header': [], 'data': []}
        self.structs = {}
        self.keyvalues = {}

        self._parse_file()
        self._parse_data()

    def _parse_data(self):
        for each in self.texts['struct']:
            defstr = ' '.join(_combine_broken_lines(each))
            o = YannyStruct(defstr)
            self.structs[o.name] = o

        names = self.structs.keys()
        for each in self.texts['data']:
            line = each.strip()
            if not line:
                continue
            tokens = each.split()
            name = tokens[0].upper()

            if name in names:
                self.structs[name].parse_data_line(each)
            else:
                # If not a struct data, it is a key-val pair.
                self.keyvalues[tokens[0]] = line[len(tokens[0]):].strip()

    def _parse_file(self):
        f = open(self.file)
        # Read header info if exist.
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith('#'):
                self.texts['header'].append(line)
            else:
                break
        # Read the rest.
        try:
            while True:
                if line.startswith('typedef'):
                    chunk = []
                    tokens = line.split()
                    if tokens[1] == 'struct':
                        chunk.append(line)
                        while not line.startswith('}'):
                            line = f.next().rstrip()
                            chunk.append(line)
                        self.texts['struct'].append(chunk)
                    elif tokens[1] == 'enum':
                        chunk.append(line)
                        while not line.startswith('}'):
                            line = f.next().rstrip()
                            chunk.append(line)
                        self.texts['enum'].append(chunk)
                elif line != '':
                    self.texts['data'].append(line)
                line = f.next().rstrip()
        except StopIteration:
            pass
        f.close()


def test():
    import os
    dir = '/usr/local/kcorrect/kcorrect/data/filters'
    for each in os.listdir(dir):
        if not each.endswith('.par'):
            continue
        d = YannyFileRead(dir +'/'+ each)
        print((each, d.structs.keys(),
               len(d.structs[d.structs.keys()[0]]['lambda'])))


if __name__ == '__main__':
    """test code"""
    test()
