import bz2
import cPickle
import gzip
import json
import io
import yaml

# default buffer for incremental read/write
DEF_BUFFER = 16384


def save_object(file_name, obj):
    """
    Serialize an object to a file with gzip compression. .gz will automatically be
    added if missing.

    :param file_name: output file name
    :param obj: object to serialize
    """
    with open_output(file_name, compress='gzip') as out_h:
        cPickle.dump(obj, out_h)


def load_object(file_name):
    """
    Deserialize an object from a file with automatic support for compression.

    :param file_name: input file name
    :return: deserialzied object
    """
    with open_input(file_name) as in_h:
        return cPickle.load(in_h)


def open_input(file_name):
    """
    Open a text file for input. The filename is used to indicate if it has been
    compressed. Recognising gzip and bz2.

    :param file_name: the name of the input file
    :return: open file handle, possibly wrapped in a decompressor
    """
    suffix = file_name.split('.')[-1].lower()
    if suffix == 'bz2':
        return bz2.BZ2File(file_name, 'r')
    elif suffix == 'gz':
        return gzip.GzipFile(file_name, 'r')
    else:
        return open(file_name, 'r')


def open_output(file_name, append=False, compress=None, gzlevel=6):
    """
    Open a text stream for reading or writing. Compression can be enabled
    with either 'bzip2' or 'gzip'. Additional option for gzip compression
    level. Compressed filenames are only appended with suffix if not included.

    :param file_name: file name of output
    :param append: append to any existing file
    :param compress: gzip, bzip2
    :param gzlevel: gzip level (default 6)
    :return:
    """

    mode = 'w' if not append else 'w+'

    if compress == 'bzip2':
        if not file_name.endswith('.bz2'):
            file_name += '.bz2'
        # bz2 missing method to be wrapped by BufferedWriter. Just directly
        # supply a buffer size
        return bz2.BZ2File(file_name, mode, buffering=65536)
    elif compress == 'gzip':
        if not file_name.endswith('.gz'):
            file_name += '.gz'
        return io.BufferedWriter(gzip.GzipFile(file_name, mode, compresslevel=gzlevel))
    else:
        return io.BufferedWriter(io.FileIO(file_name, mode))


def multicopy_tostream(file_name, *ostreams, **kwargs):
    """
    Copy an input file to multiple output streams.
    :param file_name: input file name
    :param ostreams: output streams
    :param kwargs: optional parameters: write_mode (default 'w'), compress [gzip, bzip2] default: None
    :return:
    """
    bufsize = DEF_BUFFER if 'bufsize' not in kwargs else kwargs['bufsize']

    with open(file_name, 'r') as in_h:
        done = False
        while not done:
            buf = in_h.read(bufsize)
            if not buf:
                done = True
            for oi in ostreams:
                oi.write(buf)


def multicopy_tofile(file_name, *onames, **kwargs):
    """
    Copy an input file to multiple output files.
    :param file_name: input file name
    :param onames: output file names
    :param kwargs: optional parameters: write_mode (default 'w'), compress [gzip, bzip2] default: None
    :return:
    """
    bufsize = DEF_BUFFER if 'bufsize' not in kwargs else kwargs['bufsize']
    write_mode = "w" if 'write_mode' not in kwargs else kwargs['write_mode']
    compress = None if 'compress' not in kwargs else kwargs['compress']

    out_h = None
    try:
        in_h = open(file_name, 'r')
        out_h = [open_output(oi, write_mode, compress) for oi in onames]

        done = False
        while not done:
            buf = in_h.read(bufsize)
            if not buf:
                done = True
            for oi in out_h:
                oi.write(buf)
    finally:
        if out_h:
            for oi in out_h:
                if oi:
                    oi.close()


def write_to_stream(stream, data, fmt='plain'):
    """
    Write an object out to a stream, possibly using a serialization format
    different to default string representation.

    :param stream: open stream to twrite
    :param data: object instance
    :param fmt: plain, json or yaml
    """
    if fmt == 'yaml':
        yaml.dump(data, stream, default_flow_style=False)
    elif fmt == 'json':
        json.dump(data, stream, indent=1)
    elif fmt == 'plain':
        stream.write('{0}\n'.format(data))


def read_from_stream(stream, fmt='yaml'):
    """
    Load an object instance from a serialized format. How, in terms of classes
    the object is represented will depend on the serialized information. For
    generic serialized formats, this is more than likely to involve dictionaries
    for classes with properties.

    :param stream: open stream to read
    :param fmt: yaml or json
    :return: loaded object
    """
    if fmt == 'yaml':
        return yaml.safe_load(stream)
    elif fmt == 'json':
        return json_load_byteified(stream)


"""
Code below taken from Stack Exchange question.
http://stackoverflow.com/questions/956867/how-to-get-string-objects-instead-of-unicode-ones-from-json-in-python

JSON loading with UTF-8 encoding.

The following functions returns JSON results where Unicode strings are converted to UTF-8. Potentially an
unnecessary step as Python will handle referencing these strings transparently, but I wish to keep exchanged
data tables in a single encoding.

Attribution: Mirec Miskuf
"""


def json_loads_byteified(json_text):
    return _byteify(
        json.loads(json_text, object_hook=_byteify),
        ignore_dicts=True
    )


def json_load_byteified(file_handle):
    return _byteify(
        json.load(file_handle, object_hook=_byteify),
        ignore_dicts=True
    )


def _byteify(data, ignore_dicts=False):
    # if this is a unicode string, return its string representation
    if isinstance(data, unicode):
        return data.encode('utf-8')
    # if this is a list of values, return list of byteified values
    if isinstance(data, list):
        return [_byteify(item, ignore_dicts=True) for item in data]
    # if this is a dictionary, return dictionary of byteified keys and values
    # but only if we haven't already byteified it
    if isinstance(data, dict) and not ignore_dicts:
        return {
            _byteify(key, ignore_dicts=True): _byteify(value, ignore_dicts=True)
            for key, value in data.iteritems()
            }
    # if it's anything else, return it in its original form
    return data
