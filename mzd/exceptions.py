
class ApplicationException(Exception):
    def __init__(self, message):
        super(ApplicationException, self).__init__(message)


class UnknownEnzymeException(ApplicationException):
    """All sequences were excluded during filtering"""
    def __init__(self, target, similar):
        super(UnknownEnzymeException, self).__init__(
            '{} is undefined, but its similar to: {}'.format(target, ', '.join(similar)))


class UnknownOrientationStateException(ApplicationException):
    """All sequences were excluded during filtering"""
    def __init__(self, ori):
        super(UnknownOrientationStateException, self).__init__('unknown orientation state [{}].'.format(ori))


class NoneAcceptedException(ApplicationException):
    """All sequences were excluded during filtering"""
    def __init__(self):
        super(NoneAcceptedException, self).__init__('all sequences were excluded')


class TooFewException(ApplicationException):
    """Method requires a minimum of nodes"""
    def __init__(self, minseq, method):
        super(TooFewException, self).__init__('More than {} sequences are required to apply {}'.format(minseq, method))


class NoRemainingClustersException(ApplicationException):
    def __init__(self, msg):
        super(NoRemainingClustersException, self).__init__(msg)


class NoReportException(ApplicationException):
    """Clustering does not contain a report"""
    def __init__(self, clid):
        super(NoReportException, self).__init__('No clusters contained an ordering'.format(clid))


class ZeroLengthException(ApplicationException):
    """Sequence of zero length"""
    def __init__(self, seq_name):
        super(ZeroLengthException, self).__init__('Sequence [{}] has zero length'.format(seq_name))


class ParsingError(ApplicationException):
    """An error during input parsing"""
    def __init__(self, msg):
        super(ParsingError, self).__init__(msg)
