class ATBError(BaseException):
    pass


class MissingTokenError(ATBError):
    pass


class InvalidTokenError(ATBError):
    pass


class ATBQuotaError(ATBError):
    pass
