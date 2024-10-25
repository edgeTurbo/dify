from libs.exception import BaseHTTPException


class IllegalParametersError(BaseHTTPException):
    error_code = "illegal_parameters"
    description = "Illegal parameters."
    code = 500
