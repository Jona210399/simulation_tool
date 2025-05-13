class SimulationError(Exception):
    def __init__(
        self,
        message: str,
    ):
        super().__init__(message)


class DeviceParametersIncompleteError(Exception):
    def __init__(
        self,
        message: str,
    ):
        super().__init__(message)
