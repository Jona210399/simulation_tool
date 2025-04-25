class SimulationError(Exception):
    def __init__(
        self,
        simulation_type: str,
        return_value: int,
        message: str,
    ):
        error_msg = f"{simulation_type} simulation failed with return value {return_value} and message: {message}"

        self.simulation_type = simulation_type
        self.return_value = return_value
        self.message = message

        super().__init__(error_msg)
