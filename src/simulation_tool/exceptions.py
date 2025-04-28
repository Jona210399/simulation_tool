class SimulationError(Exception):
    def __init__(
        self,
        simulation_type: str,
        return_value: int,
        message: str,
    ):
        self.simulation_type = simulation_type
        self.return_value = return_value
        self.message = message.strip().replace("\n", " ")
        error_msg = (
            f"{simulation_type} sim returned: {self.return_value} | {self.message}"
        )

        super().__init__(error_msg)
