from dataclass_wizard import JSONPyWizard

from simulation_tool.typing_ import PathLike


class JSONSerializable(JSONPyWizard):
    def to_json_file(self, json_file: PathLike) -> None:
        with open(json_file, "w") as file:
            file.write(self.to_json(indent=2))

    @classmethod
    def from_json_file(cls, json_file: PathLike) -> "JSONSerializable":
        with open(json_file, "r") as file:
            json_str = file.read()
        return cls.from_json(json_str)
