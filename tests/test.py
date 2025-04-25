from simulation_tool.templates.layer import Layer
from simulation_tool.templates.random import sample_random_parameters
from simulation_tool.templates.randomization import get_random_layer2

path_to_simss = "/path/to/simss"
path_to_session = "/path/to/session"
bandgap = 1.5
prefix = "l2"

myl = Layer.get_default_layer2(path_to_simss, bandgap).to_dict(prefix=prefix)
settings = sample_random_parameters(
    session_path=path_to_session,
    path_to_simss=path_to_simss,
    bandgap=bandgap,
)

l2_E_c = settings["l2.E_c"]
l2_E_v = settings["l2.E_v"]
myrl = get_random_layer2(path_to_session, bandgap).to_dict(prefix=prefix)


def get_keys_with_prefix(d: dict, prefix: str) -> dict:
    return {k: v for k, v in d.items() if k.startswith(prefix)}


otherl = get_keys_with_prefix(settings, f"{prefix}.")

print("myl", myl.keys())
print("otherl", otherl.keys())
print(myl.keys() == otherl.keys())

"""for (mdk, mdv), (ok, ov) in zip(myl.items(), otherl.items()):
    assert mdk == ok
    if mdv != ov:
        print(mdk, ":\tM: ", mdv, "\tO: ", ov)"""

for (mdk, mdv), (ok, ov), (mk, mv) in zip(myl.items(), otherl.items(), myrl.items()):
    assert mdk == ok
    if mdv != ov:
        print(mdk, ":\tD: ", mdv, "\tO: ", ov, "\tR: ", mv)
