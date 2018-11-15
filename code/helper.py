from models import NamedVector


def apply_change_vector(state: NamedVector, change: NamedVector):
    ret = state.copy()
    for x in state:
        ret[x] = state[x] + change[x]
    return ret
