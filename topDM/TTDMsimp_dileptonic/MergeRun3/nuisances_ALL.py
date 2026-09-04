import copy
import os

# Run 3 merged nuisance scheme.
# Per-era nuisances keep an era suffix in the dictionary key whenever their
# Combine name is era-dependent, so mkMergeYears looks for the corresponding
# Up/Down histograms in the correct era input file. Common nuisance names are
# kept once and therefore remain correlated across the merged Run 3 periods.

_PERIODS = [
    ("2022", "../Full2022v12", "nuisances_ALL.py"),
    ("2022EE", "../Full2022EEv12", "nuisances_ALL.py"),
    ("2023", "../Full2023v12", "nuisances_ALL.py"),
    ("2023BPix", "../Full2023BPixv12", "nuisances_ALL.py"),
    ("2024", "../Full2024v15", "nuisances_ALL.py"),
]

_PERIOD_TOKENS = ("2022EE", "2022", "2023BPix", "2023", "2024")


def _merge_dir_candidates():
    candidates = []
    if "__file__" in globals():
        candidates.append(os.path.dirname(os.path.abspath(__file__)))
    candidates.extend([
        os.getcwd(),
        os.path.join(os.getcwd(), "topDM", "TTDMsimp_dileptonic", "MergeRun3"),
    ])
    return candidates


def _resolve_period_file(folder, nuisance_file):
    for base_dir in _merge_dir_candidates():
        nuisance_path = os.path.normpath(os.path.join(base_dir, folder, nuisance_file))
        if os.path.exists(nuisance_path):
            return nuisance_path
    return os.path.normpath(os.path.join(_merge_dir_candidates()[0], folder, nuisance_file))


def _load_period_nuisances(folder, nuisance_file):
    namespace = {
        "samples": samples,
        "cuts": cuts,
        "os": os,
        "copy": copy,
    }
    nuisance_path = _resolve_period_file(folder, nuisance_file)
    with open(nuisance_path) as handle:
        exec(compile(handle.read(), nuisance_path, "exec"), namespace)
    return namespace["nuisances"]


def _has_period_token(name):
    return any(token in name for token in _PERIOD_TOKENS)


def _renamed_key(key, period, nuisance):
    name = nuisance.get("name", key)
    if nuisance.get("type") == "auto":
        return key
    if _has_period_token(key) or _has_period_token(name):
        return key if key.endswith("_" + period) else f"{key}_{period}"
    return key


nuisances = {}
_seen_common_names = {}

for _period, _folder, _nuisance_file in _PERIODS:
    for _key, _nuisance in _load_period_nuisances(_folder, _nuisance_file).items():
        _entry = copy.deepcopy(_nuisance)
        _new_key = _renamed_key(_key, _period, _entry)

        if _entry.get("type") == "auto":
            nuisances.setdefault(_new_key, _entry)
            continue

        _name = _entry.get("name", _new_key)
        if _new_key == _key and not _has_period_token(_name):
            # Common lnN/shape nuisances are intentionally shared across all
            # periods. Keep the first full definition and merge in any samples
            # that only appear in later eras.
            _common_key = _seen_common_names.setdefault(_name, _new_key)
            if _common_key in nuisances and isinstance(nuisances[_common_key].get("samples"), dict):
                nuisances[_common_key]["samples"].update(_entry.get("samples", {}))
            else:
                nuisances[_common_key] = _entry
        else:
            nuisances[_new_key] = _entry
