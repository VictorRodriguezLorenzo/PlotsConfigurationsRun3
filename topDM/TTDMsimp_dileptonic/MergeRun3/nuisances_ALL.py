import copy
import math
import os
import sys

# =============================================================================
# Run-3 merged nuisance configuration
#
# Rules:
#
#   1. lnN identical in all campaigns
#        -> keep as one common lnN
#
#   2. lnN not identical in all campaigns
#        -> convert to shape
#
#      mkMergeYears.py will obtain the corresponding lnN from each original
#      campaign nuisances_ALL.py and build the merged Up/Down histograms.
#
#   3. Shape uncertainty with the same Combine name across campaigns
#        -> keep one common/correlated shape
#
#   4. Shape uncertainty whose Combine name changes with the campaign
#        -> create one separate nuisance per campaign
#
#   5. rateParam
#        -> one common Run-3 rateParam, not one per campaign
#
#   6. autoStats
#        -> one common definition
#
# The dictionary KEY is important for mkMergeYears.py:
#
#   - common nuisances retain their original key
#   - campaign-specific nuisances get the campaign in the key when necessary,
#     e.g. eff_e_2024
#
# The "name" remains the actual Combine nuisance name.
# =============================================================================


# -----------------------------------------------------------------------------
# Campaigns
# -----------------------------------------------------------------------------

_PERIODS = [
    ("2022",     "../Full2022v12",     "nuisances_ALL.py"),
    ("2022EE",   "../Full2022EEv12",   "nuisances_ALL.py"),
    ("2023",     "../Full2023v12",     "nuisances_ALL.py"),
    ("2023BPix", "../Full2023BPixv12", "nuisances_ALL.py"),
    ("2024",     "../Full2024v15",     "nuisances_ALL.py"),
#    ("2025",     "../Full2025v15",     "nuisances_ALL.py"),
]

# Longer tokens first!
_PERIOD_TOKENS = (
    "2023BPix",
    "2022EE",
#    "2025",
    "2024",
    "2023",
    "2022",
)

_ALL_PERIOD_NAMES = tuple(period for period, _, _ in _PERIODS)


# -----------------------------------------------------------------------------
# Path handling
# -----------------------------------------------------------------------------

def _merge_dir_candidates():
    candidates = []

    if "__file__" in globals():
        candidates.append(
            os.path.dirname(os.path.abspath(__file__))
        )

    candidates.append(os.getcwd())

    # Remove duplicates while preserving order.
    unique = []
    for path in candidates:
        path = os.path.abspath(path)
        if path not in unique:
            unique.append(path)

    return unique


def _resolve_period_file(folder, filename):
    for base_dir in _merge_dir_candidates():
        path = os.path.normpath(
            os.path.join(base_dir, folder, filename)
        )

        if os.path.exists(path):
            return path

    attempted = [
        os.path.normpath(os.path.join(base, folder, filename))
        for base in _merge_dir_candidates()
    ]

    raise FileNotFoundError(
        "Could not find campaign configuration file:\n"
        + "\n".join(f"  - {path}" for path in attempted)
    )


# -----------------------------------------------------------------------------
# Execute one campaign nuisances_ALL.py
# -----------------------------------------------------------------------------

def _load_period_nuisances(period, folder, nuisance_file):
    nuisance_path = _resolve_period_file(folder, nuisance_file)
    nuisance_dir = os.path.dirname(nuisance_path)

    # Start from the current mkShapesRDF configuration namespace.
    #
    # This gives the campaign nuisance file access to things such as:
    #   os
    #   ROOT
    #   lumi
    #   etc.
    #
    # We explicitly override samples/cuts with the merged configuration
    # objects. This is fine for constructing the merged nuisance scheme.
    namespace = globals().copy()

    namespace.update({
        "samples": samples,
        "cuts": cuts,
        "os": os,
        "copy": copy,
        "__file__": nuisance_path,
    })

    # Needed for local imports such as:
    #
    #   from theoryNormalizations_2024_Signals import ...
    #
    added_path = nuisance_dir not in sys.path

    if added_path:
        sys.path.insert(0, nuisance_dir)

    try:
        with open(nuisance_path, "r") as handle:
            code = compile(
                handle.read(),
                nuisance_path,
                "exec",
            )

            exec(code, namespace, namespace)

    except Exception as exc:
        raise RuntimeError(
            f"Error loading {period} nuisances from:\n"
            f"  {nuisance_path}\n"
            f"{type(exc).__name__}: {exc}"
        ) from exc

    finally:
        if added_path:
            sys.path.remove(nuisance_dir)

    period_nuisances = namespace.get("nuisances")

    if not isinstance(period_nuisances, dict):
        raise RuntimeError(
            f"{period}: 'nuisances' is not a dictionary in "
            f"{nuisance_path}. Got {type(period_nuisances).__name__}."
        )

    return period_nuisances


# -----------------------------------------------------------------------------
# Small helpers
# -----------------------------------------------------------------------------

def _entry_name(key, entry):
    return entry.get("name", key)


def _contains_period_token(text):
    text = str(text)
    return any(token in text for token in _PERIOD_TOKENS)


def _strip_period_suffix(text):
    """
    Examples:

        tt2lnorm_2024       -> tt2lnorm
        DYnorm_2023BPix     -> DYnorm
        eff_e_2022EE        -> eff_e

    Only strips a trailing campaign suffix.
    """
    text = str(text)

    for token in _PERIOD_TOKENS:
        suffix = "_" + token

        if text.endswith(suffix):
            return text[:-len(suffix)]

    return text


def _ensure_period_suffix(key, period):
    suffix = "_" + period

    if key.endswith(suffix):
        return key

    return key + suffix


def _canonical_lnn_value(value):
    """
    Return an lnN as a numerical (up, down) pair.

      1.05       -> (1.05, 1/1.05)
      1.05/0.95  -> (1.05, 0.95)

    This follows the convention used by mkMergeYears.py.
    """
    if value is None:
        return (1.0, 1.0)

    if isinstance(value, (int, float)):
        up = float(value)

        if up <= 0:
            raise ValueError(f"Invalid lnN value: {value}")

        return (up, 1.0 / up)

    value = str(value).strip()

    if "/" in value:
        up, down = value.split("/", 1)

        return (
            float(up),
            float(down),
        )

    up = float(value)

    if up <= 0:
        raise ValueError(f"Invalid lnN value: {value}")

    return (
        up,
        1.0 / up,
    )


def _same_lnn_value(a, b, rel_tol=1e-10, abs_tol=1e-12):
    a_up, a_down = _canonical_lnn_value(a)
    b_up, b_down = _canonical_lnn_value(b)

    return (
        math.isclose(
            a_up,
            b_up,
            rel_tol=rel_tol,
            abs_tol=abs_tol,
        )
        and
        math.isclose(
            a_down,
            b_down,
            rel_tol=rel_tol,
            abs_tol=abs_tol,
        )
    )


def _union_samples(records):
    """
    Merge sample dictionaries, keeping the first encountered value.

    For a merged shape nuisance the values themselves are not used by
    mkMergeYears to construct the combination: the original campaign
    nuisance definitions provide the actual variations.

    The sample keys, however, are important.
    """
    merged = {}

    for record in records:
        sample_dict = record["entry"].get("samples", {})

        if not isinstance(sample_dict, dict):
            continue

        for sample, value in sample_dict.items():
            if sample not in merged:
                merged[sample] = copy.deepcopy(value)

    return merged


def _shape_sample_membership(records):
    """
    Build a neutral shape-style sample dictionary.

    Used when an original lnN is converted into a merged shape nuisance.

    Only sample membership matters for the merged configuration. The actual
    Up/Down factors are recovered from each campaign's original nuisances.py.
    """
    samples_out = {}

    for record in records:
        sample_dict = record["entry"].get("samples", {})

        if not isinstance(sample_dict, dict):
            continue

        for sample in sample_dict:
            samples_out[sample] = ["1", "1"]

    return samples_out


def _union_list_field(records, field):
    result = []

    for record in records:
        values = record["entry"].get(field)

        if not isinstance(values, (list, tuple)):
            continue

        for value in values:
            if value not in result:
                result.append(value)

    return result


def _merge_common_metadata(entry, records):
    """
    Merge fields which may legitimately differ only because the campaign
    definitions cover different cuts/samples.
    """
    merged_samples = _union_samples(records)

    if merged_samples:
        entry["samples"] = merged_samples

    for field in ("cuts", "cutspost"):
        values = _union_list_field(records, field)

        if values:
            entry[field] = values

    return entry


# -----------------------------------------------------------------------------
# Decide whether an lnN can remain a single lnN
# -----------------------------------------------------------------------------

def _can_keep_common_lnn(records):
    """
    An lnN is kept as lnN only when:

      - it exists in every Run-3 campaign
      - every affected sample appears in every definition
      - the numerical factor is identical in every campaign

    Otherwise it is converted into a shape nuisance.

    This is intentionally conservative.
    """

    records_by_period = {
        record["period"]: record
        for record in records
    }

    # Must exist in every campaign.
    if set(records_by_period) != set(_ALL_PERIOD_NAMES):
        return False

    all_samples = set()

    for record in records:
        sample_dict = record["entry"].get("samples", {})

        if not isinstance(sample_dict, dict):
            return False

        all_samples.update(sample_dict.keys())

    for sample in all_samples:

        values = []

        for period in _ALL_PERIOD_NAMES:
            sample_dict = records_by_period[period]["entry"].get(
                "samples",
                {},
            )

            # If it is absent for one campaign, it is not the same lnN
            # across Run 3.
            if sample not in sample_dict:
                return False

            values.append(sample_dict[sample])

        reference = values[0]

        for value in values[1:]:
            try:
                if not _same_lnn_value(reference, value):
                    return False
            except Exception:
                return False

    return True


def _build_common_lnn(records):
    entry = copy.deepcopy(records[0]["entry"])

    all_samples = set()

    for record in records:
        all_samples.update(
            record["entry"].get("samples", {}).keys()
        )

    merged_samples = {}

    # Values have already been checked to be equivalent.
    for sample in all_samples:
        for record in records:
            sample_dict = record["entry"].get("samples", {})

            if sample in sample_dict:
                merged_samples[sample] = copy.deepcopy(
                    sample_dict[sample]
                )
                break

    entry["type"] = "lnN"
    entry["samples"] = merged_samples

    for field in ("cuts", "cutspost"):
        values = _union_list_field(records, field)

        if values:
            entry[field] = values

    return entry


# -----------------------------------------------------------------------------
# Determine the dictionary key to use for a merged shape nuisance
# -----------------------------------------------------------------------------

def _choose_shape_key(records, key_to_names):
    keys = {
        record["key"]
        for record in records
    }

    if len(keys) != 1:
        raise RuntimeError(
            "Cannot automatically correlate nuisance "
            f"'{_entry_name(records[0]['key'], records[0]['entry'])}'.\n"
            "The same Combine nuisance name uses different dictionary "
            f"keys across campaigns: {sorted(keys)}.\n"
            "This requires an explicit manual definition."
        )

    key = next(iter(keys))

    names_for_key = key_to_names.get(key, set())

    # Simple case:
    #
    # same key always corresponds to the same Combine nuisance.
    #
    if len(names_for_key) <= 1:
        return key

    # The same dictionary key is reused with different campaign-dependent
    # Combine names.
    #
    # Example:
    #
    #   2022: eff_e -> eff_e_2022
    #   2023: eff_e -> eff_e_2023
    #   2024: eff_e -> eff_e_2024
    #
    # Each group contains one campaign, so create:
    #
    #   eff_e_2022
    #   eff_e_2023
    #   eff_e_2024
    #
    periods = {
        record["period"]
        for record in records
    }

    if len(periods) == 1:
        period = next(iter(periods))
        return _ensure_period_suffix(key, period)

    # A partially-correlated nuisance whose original key is reused with
    # several different names requires more information than mkMergeYears
    # can infer from a generic key.
    raise RuntimeError(
        "Ambiguous partially-correlated nuisance:\n"
        f"  key      = {key}\n"
        f"  periods  = {sorted(periods)}\n"
        f"  names    = {sorted(names_for_key)}\n"
        "Define this nuisance manually in the merged nuisances_ALL.py."
    )


# -----------------------------------------------------------------------------
# rateParam handling
# -----------------------------------------------------------------------------

def _build_run3_rateparams(rate_records):
    """
    Campaign definitions such as

        key  = tt2lnorm
        name = tt2lnorm_2024

    become

        key  = tt2lnorm
        name = tt2lnorm

    in the merged Run-3 configuration.
    """

    grouped = {}

    for record in rate_records:
        key = record["key"]
        entry = record["entry"]

        base_key = _strip_period_suffix(key)
        base_name = _strip_period_suffix(
            entry.get("name", key)
        )

        group_id = base_name

        grouped.setdefault(
            group_id,
            {
                "base_key": base_key,
                "base_name": base_name,
                "records": [],
            },
        )

        grouped[group_id]["records"].append(record)

    output = {}

    for info in grouped.values():
        records = info["records"]

        entry = copy.deepcopy(records[0]["entry"])

        entry["type"] = "rateParam"
        entry["name"] = info["base_name"]

        # Merge process membership.
        merged_samples = _union_samples(records)

        if merged_samples:
            entry["samples"] = merged_samples

        # All control/signal regions covered by any campaign.
        cuts_merged = _union_list_field(records, "cuts")

        if cuts_merged:
            entry["cuts"] = cuts_merged

        # rateParam ranges should be consistent.
        ranges = [
            tuple(record["entry"]["range"])
            for record in records
            if "range" in record["entry"]
        ]

        if ranges:
            reference = ranges[0]

            if any(current != reference for current in ranges[1:]):
                raise RuntimeError(
                    f"Inconsistent rateParam ranges for "
                    f"{info['base_name']}: {ranges}"
                )

            entry["range"] = list(reference)

        output[info["base_key"]] = entry

    return output


# -----------------------------------------------------------------------------
# Load all campaign nuisance dictionaries
# -----------------------------------------------------------------------------

_PERIOD_NUISANCES = {}

for _period, _folder, _nuisance_file in _PERIODS:
    _PERIOD_NUISANCES[_period] = _load_period_nuisances(
        _period,
        _folder,
        _nuisance_file,
    )


# -----------------------------------------------------------------------------
# Convert everything into records
# -----------------------------------------------------------------------------

_records = []
_rate_records = []
_auto_records = []

for _period, _folder, _nuisance_file in _PERIODS:

    for _key, _entry_original in _PERIOD_NUISANCES[_period].items():

        _entry = copy.deepcopy(_entry_original)
        _type = _entry.get("type")

        _record = {
            "period": _period,
            "key": _key,
            "entry": _entry,
        }

        if _type == "rateParam":
            _rate_records.append(_record)
            continue

        if _type == "auto":
            _auto_records.append(_record)
            continue

        _records.append(_record)


# -----------------------------------------------------------------------------
# Track whether one dictionary key corresponds to different Combine names.
#
# Important for things like:
#
#   eff_e
#     2022 -> eff_e_2022
#     2023 -> eff_e_2023
#     2024 -> eff_e_2024
# -----------------------------------------------------------------------------

_key_to_names = {}

for _record in _records:
    _key = _record["key"]
    _name = _entry_name(
        _key,
        _record["entry"],
    )

    _key_to_names.setdefault(_key, set()).add(_name)


# -----------------------------------------------------------------------------
# Group nuisances by their actual Combine name.
#
# This is the main correlation criterion.
#
# Same "name"
#     -> same nuisance parameter / correlated
#
# Different "name"
#     -> different nuisance parameters
# -----------------------------------------------------------------------------

_groups_by_name = {}

for _record in _records:
    _name = _entry_name(
        _record["key"],
        _record["entry"],
    )

    _groups_by_name.setdefault(
        _name,
        [],
    ).append(_record)


# -----------------------------------------------------------------------------
# Build merged nuisance dictionary
# -----------------------------------------------------------------------------

nuisances = {}


for _name, _group in _groups_by_name.items():

    _types = {
        record["entry"].get("type")
        for record in _group
    }

    # -------------------------------------------------------------------------
    # lnN / shape
    #
    # A nuisance may even be lnN in one campaign and shape in another.
    # The merged representation is then necessarily a shape.
    # -------------------------------------------------------------------------

    if _types.issubset({"lnN", "shape"}):

        _all_lnn = _types == {"lnN"}

        # ---------------------------------------------------------------------
        # Identical common lnN
        # ---------------------------------------------------------------------

        if _all_lnn and _can_keep_common_lnn(_group):

            _merged_key = _choose_shape_key(
                _group,
                _key_to_names,
            )

            _merged_entry = _build_common_lnn(
                _group
            )

            nuisances[_merged_key] = _merged_entry

            continue

        # ---------------------------------------------------------------------
        # Everything else becomes / remains shape.
        #
        # Includes:
        #
        #   - campaign-specific shape
        #   - common shape
        #   - campaign-dependent lnN
        #   - same correlated lnN with different numerical values
        #   - lnN/shape mixture
        # ---------------------------------------------------------------------

        _merged_key = _choose_shape_key(
            _group,
            _key_to_names,
        )

        _merged_entry = copy.deepcopy(
            _group[0]["entry"]
        )

        _merged_entry["name"] = _name
        _merged_entry["type"] = "shape"

        # If any original entry was lnN, the merged nuisance is a synthetic
        # shape. Keep only neutral shape sample membership here; the true
        # variation is recovered by mkMergeYears from the original campaigns.
        if "lnN" in _types:
            _merged_entry["samples"] = _shape_sample_membership(
                _group
            )

            # These fields have no meaning for a synthetic lnN -> shape.
            _merged_entry.pop("kind", None)
            _merged_entry.pop("mapUp", None)
            _merged_entry.pop("mapDown", None)
            _merged_entry.pop("folderUp", None)
            _merged_entry.pop("folderDown", None)
            _merged_entry.pop("AsLnN", None)

        else:
            # Original shape nuisance: keep its normal metadata and merge the
            # sample membership across campaigns.
            _merged_entry = _merge_common_metadata(
                _merged_entry,
                _group,
            )

        # Preserve the union of cut restrictions.
        for _field in ("cuts", "cutspost"):
            _values = _union_list_field(
                _group,
                _field,
            )

            if _values:
                _merged_entry[_field] = _values

        nuisances[_merged_key] = _merged_entry

        continue


    # -------------------------------------------------------------------------
    # lnU
    #
    # mkMergeYears does not currently provide the same per-campaign conversion
    # machinery for lnU as it does for lnN. Do not silently mishandle one.
    # -------------------------------------------------------------------------

    if _types == {"lnU"}:
        raise RuntimeError(
            f"lnU nuisance '{_name}' found.\n"
            "Automatic campaign-dependent lnU merging is not implemented "
            "because mkMergeYears only reconstructs campaign lnN values."
        )


    # -------------------------------------------------------------------------
    # Any unexpected nuisance type
    # -------------------------------------------------------------------------

    raise RuntimeError(
        f"Unsupported nuisance type combination for '{_name}': "
        f"{sorted(str(x) for x in _types)}"
    )


# -----------------------------------------------------------------------------
# One Run-3 rateParam per normalization
# -----------------------------------------------------------------------------

_run3_rateparams = _build_run3_rateparams(
    _rate_records
)

for _key, _entry in _run3_rateparams.items():

    if _key in nuisances:
        raise RuntimeError(
            f"rateParam key collision: {_key}"
        )

    nuisances[_key] = _entry


# -----------------------------------------------------------------------------
# autoStats: one common definition
# -----------------------------------------------------------------------------

if _auto_records:

    _stat_key = _auto_records[0]["key"]
    _stat_entry = copy.deepcopy(
        _auto_records[0]["entry"]
    )

    # Check that the important autoStats settings are consistent.
    for _record in _auto_records[1:]:

        _other = _record["entry"]

        for _field in (
            "maxPoiss",
            "includeSignal",
        ):
            if (
                _field in _stat_entry
                and
                _field in _other
                and
                str(_stat_entry[_field]) != str(_other[_field])
            ):
                raise RuntimeError(
                    f"Inconsistent autoStats setting '{_field}': "
                    f"{_stat_entry[_field]} vs {_other[_field]}"
                )

    nuisances[_stat_key] = _stat_entry


# -----------------------------------------------------------------------------
# Final validation
# -----------------------------------------------------------------------------

# No two dictionary entries should accidentally generate the same Combine
# nuisance name after all merging decisions have been made.

_final_names = {}

for _key, _entry in nuisances.items():

    if _entry.get("type") == "auto":
        continue

    if "name" not in _entry:
        continue

    _name = _entry["name"]

    _final_names.setdefault(
        _name,
        [],
    ).append(_key)


for _name, _keys in _final_names.items():

    if len(_keys) > 1:
        raise RuntimeError(
            f"Duplicate final Combine nuisance name '{_name}' "
            f"from dictionary keys: {_keys}"
        )


# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------

_n_lnn = sum(
    entry.get("type") == "lnN"
    for entry in nuisances.values()
)

_n_shape = sum(
    entry.get("type") == "shape"
    for entry in nuisances.values()
)

_n_rate = sum(
    entry.get("type") == "rateParam"
    for entry in nuisances.values()
)

_n_auto = sum(
    entry.get("type") == "auto"
    for entry in nuisances.values()
)


print("")
print("==============================================================")
print(" Run-3 merged nuisance configuration")
print("==============================================================")
print(f" Total nuisances : {len(nuisances)}")
print(f"   lnN           : {_n_lnn}")
print(f"   shape         : {_n_shape}")
print(f"   rateParam     : {_n_rate}")
print(f"   auto          : {_n_auto}")
print("")

print("\n lnN nuisances:")
for _key, _entry in nuisances.items():
    if _entry.get("type") == "lnN":
        print(
            f"   {_key:35s} -> {_entry.get('name', _key)}"
        )

print("\n shape nuisances:")
for _key, _entry in nuisances.items():
    if _entry.get("type") == "shape":
        print(
            f"   {_key:40s} -> {_entry.get('name', _key)}"
        )

print("\n Run-3 rateParams:")
for _key, _entry in nuisances.items():
    if _entry.get("type") == "rateParam":
        print(
            f"   {_key:30s} -> {_entry.get('name', _key)}"
        )

print("==============================================================")
print("")
