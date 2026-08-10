import argparse

# 1. Define the allowed components (keep these synced with your run script)
VALID_BASES = {"testing", "prod", "tucker"}
VALID_MODIFIERS = {"short", "long", "rest2", "targetAngleRest2", "2fs"}


def validate_protocol(protocol_str: str) -> str:
    """
    Validates a protocol string like 'tucker_long_rest2'.
    Raises argparse.ArgumentTypeError if invalid, stopping execution immediately.
    """
    parts = protocol_str.split("_")
    base = parts[0]
    modifiers = parts[1:]

    # Check the base
    if base not in VALID_BASES:
        raise argparse.ArgumentTypeError(
            f"Invalid base protocol '{base}'. Allowed bases: {', '.join(VALID_BASES)}"
        )

    # Check the modifiers
    for mod in modifiers:
        if mod not in VALID_MODIFIERS:
            raise argparse.ArgumentTypeError(
                f"Invalid modifier '{mod}'. Allowed modifiers: {', '.join(VALID_MODIFIERS)}"
            )

    # Prevent duplicates like 'prod_long_long'
    if len(modifiers) != len(set(modifiers)):
        raise argparse.ArgumentTypeError(
            f"Duplicate modifiers found in '{protocol_str}'."
        )

    return protocol_str
