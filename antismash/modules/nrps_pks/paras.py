import os

from parasect.core.retrain_models import retrain_model, model_needs_retraining, update_metadata_file
from parasect.core.models import ModelType

from antismash.config import ConfigType
from antismash.common.path import find_latest_database_version


def _get_model_dir(config: ConfigType) -> str:
    """ A helper to construct the absolute path to the NRPS SVM model base dir in the
        data directory.

        Arguments:
            config: the antiSMASH config

        Returns:
            the absolute path of the NRPS SVM model base dir
    """
    root = os.path.join(config.database_dir, "nrps_pks", "paras")
    version = find_latest_database_version(root)
    return os.path.join(root, version)


def prepare_data(options: ConfigType, logging_only: bool = False) -> list[str]:
    """ Ensures packaged data is fully prepared

        Arguments:
            options: antiSMASH configuration options
            logging_only: whether to return error messages instead of raising exceptions

        Returns:
            a list of error messages (only if logging_only is True)
    """
    failure_messages: list[str] = []

    try:
        model_dir = _get_model_dir(options)
        metadata_path = os.path.join(model_dir, "metadata.txt")

        if not os.path.exists(metadata_path):
            if not logging_only:
                raise
            else:
                failure_messages.append(f"Failed to locate {metadata_path}")
        else:
            # TODO: Catch any exceptions appropriately
            retrain_paras_models_if_needed(metadata_path, model_dir)

    except ValueError:
        if not logging_only:
            raise
        else:
           failure_messages.append(f"Failed to locate PARAS model dir.")

    return failure_messages


def check_prereqs(options: ConfigType) -> list[str]:
    """ Check if the metadata file is present. """
    failure_messages: list[str] = []
    if "mounted_at_runtime" in options.database_dir:  # can't prepare this one
        return failure_messages
    try:
        model_dir = _get_model_dir(options)
    except ValueError:
        model_dir = ""
    if not os.path.exists(model_dir):
        failure_messages.append(f"Failed to locate {model_dir}")
    return failure_messages


def retrain_paras_models_if_needed(metadata_path: str, model_dir: str) -> None:
    """ Retrains PARAS models if sklearn versions do not match
    Args:
        metadata_path: path to PARAS model metadata file
        model_dir: path to PARAS model dir
    """
    for model_type in ModelType.ANTISMASH_MODELS:
        if model_needs_retraining(metadata_path, model_type):
            model = retrain_model(model_type)
            model.save(model_dir)
            update_metadata_file(model_type, metadata_path)
