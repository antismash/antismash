# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Functions to load a JSON file and validate it against a JSON schema """

import json
import os
from typing import Any, Dict

import jsonschema

from antismash.common.errors import AntismashInputError


def _default_validator(validator_class: Any) -> Any:
    """ Extend a validator to insert defaults for missing properties before validation """

    # this is jsonschema's preferred method of extending validator classes
    validate_properties = validator_class.VALIDATORS["properties"]

    def add_defaults(validator: Any, properties: Dict[str, Any], instance: Any, schema: str) -> Any:
        """ Callback function for when json schema is validating a properties item """
        for property_name, subschema in properties.items():
            if "default" in subschema:
                instance.setdefault(property_name, subschema["default"])

        for error in validate_properties(validator, properties, instance, schema):
            yield error

    # add the callback
    return jsonschema.validators.extend(validator_class, {"properties": add_defaults})


# jsonschema's validators don't insert the default value, by default
# so build a draft 7 validator that extends with defaults
DefaultValidator = _default_validator(jsonschema.Draft7Validator)


def _ensure_valid(raw_json: Dict[str, Any], schema_file: str) -> Dict[str, Any]:
    """ Validate a python representation of JSON against a schema file """
    with open(schema_file, encoding="utf-8") as handle:
        schema = json.load(handle)

    # work around jsonschema breaking relative locations when subschemas reference each other
    schema_dir = os.path.dirname(schema_file)
    resolver = jsonschema.RefResolver(base_uri=f"file://{schema_dir}/subschemas/", referrer=schema)

    try:
        DefaultValidator(schema, resolver=resolver).validate(raw_json)
    except jsonschema.ValidationError as err:
        brief_message = str(err).splitlines()[0]
        path = "][".join(map(str, list(err.path)))
        if path:
            raise AntismashInputError(f"invalid sideload annotation for '[{path}]': {brief_message}")
        raise AntismashInputError(f"invalid sideload annotations: {brief_message}")
    return raw_json


def load_validated_json(data_file: str, schema_file: str) -> Dict[str, Any]:
    """ Convert the given JSON file to a python dict and validate it against the
        given schema (as a draft 7 JSON schema).
        Raises a ValueError if the data cannot be validated against the schema.

        Any separate subschemas referred to by the schema file are assumed to be
        local and in a subdirectory named "subschemas".

        Arguments:
            annotations_file: the JSON file containing data
            schema_file: the JSON file containing the schema

        Returns:
            A dictionary with the python representation of the JSON data

    """
    with open(data_file, encoding="utf-8") as handle:
        try:
            raw = json.load(handle)
        except ValueError as err:
            raise AntismashInputError(f"sideloaded data is not valid JSON: {err}") from err
        return _ensure_valid(raw, schema_file)
