# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

from io import StringIO
import unittest

from antismash.common import json
from antismash.common.json import JSONOrf
from antismash.common.test.helpers import DummyCDS


class DummyClass:
    def __init__(self, value):
        self.value = value


class TestConversions(unittest.TestCase):
    def test_orf_conversion(self):
        orf = JSONOrf.from_cds(DummyCDS(locus_tag="name"))
        assert orf.id == "name"
        for indent in [False, True]:  # catch the error of non-indenting in stdlib failing to convert
            new = JSONOrf(**json.loads(json.dumps(orf, indent=indent)))
            assert new == orf
            assert new.id == orf.id

    def test_class_without_conversion(self):
        with self.assertRaisesRegex(TypeError, "not JSON serializable"):
            json.dumps(DummyClass(5))

    def test_class_with_explicit(self):
        class WithExplicit(DummyClass):
            def to_json(self):
                return {"value": self.value}

        instance = WithExplicit(6)
        conversion = json.loads(json.dumps(instance))
        assert conversion["value"] == instance.value

    def test_class_with_dunder(self):
        class WithDunder(DummyClass):
            def __json__(self):
                return {"value": self.value}

        instance = WithDunder(7)
        conversion = json.loads(json.dumps(instance))
        assert conversion["value"] == instance.value

    def test_indents(self):
        instance = {"a": 2, "b": "val"}
        default = json.dumps(instance)
        assert default == '{"a":2,"b":"val"}'
        indented = json.dumps(instance, indent=True)
        assert indented == '{\n  "a": 2,\n  "b": "val"\n}'

    def test_sorted(self):
        instance = {"z": 2, "a": "val"}
        default = json.dumps(instance)
        assert default == '{"z":2,"a":"val"}'
        sorted_result = json.dumps(instance, sort_keys=True)
        assert sorted_result == '{"a":"val","z":2}'

    def test_key_conversion(self):
        # stdlib json conversion will convert non-string keys to strings
        # the behaviour of the override must match
        obj = {"str": 1, 7: 2}
        converted = json.loads(json.dumps(obj))
        assert converted == {"str": 1, "7": 2}

    def test_dump(self):
        handle = StringIO("")
        obj = {"a": 1, "b": 2}
        json.dump(obj, handle)
        handle.seek(0)
        assert handle.read() == json.dumps(obj)

    def test_load(self):
        obj = {"a": 1, "b": 2}
        handle = StringIO(json.dumps(obj))
        assert obj == json.load(handle)
