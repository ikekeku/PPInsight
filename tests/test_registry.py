"""Tests for the engine plugin registry."""

import os
import tempfile
import pytest

from ppinsight import registry
from ppinsight.registry import EnginePlugin, register, get, list_engines, detect_engine


# ---------------------------------------------------------------------------
# Setup / teardown – ensure defaults are loaded before we add custom plugins
# ---------------------------------------------------------------------------

@pytest.fixture(autouse=True)
def _reset_registry():
    """Snapshot and restore the registry around each test."""
    orig = dict(registry._registry)
    orig_loaded = registry._defaults_loaded
    yield
    registry._registry.clear()
    registry._registry.update(orig)
    registry._defaults_loaded = orig_loaded


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestBuiltinEngines:
    """Built-in engines (lightdock, haddock, rosetta) are auto-registered."""

    def test_list_engines_includes_builtins(self):
        engines = list_engines()
        assert "lightdock" in engines
        assert "haddock" in engines
        assert "rosetta" in engines

    def test_get_lightdock_has_runner_and_parser(self):
        p = get("lightdock")
        assert p.runner is not None
        assert p.parser is not None
        assert p.detector is not None

    def test_get_haddock_has_runner_and_parser(self):
        p = get("haddock")
        assert p.runner is not None
        assert p.parser is not None

    def test_rosetta_is_parse_only(self):
        p = get("rosetta")
        assert p.runner is None
        assert p.parser is not None

    def test_unknown_engine_raises(self):
        with pytest.raises(KeyError):
            get("does_not_exist")


class TestRegisterCustom:
    """Users can register custom engines."""

    def test_register_and_retrieve(self):
        plugin = EnginePlugin(
            name="custom_dock",
            runner=lambda *a, **kw: "/fake/dir",
            parser=lambda *a, **kw: None,
            detector=lambda d: False,
            description="A custom test engine",
        )
        register(plugin)
        assert "custom_dock" in list_engines()
        assert get("custom_dock") is plugin

    def test_register_replaces_existing(self):
        original = get("lightdock")
        replacement = EnginePlugin(name="lightdock", description="replaced")
        register(replacement)
        assert get("lightdock").description == "replaced"

    def test_register_rejects_non_plugin(self):
        with pytest.raises(TypeError):
            register({"name": "bad"})  # type: ignore

    def test_custom_engine_in_list(self):
        register(EnginePlugin(name="z_engine"))
        names = list_engines()
        assert "z_engine" in names


class TestDetectEngine:
    """detect_engine delegates to the registry."""

    def test_detect_lightdock(self, tmp_path):
        (tmp_path / "swarm_0").mkdir()
        assert detect_engine(str(tmp_path)) == "lightdock"

    def test_detect_rosetta(self, tmp_path):
        (tmp_path / "docking_scores.csv").write_text("score\n1.0\n")
        assert detect_engine(str(tmp_path)) == "rosetta"

    def test_detect_haddock(self, tmp_path):
        sub = tmp_path / "9_caprieval"
        sub.mkdir()
        (sub / "capri_ss.tsv").write_text("score\tdockq\n1.0\t0.5\n")
        assert detect_engine(str(tmp_path)) == "haddock"

    def test_detect_unknown_raises(self, tmp_path):
        (tmp_path / "random_file.txt").write_text("hello")
        with pytest.raises(ValueError, match="Cannot detect"):
            detect_engine(str(tmp_path))

    def test_custom_detector_used(self, tmp_path):
        """A custom engine's detector is picked up by detect_engine."""
        marker = tmp_path / "custom_marker.json"
        marker.write_text("{}")

        plugin = EnginePlugin(
            name="aaa_custom",  # sorts first alphabetically
            detector=lambda d: os.path.isfile(os.path.join(d, "custom_marker.json")),
            description="test",
        )
        register(plugin)
        assert detect_engine(str(tmp_path)) == "aaa_custom"


class TestCollectScoresUsesRegistry:
    """collect_scores.detect_engine() delegates to registry."""

    def test_detect_matches_registry(self, tmp_path):
        (tmp_path / "swarm_0").mkdir()
        from ppinsight.collect_scores import detect_engine as cs_detect
        assert cs_detect(str(tmp_path)) == "lightdock"
