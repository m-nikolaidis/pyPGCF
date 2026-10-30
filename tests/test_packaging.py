import re
import tomllib
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _requirement_name(requirement: str) -> str:
    return re.split(r"[<>=!~\[ ]", requirement, maxsplit=1)[0].lower()


class TestPackagingConsistency(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pyproject = tomllib.loads((ROOT / "pyproject.toml").read_text())
        cls.pixi = tomllib.loads((ROOT / "pixi.toml").read_text())
        cls.recipe = (ROOT / "recipe.yaml").read_text()

    def test_pixi_version_matches_pyproject(self):
        self.assertEqual(
            self.pixi["workspace"]["version"], self.pyproject["project"]["version"]
        )

    def test_recipe_version_matches_pyproject(self):
        match = re.search(r'^  version: "([^"]+)"', self.recipe, re.MULTILINE)
        self.assertIsNotNone(match)
        self.assertEqual(match.group(1), self.pyproject["project"]["version"])

    def test_runtime_dependencies_declared_in_pixi(self):
        pixi_names = {name.lower() for name in self.pixi["dependencies"]}
        for requirement in self.pyproject["project"]["dependencies"]:
            with self.subTest(requirement=requirement):
                self.assertIn(_requirement_name(requirement), pixi_names)

    def test_runtime_dependencies_declared_in_recipe(self):
        recipe_names = {
            _requirement_name(match)
            for match in re.findall(r"^    - ([^\s].*)$", self.recipe, re.MULTILINE)
        }
        for requirement in self.pyproject["project"]["dependencies"]:
            with self.subTest(requirement=requirement):
                self.assertIn(_requirement_name(requirement), recipe_names)

    def test_python_range_matches_pyproject(self):
        self.assertEqual(
            self.pixi["dependencies"]["python"],
            self.pyproject["project"]["requires-python"],
        )


if __name__ == "__main__":
    unittest.main()
