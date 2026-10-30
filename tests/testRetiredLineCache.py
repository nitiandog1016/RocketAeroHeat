"""核对共享历史对象的边界：必须覆盖全部链接，且当前Q41不在清理范围。"""

from pathlib import Path
import os
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import compactRetiredLineOperators as cleanup


class RetiredLineCacheTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.output = self.root / 'outputs'; self.output.mkdir()
        self.original_root, self.original_out, self.original_safe = cleanup.ROOT, cleanup.OUT, cleanup.safe
        cleanup.ROOT, cleanup.OUT, cleanup.safe = self.root, self.output, lambda p: Path(p)
        self.a = self.output / 'a.mat'; self.a.write_bytes(b'unchanged-science')
        self.b = self.output / 'b.mat'; os.link(str(self.a), str(self.b))
        stat = self.a.stat()
        self.inventory = {p.relative_to(self.root).as_posix(): dict(device=stat.st_dev, fileID=stat.st_ino)
                          for p in (self.a, self.b)}

    def tearDown(self):
        cleanup.ROOT, cleanup.OUT, cleanup.safe = self.original_root, self.original_out, self.original_safe
        self.temporary.cleanup()

    def test_complete_alias_set_accepted_once(self):
        self.assertEqual(cleanup.aliases_from_inventory(self.a, self.inventory), ['outputs/a.mat', 'outputs/b.mat'])
        self.assertEqual(self.a.read_bytes(), b'unchanged-science')

    def test_unknown_link_is_refused(self):
        os.link(str(self.a), str(self.output / 'unlisted.mat'))
        with self.assertRaisesRegex(AssertionError, '硬链接'):
            cleanup.aliases_from_inventory(self.a, self.inventory)

    def test_alias_outside_outputs_is_refused(self):
        outside = self.root / 'outside.mat'; os.link(str(self.a), str(outside))
        stat = outside.stat()
        self.inventory['outside.mat'] = dict(device=stat.st_dev, fileID=stat.st_ino)
        with self.assertRaises(ValueError):
            cleanup.aliases_from_inventory(self.a, self.inventory)

    def test_changed_alias_is_refused(self):
        self.b.unlink(); self.b.write_bytes(b'another-object')
        with self.assertRaises(AssertionError):
            cleanup.aliases_from_inventory(self.a, self.inventory)

    def test_current_state_and_unknown_process_excluded(self):
        self.assertEqual(len(cleanup.TARGETS), 20)
        names = {name for name, _ in cleanup.TARGETS}
        self.assertFalse(any('krylov41' in name or 'krylov25' in name for name in names))
        self.assertFalse(any('v101' in name or 'v102' in name or 'v103' in name for name in names))

    def test_explicit_old_engine_markers_and_missing_marker(self):
        self.assertEqual(cleanup.process_engine_ids({'engineProcessIDsFromBatchMarker': [18772]}), [18772])
        self.assertEqual(cleanup.process_engine_ids({'engineProcessIDs': [48732]}), [48732])
        with self.assertRaises(AssertionError):
            cleanup.process_engine_ids({'processID': 18772})

    def test_two_exact_source_tables_and_missing_sha(self):
        expected = {'builder.m': 'a'*64}
        self.assertEqual(cleanup.process_sources({'sourcesBefore': expected}), expected)
        self.assertEqual(cleanup.process_sources({'sourceSHA256': expected}), expected)
        with self.assertRaises(AssertionError):
            cleanup.process_sources({'sourceStableAfter': {'builder.m': True}})

    def test_constructor_dependencies_preserved_without_unrelated_snapshot(self):
        data = {'buildRocketShellLineOperator': b'function x=buildRocketShellLineOperator; x=block(); end',
                'block': b'function x=block; x=primitive(); end',
                'primitive': b'function x=primitive; x=1; end',
                'setupCutCellFlow': b'function setupCutCellFlow; end',
                'unusedResearch': b'function unusedResearch; end'}
        table = {'internalFlow/'+name+'.m': 'a'*64 for name in data}
        visited = []
        def load(relative, digest):
            visited.append(relative)
            return data[Path(relative).stem]
        result = cleanup.construction_closure(table, load)
        self.assertEqual({Path(p).stem for p in visited}, set(data)-{'unusedResearch'})
        self.assertEqual(len(result['bindings']), 4)

    def test_dynamic_constructor_dependency_is_refused(self):
        table = {name+'.m': 'a'*64 for name in ('buildRocketShellLineOperator','setupCutCellFlow')}
        with self.assertRaisesRegex(AssertionError, '动态依赖'):
            cleanup.construction_closure(table, lambda p,h: b'feval(options.method)')


if __name__ == '__main__':
    unittest.main(verbosity=2)
