"""共享归档追加维护的安全检查：保护当前状态并正确计算新增空间。"""

from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from extendMat5HistoryArchive import ROOT, object_delta, selected_result


class ExtensionSafety(unittest.TestCase):
    def test_current_states_and_traversal_are_rejected(self):
        for relative in (
            'outputs/v83_validation/krylov41/wall_velocity_krylov_evidence.mat',
            'outputs/retained_operator_results_20261030/v83_validation/krylov41/wall_velocity_krylov_evidence.mat',
            'outputs/retained_operator_results_20261030/../../../outside.mat',
            '../outside.mat',
            'outputs/retained_operator_results_20261030/v83_validation/krylov01/current_checkpoint.mat',
        ):
            with self.subTest(relative=relative), self.assertRaises((AssertionError, ValueError)):
                selected_result(relative)

    def test_only_two_existing_retired_groups_are_selected(self):
        for folder, case in (
            ('retained_operator_results_20261030', 'krylov01'),
            ('retained_graph_results_20261031', 'krylov17'),
        ):
            relative = 'outputs/%s/v83_validation/%s/wall_velocity_krylov_evidence.mat' % (folder, case)
            self.assertEqual(selected_result(relative), ROOT / relative)

    def test_existing_shared_bytes_are_not_counted_as_new(self):
        old = {'objects/old.z': {'bytes': 500, 'fileID': 1}}
        after = dict(old, **{'objects/new.z': {'bytes': 70, 'fileID': 2}})
        self.assertEqual(object_delta(old, after), ['objects/new.z'])
        self.assertEqual(object_delta(old, old), [])

    def test_missing_or_modified_shared_object_blocks_cleanup(self):
        old = {'objects/old.z': {'bytes': 500, 'fileID': 1}}
        for after in ({}, {'objects/old.z': {'bytes': 500, 'fileID': 8}},
                      {'objects/old.z': {'bytes': 499, 'fileID': 1}}):
            with self.subTest(after=after), self.assertRaises(AssertionError):
                object_delta(old, after)


if __name__ == '__main__':
    unittest.main()
