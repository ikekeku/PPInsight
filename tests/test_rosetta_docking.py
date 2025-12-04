import pytest

from ppinsight.docking import DockingPipeline



def test_smoke():
    pipeline = DockingPipeline("protein1.pdb", "protein2.pdb", n_runs=10)
    result = pipeline.run()
    print(f"Score: {result['final_score']}")
    return


# one-shot test: test againts known value
def test_oneshot():
    pipeline = DockingPipeline("protein1.pdb", "protein2.pdb", n_runs=10)
    result = pipeline.run()
    expected_score = -15.0  # hypothetical expected score
    assert abs(result['final_score'] - expected_score) < 1.0
    return

# edge test 
def test_edge_cases():
    # Test with zero runs
    with pytest.raises(ValueError, match="n_runs must be at least 1"):
        pipeline = DockingPipeline("protein1.pdb", "protein2.pdb", n_runs=0)
        pipeline.run()
    return

# pattern test
def test_pattern():
    for n_runs in [5, 10, 20]:
        pipeline = DockingPipeline("protein1.pdb", "protein2.pdb", n_runs=n_runs)
        complex_pose = pipeline.prepare()
        assert complex_pose is not None