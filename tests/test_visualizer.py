"""
Tests for the visualizer module
"""

import pytest

from ppinsight import visualizer as vis

# setup variables
models = ['../dummy_scores/model_1.csv','/home/fmclary/CSE583/PPInsight/dummy_scores/model_2.csv','/home/fmclary/CSE583/PPInsight/dummy_scores/model_3.csv']

def test_not_a_list():
    """Edge test to ensure function throws TypeError if passed a non-list object"""
    with pytest.raises(
        TypeError, match='Format your input as a list of file paths'
    ):
        vis.to_plot('../dummy_scores/model_1.csv','../dummy_scores/model_2.csv')
    return
    
def test_score_type_does_not_exist():
    """Edge test to ensure function throws LookupError if passed a score type that does not exist in one or more input files"""
    with pytest.raises(
        LookupError, match=f'The score type {score_type} does not exist for the {file} model'
    ):
        frames = pd.DataFrame(np.array([1,2,3]),columns=['a','b','c'])
        vis.compare_scores(frames, 'd', 'test plot title')
    return

def test_smoke():
    """Basic smoke test to ensure function runs"""
    vis.to_plot(models)
    return

def test_output_is_df(dfs):
    """Pattern test to ensure output of to_plot is a list of pandas dataframes"""
    def setup(output):
        output.dfs = vis.to_plot(models)
    output.assertIsInstance(dfs is pd.DataFrame)
    return

def test_graph_dummy_models():
    """One-shot test to ensure compare_scores graphs score_a for the 3 dummy models as expected"""
    vis.compare_scores(vis.to_plot(models), 'score_a', 'test plot title')
    return
