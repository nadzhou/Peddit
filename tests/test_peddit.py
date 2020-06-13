from peddit import Peditor
from peddit import parse_arguments
import pytest



@pytest.fixture
def test_make_peddit_instance(): 
    inst = Peditor("4qdi")
    
    return inst


def test_edit_pdb_file(inst):
    
