import pytest
from pathlib import Path

DATA = Path(__file__).parent / "data"

@pytest.fixture(scope="session")
def tiny_bam():
    return DATA / "Cenp_A_Aligned.out.sorted.bam"

@pytest.fixture(scope="session")
def tiny_gtf():
    return DATA / "cenpa.gtf"

@pytest.fixture(scope="session")
def tiny_interval():
    # Known region in tiny.bam where you've placed reads
    return {"chr": "5", "strand": "+", "start": 30824121, "end": 30832174}


@pytest.fixture(scope="session")
def tiny_genesid():
	return "ENSMUSG00000029177:ENSMUST00000144742"
