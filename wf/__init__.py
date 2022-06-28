"""
Predict viral taxonomy
"""

import subprocess
from pathlib import Path

from latch import large_task, workflow
from latch.types import LatchFile, LatchDir
from enum import Enum


class GenomeType(Enum):
    dna = "DNA"
    rna = "RNA"
    unknown = "Unknown"


@large_task
def batch_prediction_task(
    input_dir: LatchDir,
    genome_type: GenomeType,
) -> LatchDir:
    """
    Run VirusTaxo on the input directory.
    """

    log_file = Path(f"/root/virustaxo_log.txt")
    output_dir = Path(f"/root/VirusTaxoRun/")
    # allowed_files = ['.fasta', '.fa', '.fastq', '.fq', '.FASTA', '.FA', '.FASTQ', '.FQ']
    # input_files = [f for f in Path(input_dir).iterdir() if f.suffix in allowed_files]
    # file_paths_as_string = [f.as_posix() for f in input_files]
    # if Path(genome_type).suffix == '.conf':
    #     _assembly_cmd = [
    #         "./shasta", 
    #         "--input", 
    #         str(" ".join(file_paths_as_string)),
    #         "--genome_type", 
    #         str(Path(genome_type).resolve()),
    #         "--assemblyDirectory",
    #         str(output_dir),
    #     ]
    #     with open(log_file, "w") as f:
    #         subprocess.run(_assembly_cmd, stdout=f, stderr=f)
    # else:
    #     raise ValueError(f"{genome_type} is not a valid genome_type file.")
    return LatchDir(str(output_dir), f"latch://{output_dir}")


@workflow
def virustaxo(
    input_dir: LatchDir, 
    genome_type: GenomeType,
) -> LatchDir:
    """Prediction of viral taxonomy

    VirusTaxo
    ----

    # VirusTaxo: A tool for predicting the taxonomy of viral sequences

    __metadata__:
        display_name: VirusTaxo
        author:
            name: Abdullah Al Nahid
            email: abdnahid56@gmail.com
            github: https://github.com/nahid18
        repository:
        license:
            id: MIT

    Args:
        input_dir:
          Input directory containing FASTA files to be predicted
          __metadata__:
            display_name: Input Directory
        genome_type:
          Genome type to use for prediction
          __metadata__:
            display_name: Genome Type
    """
    return batch_prediction_task(input_dir=input_dir, genome_type=genome_type)
