"""
Predict viral taxonomy
"""

import subprocess
from pathlib import Path

from latch.types import LatchFile, LatchDir
from latch import large_gpu_task, workflow
from collections import defaultdict
from datetime import datetime
from enum import Enum
from Bio import SeqIO
import numpy as np
import pickle


class GenomeType(Enum):
    dna = "DNA"
    rna = "RNA"
    unknown = "Unknown"


@large_gpu_task
def prediction_task(
    fasta_file: LatchFile,
    model_file: LatchFile,
    k: int,
) -> (LatchFile, LatchFile):
    """
    Run VirusTaxo on the input directory.
    """

    def entropy(x):
        eps = 1e-11
        if len(x) == 1:
            return 0.0
        shift = max(x)
        x -= shift
        summation = sum(np.exp(x))
        p = np.exp(x) / summation
        entropy = sum(p * np.log(p+eps)) / np.log(len(p))
        entropy *= -1
        return entropy


    allowed_files = ['.fasta', '.fa', '.FASTA', '.FA']

    if not Path(fasta_file).suffix in allowed_files:
        raise ValueError(f"{fasta_file} is not a fasta file")

    dt = datetime.now().strftime("%Y%m%d-%H%M%S")
    log_file = Path(f"/root/VirusTaxo_log_{dt}.txt")
    outfile = Path(f"/root/VirusTaxo_Result_{dt}.tsv")

    with open(log_file, "w") as fl:
        fl.write(f"Input File: {fasta_file}\n")
        fl.write(f"Model File: {model_file}\n")
        fl.write(f"K: {k}\n")
        fl.write(f"\n")
        fl.write(f"Output File: {outfile}\n")
        fl.write(f"Log File: {log_file}\n")
        fl.write(f"\n")
        fl.write(f"Downloading {model_file}\n")

    with open(model_file, 'rb') as fm:
        model = pickle.load(fm)

    with open(outfile, "w") as fh:
        fh.write('Id\tLength\tGenus\tEntropy\n')

        for record in SeqIO.parse(fasta_file, "fasta"):
            def get_genus(seq):
                count = defaultdict(int)
                for idx in range(len(seq) - k + 1):
                    kmer = seq[idx:idx + k]
                    if kmer in model:
                        matched_genera = model[kmer]
                        for genus in matched_genera:
                            count[genus] += 1
                genus_score_tuple = list(count.items())

                if not genus_score_tuple:
                    E = 1.0
                    return -1, -1, E
                else:
                    E = entropy(np.array([x for _, x in genus_score_tuple]))
                    prediction = max(genus_score_tuple, key=lambda x: x[1])[0]
                    cnt = max(genus_score_tuple, key=lambda x: x[1])[1]
                    return prediction, cnt, E

            prediction_1, cnt_1, E_1 = get_genus(record.seq)
            prediction_2, cnt_2, E_2 = get_genus(record.seq.reverse_complement())
            length = len(record.seq)

            if prediction_1 == -1 and prediction_2 == -1:
                fh.write(f'{record.id}\t{length}\tUnclassified\t1.0\n')
            else:
                if cnt_1 > cnt_2:
                    fh.write(f'{record.id}\t{length}\t{prediction_1}\t{E_1}\n')
                else:
                    fh.write(f'{record.id}\t{length}\t{prediction_2}\t{E_2}\n')

    return (
        LatchFile(str(log_file), f"latch://{log_file}"),
        LatchFile(str(outfile), f"latch://{outfile}"),
    )


@workflow
def virustaxo(
    fasta_file: LatchFile, 
    genome_type: GenomeType = GenomeType.rna,
) -> (LatchFile, LatchFile):
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
        fasta_file:
          Input fasta file containing virus sequences
          __metadata__:
            display_name: Input Fasta File
        genome_type:
          Select the virus genome type
          __metadata__:
            display_name: Virus Genome Type
    """

    if genome_type is GenomeType.dna:
        model_file = LatchFile("s3://latch-public/virustaxo/virustaxo-dna.pkl")
        k = 21
    elif genome_type is GenomeType.rna:
        model_file = LatchFile("s3://latch-public/virustaxo/virustaxo-rna.pkl")
        k = 17
    else:
        model_file = LatchFile("s3://latch-public/virustaxo/virustaxo-all.pkl")
        k = 20

    return prediction_task(
        fasta_file=fasta_file,
        model_file=model_file,
        k=k,
    )
