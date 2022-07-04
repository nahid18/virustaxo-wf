    VirusTaxo
    ----

    ## A tool for predicting the taxonomy of viral sequences

    For taxonomic classification of viruses from metagenomic sequences, we developed VirusTaxo using diverse (e.g., 402 DNA and 280 RNA) genera of viruses. VirusTaxo has an average accuracy of 93% at genus level prediction in DNA and RNA viruses. VirusTaxo outperformed existing taxonomic classifiers of viruses where it assigned taxonomy of a larger fraction of metagenomic contigs compared to other methods. Benchmarking of VirusTaxo on a collection of SARS-CoV-2 sequencing libraries and metavirome datasets suggests that VirusTaxo can characterize virus taxonomy from highly diverse contigs and provide a reliable decision on the taxonomy of viruses.

    ## Quickstart
    1. Download test input `fasta` file from [here](https://mega.nz/file/JhAC0BRA#P1wQoYjj5mVscI-l8ADN_H723a_q2Jp4ISKpxPtGPwY) containing 1,553 SARS-CoV-2 (RNA Virus) sequences.
    2. Create an input directory on [dashboard](https://console.latch.bio/data) and upload the `fasta` file.
    3. Add `VirusTaxo` to workspace from [explore](https://console.latch.bio/explore)
    4. Go to [workflows](https://console.latch.bio/workflows) and open `VirusTaxo` tool
    5. Select the input `fasta` file and select the genome type as `RNA`
    6. Launch the workflow


    ## Input Options
    - `Input Fasta File`: Provide concatenated viral sequences in a single fasta file.
    - `Virus Genome Type`: For DNA viruses, set this to `DNA`. For RNA viruses, set this to `RNA`. If `fasta_file` contains both DNA and RNA viral sequences, then selecting `Unknown` genome type will predict the taxonomy of both DNA and RNA viruses.


    ## How to cite
    - Paper: https://doi.org/10.1016/j.ygeno.2022.110414


    ## Links:
    - Contact: https://omics-lab.com/virustaxo
    - Source Code: https://github.com/omics-lab/VirusTaxo
    - Workflow Source Code: https://github.com/nahid18/virustaxo-wf