VirusTaxo
----

## A tool for predicting the taxonomy of viral sequences

For taxonomic classification of viruses from metagenomic sequences, we developed VirusTaxo using diverse (e.g., 402 DNA and 280 RNA) genera of viruses. VirusTaxo has an average accuracy of 93% at genus level prediction in DNA and RNA viruses. VirusTaxo outperformed existing taxonomic classifiers of viruses where it assigned taxonomy of a larger fraction of metagenomic contigs compared to other methods. Benchmarking of VirusTaxo on a collection of SARS-CoV-2 sequencing libraries and metavirome datasets suggests that VirusTaxo can characterize virus taxonomy from highly diverse contigs and provide a reliable decision on the taxonomy of viruses.


## Quickstart
### Add VirusTaxo in LatchBio
1. Log into https://latch.bio.
2. Find **VirusTaxo** within **Explore** tab or, click [here](https://console.latch.bio/explore/63583/info).
3. Add VirusTaxo to your workspace (There's a button for it)
4. Go to **Workflows** tab and click on VirusTaxo.
### Provide Input and Run
1. Provide input `fasta` file. You can upload `fasta` file containing viral (single or multiple) sequences. Or, you can use test input data containing 1,553 SARS-CoV-2 genome sequences provided by us. [Download from here](https://mega.nz/file/JhAC0BRA#P1wQoYjj5mVscI-l8ADN_H723a_q2Jp4ISKpxPtGPwY).
2. Select Genome Type as **RNA** since SARS-CoV-2 is an RNA virus.
   <br/>
   *Note: If you do not know what the genome type is, select **Unknown** option. In case of DNA viruses, select **DNA**.*
3. Launch the workflow.


## Input Options
- `Input Fasta File`: Provide concatenated viral sequences in a single fasta file.
- `Virus Genome Type`: For DNA viruses, set this to `DNA`. For RNA viruses, set this to `RNA`. If the fasta file contains both DNA and RNA viral sequences, then selecting `Unknown` genome type will predict the taxonomy of both DNA and RNA viruses.


## How to cite
- Paper: https://doi.org/10.1016/j.ygeno.2022.110414


## Links:
- VirusTaxo Tool: https://console.latch.bio/explore/63583/info
- Website: https://omics-lab.com/virustaxo
- Paper Code: https://github.com/omics-lab/VirusTaxo
- Workflow Source Code: https://github.com/nahid18/virustaxo-wf


## FAQ
**Q:** Can I use `fastq` files?
<br/>
**A:** No, you cannot. But you can assemble the `fastq` reads using [Megahit](https://github.com/voutcn/megahit). Then use the Megahit contigs (in `fasta` format) in VirusTaxo.
