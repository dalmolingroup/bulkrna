# dalmolingroup/bulkrna pipeline parameters

Workflow for pre-processing, alignment and quantification of bulk RNA-Seq data

## Input/output options

Define where the pipeline should find input data and save output data.

| Parameter | Description | Type | Default |
|-----------|-----------|-----------|-----------|
| `input` | Path to comma-separated file containing information about the samples in the experiment. <details><summary>Help</summary><small>You will need to create a design file with information about the samples in your experiment before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row.</small></details>| `string` |  |
| `outdir` | The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure. | `string` |  |
| `email` | Email address for completion summary. <details><summary>Help</summary><small>Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.</small></details>| `string` |  |
| `multiqc_title` | MultiQC report title. Printed as page header, used for filename if not otherwise specified. | `string` |  |

## FastP trimming options

Define fastp parameters for trimming

| Parameter | Description | Type | Default |
|-----------|-----------|-----------|-----------|
| `reads_minlength` | reads shorter than length_required will be discarded | `integer` | 15 |
| `fastp_adapter_fasta` | specify a FASTA file to trim both read1 and read2 (if PE) by all the sequences in this FASTA file | `string` | https://gist.githubusercontent.com/jvfe/2d569adf47791b9d1e1b4ff810d410b8/raw/e970485d7fe2c9a5b501966eb26e9eb51789c86d/sample_adapters.fa |
| `fastp_save_trimmed_fail` | save reads that cannot pass the filters | `boolean` |  |
| `fastp_save_merged` | for paired-end input, merge each pair of reads into a single read if they are overlapped | `boolean` |  |
| `fastp_qualified_quality` | the quality value that a base is qualified | `integer` | 15 |
| `fastp_cut_mean_quality` | the mean quality requirement option shared by the cut_front and cut_tail sliding windows | `integer` | 15 |

## Kallisto options

Options for the alignment/quantification

| Parameter | Description | Type | Default |
|-----------|-----------|-----------|-----------|
| `transcriptome` | Path to the transcriptome FASTA file to use for reference | `string` | None |
| `index` |  | `string` | None |
| `fragment_length` | Estimated average fragment length | `integer` | 100 |
| `fragment_length_sd` | Estimated standard deviation of fragment length | `integer` | 1 |

## TXimport options

Options for importing Kallisto results into a count matrix

| Parameter | Description | Type | Default |
|-----------|-----------|-----------|-----------|
| `gtf` |  | `string` | None |

## Institutional config options

Parameters used to describe centralised config profiles. These should not be edited.

| Parameter | Description | Type | Default |
|-----------|-----------|-----------|-----------|
| `custom_config_version` | Git commit id for Institutional configs. | `string` | master |
| `custom_config_base` | Base directory for Institutional configs. <details><summary>Help</summary><small>If you're running offline, Nextflow will not be able to fetch the institutional config files from the internet. If you don't need them, then this is not a problem. If you do need them, you should download the files from the repo and tell Nextflow where to find them with this parameter.</small></details>| `string` | https://raw.githubusercontent.com/nf-core/configs/master |
| `config_profile_name` | Institutional config name. | `string` |  |
| `config_profile_description` | Institutional config description. | `string` |  |
| `config_profile_contact` | Institutional config contact information. | `string` |  |
| `config_profile_url` | Institutional config URL link. | `string` |  |

## Max job request options

Set the top limit for requested resources for any single job.

| Parameter | Description | Type | Default |
|-----------|-----------|-----------|-----------|
| `max_cpus` | Maximum number of CPUs that can be requested for any single job. <details><summary>Help</summary><small>Use to set an upper-limit for the CPU requirement for each process. Should be an integer e.g. `--max_cpus 1`</small></details>| `integer` | 16 |
| `max_memory` | Maximum amount of memory that can be requested for any single job. <details><summary>Help</summary><small>Use to set an upper-limit for the memory requirement for each process. Should be a string in the format integer-unit e.g. `--max_memory '8.GB'`</small></details>| `string` | 128.GB |
| `max_time` | Maximum amount of time that can be requested for any single job. <details><summary>Help</summary><small>Use to set an upper-limit for the time requirement for each process. Should be a string in the format integer-unit e.g. `--max_time '2.h'`</small></details>| `string` | 240.h |

## Generic options

Less common options for the pipeline, typically set in a config file.

| Parameter | Description | Type | Default |
|-----------|-----------|-----------|-----------|
| `help` | Display help text. | `boolean` |  |
| `publish_dir_mode` | Method used to save pipeline results to output directory. <details><summary>Help</summary><small>The Nextflow `publishDir` option specifies which intermediate files should be saved to the output directory. This option tells the pipeline what method should be used to move these files. See [Nextflow docs](https://www.nextflow.io/docs/latest/process.html#publishdir) for details.</small></details>| `string` | copy |
| `email_on_fail` | Email address for completion summary, only when pipeline fails. <details><summary>Help</summary><small>An email address to send a summary email to when the pipeline is completed - ONLY sent if the pipeline does not exit successfully.</small></details>| `string` |  |
| `plaintext_email` | Send plain-text email instead of HTML. | `boolean` |  |
| `max_multiqc_email_size` | File size limit when attaching MultiQC reports to summary emails. | `string` | 25.MB |
| `monochrome_logs` | Do not use coloured log outputs. | `boolean` |  |
| `hook_url` | Incoming hook URL for messaging service <details><summary>Help</summary><small>Incoming hook URL for messaging service. Currently, only MS Teams is supported.</small></details>| `string` |  |
| `multiqc_config` | Custom config file to supply to MultiQC. | `string` |  |
| `multiqc_logo` | Custom logo file to supply to MultiQC. File name must also be set in the MultiQC config file | `string` |  |
| `multiqc_methods_description` | Custom MultiQC yaml file containing HTML including a methods description. | `string` |  |
| `tracedir` | Directory to keep pipeline Nextflow logs and reports. | `string` | ${params.outdir}/pipeline_info |
| `validate_params` | Boolean whether to validate parameters against the schema at runtime | `boolean` | True |
| `show_hidden_params` | Show all params when using `--help` <details><summary>Help</summary><small>By default, parameters set as _hidden_ in the schema are not shown on the command line when a user runs with `--help`. Specifying this option will tell the pipeline to show all parameters.</small></details>| `boolean` |  |
