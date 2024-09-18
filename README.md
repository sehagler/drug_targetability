Drug Targetability
Problem Being Addressed
The software is designed to address the challenge of evaluating drug treatments within the complex context of biological networks and pathways. It recognizes that drugs can be "promiscuous," targeting multiple genes or pathways, not just the ones of interest. This complexity makes it difficult to assess the full impact of a drug solely based on individual knowledgebases. By providing metrics to evaluate drug target specificity and navigate the hierarchical nature of biological networks, the software aids in identifying candidate treatment options to increase specificity and minimize off-target effects. This software is research use only currently.

Citations
ADD CITATIONS ASSOCIATED WITH THIS PROJECT HERE

Package Overview
This software addresses the problem by providing the clinician with a set of metrics that characterize the overall effect of a drug or collection of drugs on a patient. These metrics can be used by the clinician to compare different treatment options to arrive at the treatment that should have the fewest or least serious side-effects.

DATABASES
The software performs calculations using information available in three public databases:

Cancer Targetome
The Cancer Targetome aggregates drug-target interaction and bioactivity information for FDA-approved antineoplastic drugs across four publicly available resources. It offers a novel contribution due to both the inclusion of putative target interactions encompassing multiple targets for each antineoplastic drug and the introduction of a framework for categorizing the supporting evidence behind each drug-target interaction.

Reactome Knowledgebase (Reactome)
The Reactome Knowledgebase is an open-source, open access, manually curated and peer-reviewed pathway database. Its goal is to provide intuitive bioinformatics tools for the visualization, interpretation and analysis of pathway knowledge to support basic and clinical research, genome analysis, modeling, systems biology and education.

Universal Protein Resource (Uniprot)
The Universal Protein Resource is a comprehensive resource for protein sequence and annotation data. The UniProt databases are the UniProt Knowledgebase (UniProtKB), the UniProt Reference Clusters (UniRef), and the UniProt Archive (UniParc). The UniProt consortium and host institutions EMBL-EBI, SIB and PIR are committed to the long-term preservation of the UniProt databases.

Licenses
The source code to the is made available under the GNU General Public License 3.0

The license information for the Cancer Targetome is available at The-Cancer-Targetome Repository

The license information for the Reactome Knowledgebase is available at Reactome License Agreement

The license information for the Universal Protein Resource is available at Uniprot License & Disclaimer

Installation
Download the <targetability_project_name> project and copy to the desired location.

Install the Reactome flat files:

2.1. Go to Reactome Downloads

2.2. Go to Identifier mapping files > All levels of the pathway hierarchy > UniProt to All pathways to get the file UniProt2Reactome_All_Levels.txt by right-clicking and selecting ‘Save link as ...’ and save the file into the reactome_flat_files directory

2.3. Go to Pathways > Complete List of Pathways to get the file ReactomePathways.txt by right-clicking and selecting ‘Save link as ...’ and save the file into the reactome_flat_files directory

2.4. Go to Pathways > Pathways hierarchy relationship to get the file ReactomePathwaysRelation.txt by right-clicking and selecting ‘Save link as ...' and save the file into the reactome_flat_files directory

Install the Targetome flat files:

3.1. Go to The-Cancer-Targetome/results_070617/

3.2. Go to Targetome_Evidence_TIPS_101017.txt and select ‘Download raw file’, note you must select Targetome_Evidence_TIPS_101017.txt first and this will take you to the screen that allows you to download. Save the file into the targetome_flat_files directory

Initialize the <targetability_project_name> PKL-file:

4.1. Open a console to the <targetability_project_name> directory and run the command:

targetability --update

The Targetome and Reactom flat files should be regularly updated as iterative improvements become available. Any time a flat file is updated, the targetability --update command must be run again.

PAPER
A more detailed treatment of the construction and use of the metrics calculated by this software is provided in . Only a high-level summary is provided in the following.

ADD LINK TO PAPER HERE

REACTOME TREES
The Reactome Knowledgebase consists of units called pathways which are organized hierarchically in binary parent-child relationships forming a lattice where a pathway is a grouping of interlinked reactions where a reaction is a process that convertes input molecules and/or complexes into output molecules and/or complexes. This lattice may be decomposed into a number of trees: Autophagy, Cell Cycle, Cell-Cell communication, Cellular responses to stimuli, Chromatin organization, Circadian Clock, DNA Repair, DNA Replication, Developmental Biology, Digestion and absorption, Disease, Drug ADME, Extracellular matrix organization, Gene expression (Transcription), Hemostasis, Immune System, Metabolism, Metabolism of RNA, Metabolism of proteins, Muscle contraction, Neuronal System, Organelle biogenesis and maintenance,Programmed Cell Death, Protein localization, Reproduction, Sensory Perception, Signal Transduction, Transport of small molecules, and Vesicle-mediated transport.

TREE AND COVERAGE METRICS
Each vertex (pathway) of a tree involves a number of targets. Subtrees of a tree can be constructed by selecting a collection of targets and constructing the subtree containing only vertices (pathways) involving targets in that collection. Moreover, using the Cancer Targetome, a collection of drugs can be associated with the collection of targets targeted by the drugs up to some max assay value. Performing this construction on all the Reactome trees, a collection of drugs can be characterized by appropriate subtrees of the Reactome trees.

The software facilitates this kind of reasoning using subtrees of the Reactome trees by calculating values for the following metrics:

Tree metrics:
P The number of unique targets involved in the Reactome pathways contained in the tree. |L| The number of lineages making up the tree where a lineage is a maximal subtree whose root vertex is a child of the root vertex of the original tree. |Pi| The number of paths making up the tree where a path is a sequence of vertices from the root vertex of the tree to a leaf vertex of the tree where each vertex (except for the root vertex) is a child of the previous vertex. eta A summary measure of the structure of the tree.

Coverage metrics:
KP The fraction of unique targets in the tree that are targeted by the collection of drugs K|L| The fraction of lineages in the tree that are targeted by the collection of drugs k|L| The ratio of the number of lineages targed by the collection of drugs to the number of lineages involving at least one member of a collection of targets. K|Pi| The fraction of paths in the tree that are targeted by the collection of drugs k|Pi| The ratio of the number of paths targed by the collection of drugs to the number of paths involving at least one member of a collection of targets.

Use Cases
Targetability --update
Example: targetability --update

The update command generates a new PKL-file from the flat data files. The update command should be run first when the software is install, and run again whenever any of the flat data files are updated. Except when running the update command the software does not read the flat data files directly, but only reads what is in the PKL-file.

Targetability --drug_coverage (--max_assay_value) (--filename)
Example: targetability --drug_coverage "Imatinib Mesylate" Vemurafenib NODRUG

The drug_coverage command with only a list of drugs generates the drug coverage metric data for the drug combination Imatinib Mesylate, Vemurafenib, and NODRUG. Note that Imatinib Mesylate must be entered using quotation marks to be read as a the name of a single drug, otherwise the code will interpret Imatinib and Mesylate as being the names of two distinct drugs. NODRUG has been included to illustrate how the software handles a drug name that does not appear in the Cancer Targetome. No values have been specified for the optional command line arguments --max_assay_value and --filename. This makes the max assay value take on the default value of 100, and the output files are written using the default filenames drug_coverage_YYYYmmdd_HHMMSS.csv and drug_coverage_YYYYmmdd_HHMMSS.txt.

Running the command generates the following console output:

--target_coverage not specified --max_assay_value not specified --filename not specified Targetome query failed to find drug name NODRUG Fetched: 13 / 13

The output indicates that the line arguments --target_coverage, --max_assay_value, and --filename have not been specified. The absence of a value for the --target_coverage argument means that the code will not generate values for the metrics k|L| and k|Pi| as these values are determined for specific sets of targets. The absence of values for the --max_assay_value and --filename line arguments means that the default values for these will be used. The output also indicates that the drug name NODRUG could not be found in the Cancer Targetome. Finally the output indicates that a query to Uniprot fetched 13 items.

The two files drug_coverage_YYYYmmdd_HHMMSS.csv and drug_coverage_YYYYmmdd_HHMMSS.txt should now appear in the directory from which the command was run. The two files contain the same information. In drug_coverage_YYYYmmdd_HHMMSS.txt the data are formatted to make them human readable, while drug_coverage_YYYYmmdd_HHMMSS.csv presents the data in a comma-separated value format so that they may be easily read into another software application.

The generated drug_coverage_YYYYmmdd_HHMMSS.txt contains the following information:

example/drug_coverage_20240502_134548.csv

example/drug_coverage_20240502_134548.txt

The data are interpreted as follows:

INPUT_DRUGS lists the drugs input using line arguments.

INPUT_MAX_ASSAY_VALUE gives the max assay value being used. In this case it is the default value of 100.

INPUT_DRUGS_IN_TARGETOME lists the drugs in INPUT_DRUGS that were found to appear in the Cancer Targetome. In this case since NODRUG does not appear in the Cancer Targetome it does not appear here.

TARGETS_OF_INPUT_DRUGS_IN_TARGETOME lists all targets of the collection of drugs in INPUT_DRUGS_IN_TARGETOME with an assay value indictaed in the Cancer Targetome to be less than or equal to INPUT_MAX_ASSAY_VALUE for at least one drug in INPUT_DRUGS_IN_TARGETOME.

TARGETS_OF_INPUT_DRUGS_IN_TARGETOME_UNIPROT_ACCESSION_IDS lists the Uniprot Accessions Ids for all targets in TARGETS_OF_INPUT_DRUGS_IN_TARGETOME. The entries in TARGETS_OF_INPUT_DRUGS_IN_TARGETOME and TARGETS_OF_INPUT_DRUGS_IN_TARGETOME_UNIPROT_ACCESSION_IDS should appear in the same order.

DRUG_COVERAGE provides a table of the values of the metrics P, |L|, |Pi|, eta and KP, K|L|, K|Pi|, eta measuring the coverage of the collection of drugs in INPUT_DRUGS_IN_TARGETOME across the trees extracted from Reactome. No values for k|L| or k|Pi| are provided as no list of targets was input by the user. The first four metrics P, |L|, |Pi|, eta provide summary information about the structures of the trees. The next four metrics KP, K|L|, K|Pi|, eta provide summary information on the subtree of each tree consisting only of those vertices containing at least one target in TARGETS_OF_INPUT_DRUGS_IN_TARGETOME. These are the subtrees of each tree targeted by the collection of drugs in INPUT_DRUGS_IN_TARGETOME. Note that entries of U in the table indicate an undefined value.

Targetability --drug_coverage --target_coverage (--max_assay_value) (--filename)
Example: targetability --drug_coverage "Imatinib Mesylate" Vemurafenib NODRUG --target_coverage BHLHe39 EGFR NOTARGET

The drug_coverage command with a list of drugs and a list of targets generates the drug coverage metric data for the drug combination Imatinib Mesylate, Vemurafenib, and NODRUG targeting the target combination BHLHe39 EGFR NOTARGET. Note that Imatinib Mesylate must be entered using quotation marks to be read as a the name of a single drug, otherwise the code will interpret Imatinib and Mesylate as being the names of two distinct drugs. NODRUG has been included to illustrate how the software handles a drug name that does not appear in the Cancer Targetome. NOTARGET has been included to illustrate how the software handles a target name that does not appear in Uniprot. No values have been specified for the optional command line arguments --max_assay_value and --filename. This makes the max assay value take on the default value of 100, and the output files are written using the default filenames drug_target_coverage_YYYYmmdd_HHMMSS.csv and drug_target_coverage_YYYYmmdd_HHMMSS.txt.

Running the command generates the following console output:

--max_assay_value not specified --filename not specified Targetome query failed to find drug name NODRUG Uniprot request failed to find gene name NOTARGET Fetched: 1 / 1 Fetched: 12 / 12 Fetched: 1 / 1

The output indicates that the line arguments --max_assay_value, and --filename have not been specified. The absence of values for the --max_assay_value and --filename line arguments means that the default values for these will be used. The output also indicates that the drug name NODRUG could not be found in the Cancer Targetome and that the target name NOTARGET could not be found in Uniprot. Finally the output indicates that three queries to Uniprot fetched 1, 12, and 1 items, respectively.

The two files drug_target_coverage_YYYYmmdd_HHMMSS.csv and drug_target_coverage_YYYYmmdd_HHMMSS.txt should now appear in the directory from which the command was run. The two files contain the same information. In drug_target_coverage_YYYYmmdd_HHMMSS.txt the data are formatted to make them human readable, while drug_target_coverage_YYYYmmdd_HHMMSS.csv presents the data in a comma-separated value format so that they may be easily read into another software application.

The generated drug_target_coverage_YYYYmmdd_HHMMSS.txt contains the following information:

example/drug_target_coverage_20240502_134613.csv

example/drug_target_coverage_20240502_134613.txt

These data are interpreted as follows:

INPUT_DRUGS lists the drugs input using line arguments.

INPUT_TARGETS lists the targets input using line arguments.

INPUT_MAX_ASSAY_VALUE gives the max assay value being used. In this case it is the default value of 100.

INPUT_DRUGS_IN_TARGETOME lists the drugs in INPUT_DRUGS that were found to appear in the Cancer Targetome. In this case since NODRUG does not appear in the Cancer Targetome it does not appear here.

INPUT_TARGETS_IN_UNIPROT lists the targets in INPUT_TARGETS that where found to appear in Uniprot. In this case since NOTARGET does not appear in Uniprot it does not appear here. In addition the target BHLHe39 appears here under the name MYC. This is because Uniprot uses MYC as the preferred name for the target also designated by BHLHe39.

INPUT_DRUGS_TARGETING_INPUT_TARGETS lists the drugs in INPUT_DRUGS_IN_TARGETOME that target any of the targets in INPUT_TARGETS_IN_UNIPROT with an assay_value indicated in the Cancer Targetome to be less than or equal to INPUT_MAX_ASSAY_VALUE.

INPUT_TARGETS_TARGETED_BY_INPUT_DRUGS lists the targets in INPUT_TARGETS_IN_UNIPROT that are targeted by any of the drugs in INPUT_DRUGS_IN_TARGETOME with an assay_value indicated in the Cancer Targetome to be less than or equal to INPUT_MAX_ASSAY_VALUE.

TARGETS_OF_INPUT_DRUGS_IN_TARGETOME lists all targets of the collection of drugs in INPUT_DRUGS_IN_TARGETOME with an assay value indictaed in the Cancer Targetome to be less than or equal to INPUT_MAX_ASSAY_VALUE for at least one drug in INPUT_DRUGS_IN_TARGETOME.

ADDITIONAL_TARGETS_NOT_IN_INPUT_TARGETS lists any targets not in INPUT_TARGETS_IN_UNIPROT that are targeted with an assay value indictaed in the Cancer Targetome to be less than or equal to INPUT_MAX_ASSAY_VALUE for at least one drug in INPUT_DRUGS_IN_TARGETOME.

INPUT_TARGETS_TARGETED_BY_INPUT_DRUGS_UNIPROT_ACCESSION_IDS lists the Uniprot Accessions Ids for all targets in INPUT_TARGETS_TARGETED_BY_INPUT_DRUGS. The entries in INPUT_TARGETS_TARGETED_BY_INPUT_DRUGS and INPUT_TARGETS_TARGETED_BY_INPUT_DRUGS_UNIPROT_ACCESSION_IDS should appear in the same order.

ADDITIONAL_TARGETS_NOT_IN_INPUT_TARGETS_UNIPROT_ACCESSION_IDS lists the Uniprot Accessions Ids for all targets in ADDITIONAL_TARGETS_NOT_IN_INPUT_TARGETS. The entries in ADDITIONAL_TARGETS_NOT_IN_INPUT_TARGETS and ADDITIONAL_TARGETS_NOT_IN_INPUT_TARGETS_UNIPROT_ACCESSION_IDS should appear in the same order.

DRUG_TARGET_COVERAGE provides a table of the values of the metrics P, |L|, |Pi|, eta and KP, K|L|, k|L|, K|Pi|, k|Pi| eta measuring the coverage of the collection of drugs in INPUT_DRUGS_IN_TARGETOME across the trees extracted from Reactome with information about how much of the coverage affects the collection of targets in INPUT_TARGETS_IN_UNIPROT. The first four metrics P, |L|, |Pi|, eta provide summary information about the structures of the trees. Four metrics of the next six metrics KP, K|L|, K|Pi|, eta provide summary information on the subtree of each tree consisting only of those vertices containing at least one target in TARGETS_OF_INPUT_DRUGS_IN_TARGETOME. These are the subtrees of each tree targeted by the collection of drugs in INPUT_DRUGS_IN_TARGETOME. The remaining two metrics k|L| and k|Pi| provide summary information about what portion of the subtrees involve the targets in INPUT_TARGETS_IN_UNIPROT. Note that entries of U in the table indicate an undefined value.

Targetability --drug_tree (--max_assay_value) (--filename)
Example: targetability --drug_tree "Imatinib Mesylate" Vemurafenib NODRUG

The drug_tree command with a list of drugs generates the drug tree metric data for the drug combination Imatinib Mesylate, Vemurafenib, and NODRUG. Note that Imatinib Mesylate must be entered using quotation marks to be read as a the name of a single drug, otherwise the code will interpret Imatinib and Mesylate as being the names of two distinct drugs. NODRUG has been included to illustrate how the software handles a drug name that does not appear in the Cancer Targetome. No values have been specified for the optional command line arguments --max_assay_value and --filename. This makes the max assay value take on the default value of 100, and the output files are written using the default filenames drug_tree_YYYYmmdd_HHMMSS.csv and drug_tree_YYYYmmdd_HHMMSS.txt.

Running the command generates the following console output:

--max_assay_value not specified --filename not specified Targetome query failed to find drug name NODRUG Fetched: 13 / 13

The output indicates that the line arguments --max_assay_value, and --filename have not been specified. The absence of values for the --max_assay_value and --filename line arguments means that the default values for these will be used. The output also indicates that the drug name NODRUG could not be found in the Cancer Targetome. Finally the output indicates that a query to Uniprot fetched 13 items.

The two files drug_tree_YYYYmmdd_HHMMSS.csv and drug_tree_YYYYmmdd_HHMMSS.txt should now appear in the directory from which the command was run. The two files contain the same information. In drug_tree_YYYYmmdd_HHMMSS.txt the data are formatted to make them human readable, while drug_tree_YYYYmmdd_HHMMSS.csv presents the data in a comma-separated value format so that they may be easily read into another software application.

The generated drug_tree_YYYYmmdd_HHMMSS.txt contains the following information:

example/drug_tree_20240502_134649.csv

example/drug_tree_20240502_134649.txt

These data are interpreted as follows:

INPUT_DRUGS lists the drugs input using line arguments.

INPUT_MAX_ASSAY_VALUE gives the max assay value being used. In this case it is the default value of 100.

INPUT_DRUGS_IN_TARGETOME lists the drugs in INPUT_DRUGS that were found to appear in the Cancer Targetome. In this case since NODRUG does not appear in the Cancer Targetome it does not appear here.

TARGETS_OF_INPUT_DRUGS_IN_TARGETOME lists all targets of the collection of drugs in INPUT_DRUGS_IN_TARGETOME with an assay value indictaed in the Cancer Targetome to be less than or equal to INPUT_MAX_ASSAY_VALUE for at least one drug in INPUT_DRUGS_IN_TARGETOME.

TARGETS_OF_INPUT_DRUGS_IN_TARGETOME_UNIPROT_ACCESSION_IDS lists the Uniprot Accessions Ids for all targets in TARGETS_OF_INPUT_DRUGS_IN_TARGETOME. The entries in TARGETS_OF_INPUT_DRUGS_IN_TARGETOME and TARGETS_OF_INPUT_DRUGS_IN_TARGETOME_UNIPROT_ACCESSION_IDS should appear in the same order.

DRUG_TREE provides a table of the values of the metrics P, |L|, |Pi|, eta and |L|, |Pi|, eta measuring the coverage of the collection of drugs in INPUT_DRUGS_IN_TARGETOME across the trees extracted from Reactome. The first four metrics P, |L|, |Pi|, eta provide summary information about the structures of the trees. The next three metrics |L|, |Pi|, eta provide summary information on the subtree of each tree consisting only of those vertices containing at least one target in TARGETS_OF_INPUT_DRUGS_IN_TARGETOME. These are the subtrees of each tree targeted by the collection of drugs in INPUT_DRUGS_IN_TARGETOME. Note that entries of U in the table indicate an undefined value.

Targetability --target_tree (--max_assay_value) (--filename)
Example: targetability --target_tree BHLHe39 EGFR NOTARGET

The target_tree command with a list of targets generates the target tree metric data for the target combination BHLHe39, EGFR, and NOTARGET. NOTARGET has been included to illustrate how the software handles a target name that does not appear in Uniprot. No values have been specified for the optional command line arguments --max_assay_value and --filename. This makes the max assay value take on the default value of 100, and the output files are written using the default filenames target_tree_YYYYmmdd_HHMMSS.csv and target_tree_YYYYmmdd_HHMMSS.txt.

Running the command generates the following console output:

--max_assay_value not specified --filename not specified Uniprot request failed to find gene name NOTARGET Fetched: 2 / 2

The output indicates that the line arguments --max_assay_value, and --filename have not been specified. The absence of values for the --max_assay_value and --filename line arguments means that the default values for these will be used. The output also indicates that the target name NOTARGET could not be found in Uniprot. Finally the output indicates that a query to Uniprot fetched 2 items.

The two files target_tree_YYYYmmdd_HHMMSS.csv and target_tree_YYYYmmdd_HHMMSS.txt should now appear in the directory from which the command was run. The two files contain the same information. In target_tree_YYYYmmdd_HHMMSS.txt the data are formatted to make them human readable, while target_tree_YYYYmmdd_HHMMSS.csv presents the data in a comma-separated value format so that they may be easily read into another software application.

The generated target_tree_YYYYmmdd_HHMMSS.txt contains the following information:

example/target_tree_20240502_134716.csv

example/target_tree_20240502_134716.txt

These data are interpreted as follows:

INPUT_TARGETS lists the targets input using line arguments.

INPUT_MAX_ASSAY_VALUE gives the max assay value being used. In this case it is the default value of 100.

INPUT_TARGETS_IN_UNIPROT lists the targets in INPUT_TARGETS that where found to appear in Uniprot. In this case since NOTARGET does not appear in Uniprot it does not appear here. In addition the target BHLHe39 appears here under the name MYC. This is because Uniprot uses MYC as the preferred name for the target also designated by BHLHe39.

INPUT_TARGETS_IN_UNIPROT_UNIPROT_ACCESSION_IDS lists the Uniprot Accessions Ids for all targets in INPUT_TARGETS_IN_UNIPROT. The entries in INPUT_TARGETS_IN_UNIPROT and INPUT_TARGETS_IN_UNIPROT_UNIPROT_ACCESSION_IDS should appear in the same order.

TARGET_TREE provides a table of the values of the metrics P, |L|, |Pi|, eta and |L|, |Pi|, eta measuring the coverage of the collection of targets in INPUT_TARGETS_IN_UNIPROT across the trees extracted from Reactome. The first four metrics P, |L|, |Pi|, eta provide summary information about the structures of the trees. The next three metrics |L|, |Pi|, eta provide summary information on the subtree of each tree consisting only of those vertices containing at least one target in INPUT_TARGETS_IN_UNIPROT. Note that entries of U in the table indicate an undefined value.

Optional Line Arguments
Targetability --filename
Example: targetability --drug_coverage "Imatinib Mesylate" Vemurafenib --filename FILENAME

The output files are written using the filenames FILENAME.csv and FILENAME.txt.

Targetability --max_assay_value
Example: targetability --drug_coverage "Imatinib Mesylate" Vemurafenib --max_assay_value MAX_ASSAY_VALUE

The max assay value take on the value of MAX_ASSAY_VALUE. The software will only return results for drug-target interactions with an assay value in the Cancer Targetome less than or equal to the max assay value.
