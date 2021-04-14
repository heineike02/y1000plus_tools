# y1000plus_tools

This code provides the scripts used to analyse protein and promoter conservation in budding yeast for [1] based on data from [2]. 

## Getting Started

The main functions are contained in y1000_plus_tools.py

Scripts to generate analysis and figures are in the scripts folder
The scripts used for ref [1] are listed in PKA_evolution_fig_notes.txt
Also in the scripts folder is notebook_setup.py which you can load at the beginning of a script with: 

%load notebook_setup.py

It runs the std_libraries.py file which loads all necessary libraries for the provided scripts. 

Some data is the folder y1000plus_tools_data which can be obtained from figshare.  

The following from [2] can be found at https://figshare.com/articles/dataset/Tempo_and_mode_of_genome_evolution_in_the_budding_yeast_subphylum/5854692

Folders:

0_332yeast_genomes
orthomcl_output 

Files:

332_2408OGs_time-calibrated_phylogeny_species-names_updated.newick
343taxa_2408OGs_long_seqIDs.txt
343taxa_proteins.fasta
343taxa_speicies-name_clade-name_color-code


## References

[1] Heineike, B., and El-Samad, H. (2021). Paralogs in the PKA regulon traveled different evolutionary routes to divergent expression in budding yeast. Front. Fungal Biol. 2. doi:10.3389/ffunb.2021.642336.

[2] Shen, X.-X., Opulente, D. A., Kominek, J., Zhou, X., Steenwyk, J. L., Buh, K. V., et al. (2018). Tempo and Mode of Genome Evolution in the Budding Yeast Subphylum. Cell 175, 1533-1545.e20. doi:10.1016/j.cell.2018.10.023.


## Authors

* **Benjamin Heineike** [heineike02](https://github.com/heineike02)


## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
