# clust-assess
Assessing pling clusters of potentially related plasmids

# Objective

Given a set of plasmids that are suspected on being involved in an outbreak, are these plasmids more closely related to each other (lower intra-cluster variability) than others in public databases based on the pling subcommunity?

## Inputs

A multi-fasta file containing the plasmids of interest, and a database (pling PLSDB)

## Outputs

A results file listing the within-subcommunity variability of each input plasmids relative to the within-subcommunity variability of closest databases plasmids. 

If the within-subcommunity variability of the input plasmids is much lower than within-cluster variability, this indicates the plasmids may be closely related and part of an outbreak. If the within-subcommunity variability of the input plasmids is much higher than within-cluster variability, this indicates that the plasmids are not closely related. If the subcommunities for all input plasmids are different, then the plasmids are not related. 

# Tasks 
- [ ] get dataset from galaxy onto climb - Nathan
- [ ] PLSDB pling onto CLIMB 
- [ ] run ACCIO on Galaxy test dataset
- [ ] split PLSDB file prior to ACCIO
    - [ ] like pling sourmash pre-filtering, separate it by community or subset PLSDB regionally 
- [ ] testing pling add with `-asyn` to not re-compute clusters and run fast 
- [X] make GitHub page 
- [ ] writing workflow - Robin
- [ ] 
