## **Important**
This project is still in progress and will be continually updated.

# Genome scale metabolic model of *Xanthobacter sp. SoF1*

The novel food, which is subject to the present novel food application, and has been submitted in accordance with Regulation (EU) 2015/2283, is a microbial protein-rich powder comprised of biomass of inactivated bacterial cells of a non-genetically modified strain of *Xanthobacter sp.*
The production of the novel food is based on the cultivation of *Xanthobacter sp. SoF1* in a continuous bioprocess. It is well characterized and carefully monitored. Quality assurance programs have been established to ensure the product meets HACCP (Hazard Analysis and Critical Control Points) principles
and release specifications. The composition of the novel food has been characterized and product specifications have been established. Detailed physical, chemical, and microbiological product analysis, as well as stability tests show
the safety and long-term stability of the novel food.
The novel food is intended to be used as an ingredient in different food groups. The novel food is not intended to replace another food but may partially replace the consumption of meat in the non-vegetarian population and the consumption of meat imitates in the vegetarian population. 

## **Updates** (14/06/2024)
- The growth rate parameter was recalibrated for accuracy.
- The complete Calvin cycle was integrated into the model.
- Additional significant metabolic pathways were incorporated.

# Annotation 
Genome annotation serves as a critical initial step in deciphering the genetic blueprint of organisms, providing fundamental insights into gene function and metabolic potential. Our approach to annotating the Xanthobacter sp. SoF1 genome involved a multi-faceted strategy integrating different bioinformatic tools and databases. 
Annotation tools utilize sequence alignment algorithms to compare genetic sequences with reference databases, predict genes based on features such as open reading frames and splice sites, and assign putative functions to proteins through homology searches and domain predictions. Initially, we used Bakta (Schwengers et al., 2021), a widely used annotation tool, to identify 4520 genes in the genome. Recognizing the complexity of genomic functions, we extended the scope of our analysis to include data from EggNog (Cantalapiedra et al., 2021), RAST (Overbeek et al., 2014) and Prokka (Seemann, 2014).
The additional insights gained from these databases significantly improved the scope and accuracy of our annotation. Specifically, we identified 113 previously unrecognized gene functions. To further enrich our annotation, we used a machine learning tool called DeepEC Transformer. This allowed predicting 824 new EC numbers, critical for elucidating enzyme functions and pathways. By integrating the data of DeepEC Transformer, we strengthened our annotation and ensured a more nuanced representation of the metabolic potential of Xanthobacter sp. SoF1.

# Draft Metabolic reconstruction  

We drafted a metabolic reconstruction from the annotated Xanthobacter sp. SoF1 genome using CarveMe (Machado et al., 2018), leveraging its top-down reconstruction approach. This methodology commences with a manually curated universal metabolic model, addressing common issues encountered in reconstructions such as missing reactions and incomplete pathways. Unlike traditional bottom-up methods, CarveMe infers an organism's metabolic capabilities solely from genetic evidence, enhancing the efficiency and accuracy of our draft model. Through the 'carving' process, this universal model is tailored to specific organisms, ensuring the preservation of structural properties while streamlining reconstruction. The draft model contained 2258 reactions, 1534 metabolites and 941 metabolic genes.

# Manual Curation 

The draft reconstruction obtained a medium MEMOTE score of 67% (Lieven et al., 2020), indicating the need for further refinement. This assessment tool evaluates various parameters, including annotation tests, general quality metrics, and stoichiometric consistencies, in accordance with community standards . The entire model underwent rigorous charge and mass balance curation procedures, resulting in a total of 431 reactions being successfully balanced. This process was accompanied by the removal of exchange reactions that induced unrealistic growth rate. The biomass composition underwent curation based on available experimental data, leading to a significant enhancement in model quality, evidenced by a MEMOTE score of 88%. Nevertheless, further refinement is necessary to incorporate missing pathways in the model.

