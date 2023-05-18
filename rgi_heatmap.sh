#!/bin/bash

# Yellow represents a perfect hit, teal represents a strict hit, purple represents no hit.

rgi heatmap \
    --input rgi_heatmap \
    --cluster genes \
    --output rgi_heatmap/RGI_heatmap



# NOTE: 
#   with --category drug_class
#       AMR genes categorised by Drug Class. Yellow represents a perfect hit, teal represents a strict hit, purple represents no hit. 
#       Genes with asterisks (*) appear multiple times because they belong to more than 
#       one Drug Class category in the antibiotic resistance ontology (ARO).


#   with --category resistance_mechanism
#       AMR genes categorised by Resistance Mechanism. Yellow represents a perfect hit, teal represents a strict hit, 
#       purple represents no hit. Genes with asterisks (*) appear multiple times because they belong to more than
#       one Resistance Mechanism category in the antibiotic resistance ontology (ARO)

#   with --category gene_family
#       AMR genes categorised by AMR Gene Family. Yellow represents a perfect hit, teal represents a strict hit, 
#       purple represents no hit. Genes with asterisks (*) appear multiple times because they belong to more than 
#       one AMR Gene Family category in the antibiotic resistance ontology (ARO).
