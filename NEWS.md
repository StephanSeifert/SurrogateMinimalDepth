# SurrogateMinimalDepth 0.1.10
        - fixed bug not including first and last surrogate
	      - added random selection of surrogates when adj_agree are the same
	      - fixed bug in var.relation example
	      - added build.clusters function to obtain variable groups based on mean adjusted agreement 

# SurrogateMinimalDepth 0.1.9	
        - some specifics to set s in var.select.smd were adapted
		    - added save.memory parameter (to build the forest with ranger) to var.select.smd and var.select.md  
		    
# SurrogateMinimalDepth 0.1.8	
        - MD now executes var.select.smd with s = 0 instead of using a separate function
		    - The function reduce.surrogates is included. 
		    - Included parameters create.forest and forest to use var.select.smd for existing forests (e.g. created by reduce.surrogates)
		    - Changed the value "trees" to "forest" containing trees and variable names 
		    - The C-code is updated to enable multicore analysis

# SurrogateMinimalDepth 0.1.7	
        - Included parameter save.ranger in var.select.md and var.select.smd to save ranger object

# SurrogateMinimalDepth 0.1.6	
        - Adapted threshold for low depth trees in var.select.smd

# SurrogateMinimalDepth 0.1.5 	
        - Added s as value in var.select.smd 
		    - Implemented survival function for var.select.smd and var.select.md

# SurrogateMinimalDepth 0.1.4 	
        - All errors and comments from coauthors implemented: first version uploaded on github
