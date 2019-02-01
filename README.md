# finding_domain_interactions
Given a list of paired protein interactions, this program will identify all possible domain combinations and reports the frequency of each domain combination.

In more detail, this program will:

About this file:
Input file: (this example input file uses PFAM identifiers, but this program can also use InterPro ids)

Protein X    Xdomains                           Protein Y           YDomains
P10408	     PF02810 PF07517 PF01043 PF07516 	PF31119	            PF01553 PF00501 
P0AA37	     PF00849 	                        PF0A8P1	            PF03588 
P14294	     PF01131 PF01751 	                PF46133	            PF03806 



Each row contains a pair of interacting proteins (Protein X and Protein Y) and the corresponding domains for 
proteinX and protein Y

This program aims to:
     1. Make a nested list of XDomains and YDomains
     2. Find all possible domain combinations for an interacting pair
          Ex: for the third interaction (P14294 and P46133)
                Possible domain interactions are:
                    PF01131 - PF46133
                    PF01751 - PF46133
     3. Remove duplicate combinations from a list inside a nested list
                Ex: PF01131 - PF46133
                    PF46133  - PF01131  (Will remove this combination pair)
     4. Count the number of times each domain combination appears
            Ex: PF01131 - PF46133  4
                PF01751 - PF46133  1
     5. Add reverse combination counts, for example:
            Ex: PF01131 - PF46133  4
                PF46133 - PF01131  3
            
            So you add the reverse combination to the original count so:
                Ex: PF01131 - PF46133  7
            
     6. Shuffle one of domain lists (Ex: XDomains)
     7. Then do processes 2-4 again
     8. Do steps 5-6 a 1000 times, while keeping a dictionary that looks like this:
                        {PF01131 - PF46133: [245250, 9760]}
                        245250 - the sum of the total number of times a combination appears in a randomization
                        9760 - the total number of times the combination appeared in the 1000 randomizations

     9. Take the average number of times a combination appeared (Ex: 245250/9760 = 25)
     10. Do steps 5-8 for the other domain lists (YDomains)
     11. For each domain pair, get the number of times the combination appeared for the non-shuffled combinations (
         (step 4) and the number of times the combination appeared for the shuffled combinations (step 8) and take
         log(non-shuffled/shuffled).
     12. Export results
