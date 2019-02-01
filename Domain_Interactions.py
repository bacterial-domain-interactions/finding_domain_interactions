# Nikhita Puthuveetil
import pandas as pd
from collections import Counter
import itertools
from random import sample
from math import log10
import copy

""" About this file:
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
"""

# Functions ------------------------------------------------------------------------------------------


def domain_combinations(proteinX, proteinY):
    """
     Find all possible domain combinations for one interacting pair
    :param proteinX: a list of domains from one protein in a pair of interacting proteins Ex: ["A", "B"]
    :param proteinY: a list of domains from the other protein in a pair of interacting proteins Ex: ["C"]
    :return: returns a list of all possible domain combinations for a pair of interacting proteins
             Ex: ["A-C", "B-C"]
    """
    domain_combo = []

    for domains1 in proteinX:
        for domains2 in proteinY:
            combo = domains1 + "-" + domains2
            domain_combo.append(combo)
    return domain_combo


# function to find the total domain combinations
def total_combinations(domainslist1, domainslist2):
    """
    Find all possible domain combinations for multiple interactions
    :param domainslist1: nested list of domains for a column of interactant in a pair of interacting proteins
           Ex: XDomains in input file in header comment
    :param domainslist2: nested list of domains for a column of the other interactant in a pair of interacting proteins
           Ex: YDomains in input file in header comment
    :return: returns a nested list of all possible domain combinations for interacting proteins
    """
    total_combo = []
    for i in range(len(domainslist1)):
        combination = domain_combinations(domainslist1[i], domainslist2[i])
        total_combo.append(combination)
    return total_combo


def remove_reverse_duplicates(total_combinations_list):
    """
    Calculating domain combinations can yield duplicates, such as when calculating combinations for domains
    ["A", "B"] and ["A"] will yield ["A-B", "A-A", "B-A", "B-B"], this method will remove "B-A"
    :param total_combinations_list: returns a nested list of all possible domain combinations for interacting proteins
    :return: returns a nested list of all possible domain combinations for a pair of interacting proteins with
             duplicates removed
    """
    for combinations in total_combinations_list:
        for individ_combos in combinations:
            # if dealing with PFAM ids
            if individ_combos.startswith("PF"):
                reverse_PFAM = individ_combos[8:] + "-" + individ_combos[0:7]
                # if reverse is in the list and to ignore combinations like PF02134-PF02134
                if reverse_PFAM in combinations and individ_combos[0:7] != individ_combos[8:]:
                    combinations.remove(individ_combos)
            # if dealing with interPro ids
            elif individ_combos.startswith("IPR"):
                reverse_interPro = individ_combos[10:] + "-" + individ_combos[0:9]
                if reverse_interPro in combinations and individ_combos[0:9] != individ_combos[10:]:
                    combinations.remove(individ_combos)

    return total_combinations_list


def combo_counter(combo_list):
    """
    Will count the number of times each combinations appears in combo_list
    :param combo_list: returns a nested list of all possible domain combinations for interacting proteins
    :return: a dictionary of the counts for each of the combinations in combo_list
    """
    # combine lists of lists to one list
    combined_list = itertools.chain.from_iterable(combo_list)
    combined_list = list(combined_list)
    count = Counter(combined_list)
    return count


def add_reverse_combos(count_dict):
    """
    Once a dictionary of counts for each domain combination has been made, there will be duplicates. Take for example, the domain combination
    PF00486-PF00072 which has a count of 13. Its reverse combination, PF00072-PF00486 has a count of 7. Since the order of the domains do not
    matter, we treat both of these combinations as the same so need to add their counts, so PF00486-PF00072 will have a count of 20.
    This method aims to go through the dictionary of counts, see if the reverse combination has a count of more than 0, if it does,
    then it will add it to the original combination and delete the reverse combination.
    Additionally, the method will remove any reverse combinations that has a count of 0.
    :param count_dict: dictionary of counts
    :return: edit_count: dictionary of counts with the reverse combos added/deleted
    """
    edit_count = copy.deepcopy(count_dict)
    for keys in count_dict.keys():
        # if keys in dict are PFAM ids
        if keys.startswith("PF"):
            # if the reverse combination has counts and you want to skip over entries like this PF00072-PF00072
            if count_dict[keys[8::] + "-" + keys[0:7]] > 0 and count_dict[keys] != count_dict[keys[8::] + "-" + keys[0:7]]:
                # if the reverse combination has counts, then add it to original combination and delete the reverse combination
                edit_count[keys] = count_dict[keys] + count_dict[keys[8::] + "-" + keys[0:7]]
                del edit_count[keys[8::] + "-" + keys[0:7]]

            # if reverse combo has no counts, then delete it
            if count_dict[keys[8::] + "-" + keys[0:7]] == 0:
                del edit_count[keys[8::] + "-" + keys[0:7]]

        # if keys in dict are InterPro ids
        if keys.startswith("IPR"):
            if count_dict[keys[10::] + "-" + keys[0:9]] > 0 and count_dict[keys] != count_dict[keys[10::] + "-" + keys[0:9]]:
                edit_count[keys] = count_dict[keys] + count_dict[keys[10::] + "-" + keys[0:9]]
                del edit_count[keys[10::] + "-" + keys[0:9]]

            if count_dict[keys[10::] + "-" + keys[0:9]] == 0:
                del edit_count[keys[10::] + "-" + keys[0:9]]

    return edit_count


def randomize_a_1000(domain_list1, domain_list2):
    """
    Randomization (see steps 5-8)
    :param domain_list1: nested list of domains for a column of interactant in a pair of interacting proteins (list you want to shuffle)
    :param domain_list2: nested list of domains for a column of the other interactant in a pair of interacting proteins (list remains unchanged)
    :return: a dictionary of the average counts for each domain combination
    """
    my_count = {}
    for times in range(1000):
        randomized_domain = sample(domain_list1, len(domain_list1))
        randomized_domain_combo = remove_reverse_duplicates(total_combinations(randomized_domain, domain_list2))
        randomized_count = add_reverse_combos(combo_counter(randomized_domain_combo))

        for k, v in randomized_count.items():
            # if combination is new (not already in the dictionary), make a new key for it
            if k not in my_count.keys():
                my_count[k] = [v, 1]

            else:
                my_count[k][0] += v  # the sum of all the counts of a particular combinations
                my_count[k][1] += 1  # the number of times the combination appears

    # Take the average of the sum of the counts and the number of times the combination appears
    for k, v in my_count.items():
        # my_count[k][0] = my_count[k][0]/my_count[k][1]
        my_count[k][0] = my_count[k][0] / 1000
        my_count[k] = my_count[k][0]

    return my_count


def consolidate_my_information(original_dict, shuffled_dict):
    """
    Now that the counts for the domain pairs have been calculated (and the randomization has been done), this function
    will consolidate the information into one dictionary and then normalize the counts and then see if a given domain
    pair can be found in the 3did database
    :param original_dict: dictionary of counts with the reverse domain combos added/deleted
    :param shuffled_dict: a dictionary of the average counts for each domain combination
    :return:
    """
    # Open the data from 3did database
    all_3did = []
    with open("3did_validation.dat", "r") as fh:
        for lines in fh:
            if lines.startswith("#=ID"):
                lines = lines.replace("@Pfam", "")
                start = lines.find("(")
                all_3did.append(lines[start:].replace("(", "").replace(")", "").split())

    complete_dict = {}
    log_keys = 0
    my_3did = ""
    for keys in shuffled_dict.keys():
        if keys in original_dict.keys():

            # Take the log values
            if shuffled_dict[keys] != 0:
                division = original_dict[keys] / shuffled_dict[keys]
                if division == 0.0:
                    log_keys = 0
                else:
                    log_keys = log10(division)

            # Find in 3did
            found = False

            # If domains are labeled with PFAM IDS
            if keys.startswith("PF"):
                domain1 = keys[0:7]
                domain2 = keys[8:]

            # If domains are labeled with InterPro IDs
            elif keys.startswith("IPR"):
                domain1 = keys[0:9]
                domain2 = keys[10:]

            for data in all_3did:
                if domain1 in data[0] and domain2 in data[1]:
                    my_3did = "Yes"
                    found = True

                elif domain1 in data[1] and domain2 in data[0]:
                    my_3did = "Yes"
                    found = True

            if found is False:
                my_3did = "No"

            complete_dict[keys] = [original_dict[keys], shuffled_dict[keys], log_keys, my_3did]

    return complete_dict


def export_combos(domain, complete_dict, species_name):
    """
    Exports counts from the complete dictionary of values (which includes the domain counts from the original_dict,
    shuffled_dict, the log_values, and if the domain pair is in 3did
    :param domain: a String indicating which domain is shuffled in shuffled_dict (ex: X)
    :param complete_dict: a complete dictionary of values containing counts, log values, and if the pair is in 3did
    :param name of the species the proteins are from
    :return: a csv file with the domain combination and the corresponding counts in from the nonshuffled dictionary,
             shuffled dictionary, the log values, and if the pair is in the 3did database
    """
    outputfile = species_name + "_" + domain + "results.csv"
    with open(outputfile, "w") as writer:
        header = "Combination" + "," + "Original Count" + "," + "Shuffled " + domain + "Count" + \
                 "," + "Log(original/shuffled" + domain + ")" + "," + "Domain in 3did" + "\n"
        writer.write(header)
        for entry in complete_dict.keys():
            row = entry + "," + str(complete_dict[entry][0]) + "," + str(complete_dict[entry][1]) + "," \
                    + str(complete_dict[entry][2]) + "," + str(complete_dict[entry][3]) + "\n"
            writer.write(row)


def get_my_domain_combinations(file_name, species_name):
    """

    :param file_name: a csv file containing a paired protein list with their domains (see info at the very top of this program)
    :param species_name: name of the species the proteins are coming from
    :return: a csv file with the domain combination and the corresponding counts in from the nonshuffled dictionary,
             shuffled dictionary, the log values, and if the pair is in the 3did database
    """
    df = pd.read_csv(file_name)

    """
    Make a nested list of XDomains and YDomains
    In the XDomains column, take each row and domains for each protein into a list
    X_domain will be a nested list (each list containing the domains for a protein) (Same will be done for YDomains)
    """
    X_domain = []
    for i in df.XDomains:
        X_split = i.split()
        X_domain.append(X_split)

    Y_domain = []
    for i in df.YDomains:
        Y_split = i.split()
        Y_domain.append(Y_split)

    """Find domain combinations for lists of domains for X and Y (Steps 2-4)"""
    all_combos = remove_reverse_duplicates(total_combinations(X_domain, Y_domain))
    original_count = combo_counter(all_combos)
    edited_original_count = add_reverse_combos(original_count)

    """Shuffle the first domain list"""
    shuffled_X_count = randomize_a_1000(X_domain, Y_domain)

    """Shuffle the second domain list"""
    shuffled_Y_count = randomize_a_1000(Y_domain, X_domain)

    """Export results"""
    export_combos("X", consolidate_my_information(edited_original_count, shuffled_X_count), species_name)
    export_combos("Y", consolidate_my_information(edited_original_count, shuffled_Y_count), species_name)
