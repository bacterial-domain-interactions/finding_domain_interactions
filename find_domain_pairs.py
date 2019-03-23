from Domain_Interactions import get_my_domain_combinations
from Sort_Uniprot import give_me_my_domain_list
from Combining_paired_domains import get_me_my_paired_domain_list

name = input("What is the name of the organism you are looking at? Example answers: e.coli, strep?")

paired_list = input("Please enter your paired file: ")
bait_list = input("Please input your list of bait proteins with their domains: ")
prey_list = input("Please input your list of prey proteins with their domains: ")
domain_id = input("Are you using PFAM ids or INTERPRO ids? If using PFAM, type PFAM. If using INTERPRO, type INTERPRO: ")
export_file_name = name + "_" + domain_id + "_"

print()
print("Your domain frequency list is being computed...")

bait_domains = give_me_my_domain_list(paired_list, "protein 1", bait_list, domain_id)
prey_domains = give_me_my_domain_list(paired_list, "protein 2", prey_list, domain_id)
paired_file = get_me_my_paired_domain_list(paired_list, bait_domains, prey_domains)
get_my_domain_combinations(paired_file, export_file_name)

print()
print("Your domain frequency list has been computed! Your results are in csv file ending with Xresults or ending in Yresults")